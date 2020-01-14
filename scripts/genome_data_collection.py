import os
import logging

from pybedtools import BedTool
from pybedtools.featurefuncs import five_prime
from pybedtools.cbedtools import create_interval_from_list

from m2e.func_utils import get_cpgs, get_genome, get_gff
from m2e.config import configs, CHROMOSOMES


### inputs
GENOME_URL = configs["urls"]["genome"]
GFF_URL = configs["urls"]["gff"]
CPG_URL = configs["urls"]["cpg"]
GENOMICS_DIR = configs["dirs"]["genomics"]
GENOME_GFF_PATH = os.path.join(GENOMICS_DIR, configs["names"]["genome_complete_gff"])
GENOME_SEQ_PATH = os.path.join(GENOMICS_DIR, configs["names"]["genome_seq"])

UPSTREAM_LENGTH = configs['params']['upstream_len']

### outputs
GENES_GFF_PATH = os.path.join(GENOMICS_DIR, configs['names']['genome_genes_gff'])
PROMS_GFF_PATH = os.path.join(GENOMICS_DIR, configs['names']['proms_seq'])
PROMS_SEQ_PATH = os.path.join(GENOMICS_DIR, configs['names']['proms_seq'])

# log
logging.basicConfig(filename=os.path.join(configs['dirs']['log'], 'genomics_data_collection_logs.txt'), 
                    level=logging.INFO, 
                    format='%(asctime)s %(message)s', 
                    datefmt='%m/%d/%Y %I:%M:%S %p')


def add_geneIDs(bed):
    """ Add gene IDs to name field of gff for rendering fasta/tab sequence files with gene ids in headers."""
    new_intervals = []
    for interval in bed:
        new_fields = interval.fields[:]
        new_fields[2] = interval.name.split(":")[1] + "_promoter"
        new_interval = create_interval_from_list(new_fields)
        new_intervals.append(new_interval)
    return BedTool(new_intervals)


if __name__ == "__main__":
    

    # check availability of genome sequence and genome annotations(gff)
    if not os.path.isfile(GENOME_SEQ_PATH):
        get_genome(GENOME_URL)

    if not os.path.isfile(GENOME_GFF_PATH):
        get_gff(GFF_URL)

    if not os.path.isfile(CPG_PATH):
        get_cpgs(CPG_URL)

    # derive promoter gff and extract sequences from genome
    if os.path.isfile(PROMS_GFF_PATH):
        proms_bed = BedTool(PROMS_GFF_PATH)
        logging.info("Found proms gff at " + PROMS_GFF_PATH)
    else:
        # point to genome gff3
        genome_bed = BedTool(GENOME_GFF_PATH)
    
        # filter for genes
        genes_bed = ( genome_bed.filter(lambda x: (x[2] == 'gene') and (x.chrom in CHROMOSOMES))
                      .saveas(GENES_GFF_PATH) )
        logging.info("Extracted # genes = " + str(len(genes_bed)) + ", saved at " + GENES_GFF_PATH)
    
        # extract promoters from genes features
        proms_bed = genes_bed.each(func=five_prime, upstream=UPSTREAM_LENGTH, downstream=0)
        proms_bed = add_geneIDs(proms_bed)
        proms_bed.saveas(PROMS_GFF_PATH)
        logging.info("Extracted # promoters = " + str(len(proms_bed)) + ", saved at " + PROMS_GFF_PATH)
    
    # extract sequences defined by promoter features
    proms_seqs = proms_bed.sequence(fi=GENOME_SEQ_PATH, fo=PROMS_SEQ_PATH, s=True, fullHeader=True, name=True)
    logging.info("Extracted prom seqs. Saved at " + PROMS_SEQ_PATH)