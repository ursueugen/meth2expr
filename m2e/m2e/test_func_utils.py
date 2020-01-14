import unittest
import os

from func_utils import *

DIR = '../data/'
GENOMICS_DIR = DIR + 'genomics/'

class TestDataCollection(unittest.TestCase):


    def test_get_gff(self):
        
        url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz"
        dir = GENOMICS_DIR

        path = get_gff(url, dir)
        self.assertTrue(os.path.isfile(path))


    def test_get_genome(self):

        url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz"
        dir = GENOMICS_DIR

        path = get_gff(url, dir)
        self.assertTrue(os.path.isfile(path))


    def test_get_cpgs(self):

        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL13534&id=11288&db=GeoDb_blob92"
        dir = GENOMICS_DIR

        path = get_cpgs(url, dir)
        self.assertTrue(os.path.isfile(path))


if __name__ == '__main__':
    unittest.main()