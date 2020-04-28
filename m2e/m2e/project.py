"""
Project class for handling project data loading.
"""

import warnings

import numpy as np
from pathlib import Path
import pandas as pd


PROJECTS_DIR = "/data/eugen/tcga/projects/"


class Project(object):
    """
    Project data loader and descriptor.
    
    :attributes:
    
    :public:
        get_case_data()
    
    :private:
        _collect_samples()
        _collect_metadata()
        _get_case_datapaths()
    """
    
    def __init__(self, name, projects_dir=PROJECTS_DIR):
        '''
        :params
            name -- project name
            path -- project dir
        '''
        self.name = name
        self.dir = projects_dir + name
        
        self.meth_path = None
        self.meth_fpath = None
        
        self.samples = {}
        self.sample_ids = []
        self.sample_paths = {}
        self.cases = {}
        self.case_ids = []
        self.meta = {}
        
        # Get samples ids and paths
        self.collect_samples_()
        self.collect_metadata_()
        
        
    def collect_samples_(self):
        """
        Extracts samples' ids and paths from project directory.
        """
        
        self.sample_paths = {"methylation": {},
                              "expression": {}}
        
        # Extract methylation
        self.meth_path = Path(self.dir) / 'data/methylation' / self.name / 'harmonized/DNA_Methylation/Methylation_Beta_Value'
        self.expr_path = Path(self.dir) / 'data/expression' / self.name / 'harmonized/Transcriptome_Profiling/Gene_Expression_Quantification'
        self.sample_paths = {datatype: {str(f).split("/")[-1]: str(f) 
                                        for f in path.iterdir()}
                             for datatype, path in [('methylation', self.meth_path), ('expression', self.expr_path)]}
        
        # In path there is an experiment ID, not sample id
        self.samples = {f: list(self.sample_paths[f].keys()) for f in ['methylation', 'expression']}
        self.sample_ids = list(self.samples.keys())
    
    
    def collect_metadata_(self):
        """
        Gets all necessary metadata of project.
        Includes ids and case ids: necessary to match expression and methylation measurements.
        """
        
        metadata_paths = [(Path(self.dir) / (self.name + "_" + f + ".csv"))
                              for f in ['methylation', 'expression'] ]
        
        ### Access metadata files
        cols = ['file_id', 'file_name', 'cases']
        meth, expr = list(map(lambda p: pd.read_csv(p, sep='\t', usecols=cols), metadata_paths))
        
        
        ### Prune the TCGA barcode of case ids up to sample level, resulting in Project-TSS-ParticipantID-<SampleID><vial>
        barcode_splitter = lambda x: "-".join(x.split('-')[:4])
        meth['cases'] = meth['cases'].apply(barcode_splitter)
        expr['cases'] = expr['cases'].apply(barcode_splitter)

                          
        meth = meth.set_index("file_id", drop=False)
        expr = expr.set_index("file_id", drop=False)
        
        meta = {'methylation': meth.to_dict(orient='index'),
                'expression': expr.to_dict(orient='index') }
        
        # Assert equivalency of ids in metadata files and available files
        assert set(meta['methylation'].keys()) == set(self.samples['methylation'])
        assert set(meta['expression'].keys()) == set(self.samples['expression'])
        
        self.meta = meta
        
        
        ### Collect metadata for cases
        
        # There may be case duplicates, will drop
        meth_dup = meth['cases'][meth.duplicated(['cases'], keep=False)]
        expr_dup = expr['cases'][expr.duplicated(['cases'], keep=False)]
        dup = set(meth_dup).union(set(expr_dup))  # set of duplicated cases in meth or expr
        if (len(dup) > 0):
            warnings.warn("Droping duplicated cases entries.", UserWarning)
        meth = meth.loc[~(meth['cases'].isin(dup))]
        expr = expr.loc[~(expr['cases'].isin(dup))]
        
        meth = meth.set_index("cases")
        expr = expr.set_index("cases")
        
        # recheck equivalency of cases for methylation and expression (should be true by constraint of download script)
        if (set(meth.index) != set(expr.index)):
            warnings.warn("There is non-equivalence of cases in methylation and expression.")
            common_cases = set(meth.index).intersection(set(expr.index))
            meth = meth.loc[common_cases]
            expr = expr.loc[common_cases]
            
        
        cases = {'methylation': meth.to_dict(orient='index'),
                 'expression': expr.to_dict(orient='index') }
        self.old_cases = cases
            
            
        # make cases as keys { case -> {methylation -> {}; expression  -> {}} }
        new_dict = {case: {'methylation': cases['methylation'][case],
                           'expression': cases['expression'][case]}
                    for case in cases['methylation'].keys() }
        
        self.cases = new_dict
        self.case_ids = list(self.cases.keys())
        
        
    def get_case_datapaths_(self, case: str):
        '''Getter method for meth, expr or both types data.'''
        
        assert case in self.cases.keys()
        
        self.meth_fpath = Path(self.meth_path) / self.cases[case]['methylation']['file_id'] / self.cases[case]['methylation']['file_name']
        self.expr_fpath = Path(self.expr_path) / self.cases[case]['expression']['file_id'] / self.cases[case]['expression']['file_name']
        
        return {'methylation': self.meth_fpath, 'expression': self.expr_fpath}
    
    
    def get_case_data(self, case: str, genes = None, cpgs = None, get_expr=True, get_meth=True) -> dict:
        """
        Get case methylation and expression dataframe.
        
        :params
            case -- case id;
            genes -- list of genes to get data for. if None: gets data for all genes;
            cpgs -- list of cpgs to get data for. if None: gets data for all cpgs;
        """
        
#         get_expr, get_meth = True, True
        paths = self.get_case_datapaths_(case)
        
        # We ignore rest columns, since terribly redundant. They describe CpG features, which are assumed to be constant for every
        # experiment (to be tested). That's the reason for the inefficient data storage. After testing invariance,
        # methylation data files have to cleared (i.e. remove redundant columns).
        meth_cols = ['Beta_value']  # 'Composite Element REF' becomes index
        
        # load both meth and expr for now
        # if dtype == 'methylation':
        #     get_expr = False
        # elif dtype == 'expression':
        #     get_meth = False
        
        if get_expr:
            if genes is not None:
                expr = pd.read_csv(paths['expression'], sep='\t', index_col=0, header=None)
                expr.index = list(map(lambda x: x.split(".")[0], expr.index))  # cut off version of ids e.g. ENS0000.1 -> ENS0000
                expr = expr.loc[genes]
            else:
                expr = pd.read_csv(paths['expression'], sep='\t', index_col=0, header=None)
        else:
            expr = None
        
        if get_meth:
            if cpgs is not None:
                meth = pd.read_csv(paths['methylation'], sep='\t', index_col=0, header=0)[meth_cols].loc[cpgs]
            else:
                meth = pd.read_csv(paths['methylation'], sep='\t', index_col=0, header=0)[meth_cols]
        else:
            meth = None
        
        return {'methylation': meth, 'expression': expr}