"""
Script for running the experiment on predicting gene expression
 from methylation of CpG sites based on sites found to have the
 highest correlation with individual genes' expressions.
"""

from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.model_selection import RandomizedSearchCV

from m2e.project import Project


# input
NUM_METH = 50  # top 50 with predictive power
NUM_GENES = 100  # num samples = (this) * nr_cases
SEARCH_SPACE = None

# const
DATA_PATH = Path("data/")
LOOKUP_PATH = DATA_PATH / "genomics/gene_id_lookup.csv"
CPG_CORR_PATH = DATA_PATH / "broad_tcga/analysis/gdac.broadinstitute.org_STAD-TP.Correlate_Methylation_vs_mRNA.Level_4.2016012800.0.0/Correlate_Methylation_vs_mRNA_STAD-TP_matrix.txt"



def get_random_genes(num: int) -> List[str]:
    """
    Return: list of random ncbi gene ids of length num.
    """
    df = pd.read_csv(LOOKUP_PATH)
#     col = "NCBI gene ID"
    col = "Gene stable ID"
    ids = df[col].to_list()
    choice = np.random.choice(ids, size=num, replace=False)
    return choice.tolist()
    
def get_top_cpgs(num: int) -> List[str]:
    """
    Return: list of cpgs with num highest absolute coefficients
            for gene expression correlation from the firehose analysis.
    """
    df = pd.read_csv(CPG_CORR_PATH, sep='\t')
    cpg_col = 'Meth_Probe'
    corr_col = 'Corr_Coeff'
    df['abs_corr_coef'] = abs(df[corr_col])
    df.sort_values('abs_corr_coef', ascending=False, inplace=True)
    return df[cpg_col][:num].to_list()


def build_model() -> lgb.LGBMRegressor:
    """
    Returns: model with defaults, to be optimized.
    """
    model = lgb.LGBMRegressor()
    return model

def hyperparam_opt(search: RandomizedSearchCV,
                   X_train: pd.DataFrame, 
                   y_train: pd.Series) -> lgb.LGBMRegressor:
    """
    Args:
        search: search instance with scikit-learn Search API.
        {X,y}_train
    Returns:
        instance of model configured with the best found parameters.
    """
    search.fit(X_train, y_train)
    return search.best_estimator_


def get_all_projects() -> List[str]:
    """
    Returns: list of projects names.
    """
    projects_dir = Path("..") / "tcga/projects/"
    assert projects_dir.is_dir()
    projects = [p.name for p in projects_dir.glob("*/")]
    return projects
    
def get_projects_cases(projects: List[str]) -> pd.DataFrame:
    """
    Returns: list of cases from the selected projects.
    """
    cases = {}
    for p in [Project(i) for i in projects]:
        for c in p.cases.keys():
            cases[c] = p.name
    return pd.DataFrame.from_dict(cases, orient="index")

def split_train_test_cases(projects: List[str]) -> (
        Tuple[List[str], List[str]]):
    """
    Splits cases in train and test sets.
    
    Returns: 2-tuple of train and test list.
    """
    cases = get_projects_cases(projects)
    case_counts = cases[0].value_counts()
    projs_test = case_counts.index[-17:]
    
    cases_test = cases.index[cases[0].isin(projs_test)].to_list()
    cases_train = cases.index[~cases[0].isin(projs_test)].to_list()
    return cases_train, cases_test

def aggregate_case_meth_expr(df_expr: pd.DataFrame,
                            df_meth: pd.DataFrame) -> (
                            pd.DataFrame):
    """
    Aggregates expr and meth dataframes per case.
    
    Args:
        df_expr: expression dataframe of format (genes x expression)
        df_meth: methylation dataframe of format (cpgSite x methylation)
    
    Returns:
        dataframe of format (gene x [expression, cpg1, cpg2, ...])
    """
    assert df_expr.shape[1] == 1
    assert df_meth.shape[1] == 1
    
    array = np.zeros(shape=(df_expr.shape[0], 1 + df_meth.shape[0]))
    
    array[:, 0] = df_expr.iloc[:, 0]
    array[:, 1:] = df_meth.iloc[:, 0]
    
    df = pd.DataFrame(data=array, 
                      index=df_expr.index, 
                      columns=([df_expr.columns[0]]
                               + df_meth.index.to_list()))
    return df

def build_dataset(projects: List[str], genes: List[str],
                   cpgs: List[str]) -> pd.DataFrame:
    """
    Builds the dataset used for the experiment.
    
    Args:
        projects: list of projects to include in dataset.
        genes: list of genes for which to extract expression.
        cpgs: list of CpG sites for which to extract methylation.
        
    Returns:
        DataFrame with (case_id, gene_id) x (gene_expression, cpg_meth1, ...)
    """
#     raise NotImplementedError("Multithreading for IO-bounded runtime")
    
    dfs = []
    for p_name in projects:
        proj = Project(p_name)
        
        for case in proj.cases.keys():    
            case_data = proj.get_case_data(case, genes=genes, cpgs=cpgs)
            df_expr = case_data['expression']
            df_meth = case_data['methylation']
            
            assert df_expr.shape[0] == len(genes)
            assert df_meth.shape[0] == len(cpgs)
            
            df_case = aggregate_case_meth_expr(df_expr, df_meth)
            
            # add multilevel index
            df_case.index = pd.MultiIndex.from_tuples(
                [(case, gene) for gene in df_case.index],
                names=['case', 'gene'])
            
            # add colname for expression
            assert df_case.columns.to_list()[1:] == cpgs
            df_case.columns = ['expression'] + cpgs
            
            dfs.append(df_case)
    
    df_final = pd.concat(dfs, axis=0, join='inner')
    assert df_final.columns.to_list() == ['expression'] + cpgs
    
    return df_final
    
    
if __name__ == '__main__':
    
    genes = get_random_genes(NUM_GENES)
    cpgs = get_top_cpgs(NUM_METH)
    
    tcga_projs = get_all_projects()
    cases_train, cases_test = split_train_test_cases(tcga_projs)
    
    df = build_dataset(tcga_projs, genes, cpgs)
    df_train = df.iloc[df.loc['case'].isin(cases_train), :].copy()
    df_test = df.iloc[df.loc['case'].isin(cases_test), :].copy()
    
    X_train = df_train.iloc[:, 1:]
    y_train = df_train.iloc[:, 0]
    X_test = df_test.iloc[:, 1:]
    y_test = df_test.iloc[:, 0]
    
    model = build_model()
    model = hyperparam_opt(SEARCH_SPACE, X_train, y_train)
    model.fit(X_train, y_train)
    
    
#     results = {"R2", "MSE", "MAE", "?"...}
#     # can iterate and have distr