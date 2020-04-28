import pdb
import unittest

import pandas as pd
from sklearn.model_selection import RandomizedSearchCV
from sklearn.datasets import make_regression

from expr_meth_pred import *


X, y = make_regression(n_samples=1000, n_features=100,
                       n_informative=35, n_targets=1,
                       bias=-0.33, effective_rank=23)
X = pd.DataFrame(data=X)
y = pd.Series(data=y)


def test_get_random_genes():
    genes1 = get_random_genes(10)
    genes2 = get_random_genes(100)
    
    assert len(genes1) == 10
    assert type(genes1) == list
    
    assert len(genes2) == 100
    assert type(genes2) == list

def test_get_top_cpgs():
    cpgs1 = get_top_cpgs(10)
    cpgs2 = get_top_cpgs(100)
    
    assert type(cpgs1) == type(cpgs2) == list
    assert len(cpgs1) == 10
    assert len(cpgs2) == 100
    assert cpgs1 == cpgs2[:10]
    assert (cpgs1[:3] == cpgs2[:3] 
            == ["cg14830003", "cg07280731", "cg26475649"])

def test_build_model():
    model = build_model()

def test_hyperparam_opt():
    model = build_model()
    
    param_distributions = {
        "num_leaves": [25,26,30,31, 35],
        "n_estimators": [80, 85, 90, 100, 110, 120],
    }
    
    scoring = "neg_mean_absolute_error" 
    
    search = RandomizedSearchCV(model, 
                                param_distributions, 
                                n_iter=3,
                                scoring=scoring,
                                n_jobs=1)
    
    model = hyperparam_opt(search, X, y)

def test_get_all_projects():
    projs = get_all_projects()
    assert len(projs) == 33
    assert "TCGA-KIRP" in projs
    
def test_case_extraction():
    projs = get_all_projects()
    df_cases = get_projects_cases(projs)
    assert type(df_cases) == pd.DataFrame
    assert df_cases.shape[0] == 9051

def test_split_train_test_cases():
    projs = get_all_projects()
    cases_train, cases_test = split_train_test_cases(projs)
    assert type(cases_train) == type(cases_test) == list
    assert len(cases_test) == 1810

def test_aggregate_case_meth_expr():
    x = pd.DataFrame(data=[[1],[2],[3],[4]], index=['g1','g2','g3','g4'])
    y = pd.DataFrame(data=[[-1],[-2],[-3]], index=['c1','c2','c3'])
    z = aggregate_case_meth_expr(x, y)
    assert type(z) == pd.DataFrame
    assert z.equals(pd.DataFrame(data=[
        [1.0, -1.0, -2.0, -3.0],
        [2.0, -1.0, -2.0, -3.0],
        [3.0, -1.0, -2.0, -3.0],
        [4.0, -1.0, -2.0, -3.0],
    ], columns=[
        0, 'c1', 'c2', 'c3'
    ], index=[
        'g1','g2','g3','g4'
    ]))
    
def test_build_dataset():
#     projs = get_all_projects()[:1]
    projs = ["TCGA-CESC", "TCGA-UCS"]
    genes = get_random_genes(10)
    cpgs = get_top_cpgs(3)
    df = build_dataset(projs, genes, cpgs)
    assert type(df) == pd.DataFrame