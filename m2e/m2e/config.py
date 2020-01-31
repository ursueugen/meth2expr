"""
"""


import os
import json

#fn = '/home/dicer/Documents/Gerontomics/code/meth2expr/m2e/m2e/basicConfig.json'
fn = os.getcwd() + "/../m2e/m2e/basicConfig.json"

with open(fn, "r") as f:
    configs = json.load(f)