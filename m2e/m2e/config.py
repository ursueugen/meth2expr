"""
"""


import json


CHROMOSOMES = list(map(lambda x: str(x), range(1, 23))) + ["X", "Y", "MT"]
with open("basicConfig.json", "r") as f:
    configs = json.load(f)