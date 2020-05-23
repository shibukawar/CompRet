import sys, os, copy, linecache, time, pickle
import networkx as nx
from route import Route, RouteManager
#import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit import rdBase
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from networkx.drawing.nx_pydot import read_dot
from networkx.algorithms import dag_longest_path_length
from networkx.algorithms.matching import is_perfect_matching
from typing import List, Dict
from itertools import permutations

if __name__ == "__main__":
    argv = sys.argv
    reactions = argv[1]
    result_dir = argv[2]
    n = int(argv[3])
    with open(os.path.join(result_dir, "route_manager.pickle"), 'rb') as f:
        route_manager = pickle.load(f)
    print(len(route_manager.route_list))
    route_manager.convert_route_to_smarts(result_dir, os.path.join(result_dir, "route/ranking.txt"), n)
