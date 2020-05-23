import sys, os, copy, linecache, time, pickle
import networkx as nx
from extractor.route import Route, RouteManager
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
from extractor.references import amphetamine, _cetirizine, cetirizine, zolpidem
from extractor.run import start_from_pickle, initial_computation

if __name__ == "__main__":
    argv = sys.argv
    reactions = argv[1]
    result_dir = argv[2]
    tree_size = int(argv[3])
    threshold = int(argv[4])
    mode = argv[5]
    print("{}", reactions)
    print("{}", result_dir)
    print("{}", tree_size)
    print("{}", threshold)
    print(mode)
    #method = argv[5]
    true_route = _cetirizine()
    true_route = cetirizine()
    for n in list(true_route.nodes):
        print(true_route.nodes[n])
    #true_route = zolpidem()
    if os.path.exists(os.path.join(result_dir,"route_manager_" + str(tree_size) + ".pickle")):
        route_manager = start_from_pickle(result_dir, tree_size)
    else:
        print("initial computation")
        route_manager = initial_computation(reactions, result_dir, tree_size, threshold, true_route, ["scscore", "step_num"], mode)
    route_manager.rank_and_convert(n=threshold, suffix=str(tree_size))
