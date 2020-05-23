import sys, os, copy, linecache, time, pickle
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from route import Route, RouteManager
from references import cetirizine
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw, MolFromSmiles, MolFromSmarts, AllChem
from rdkit import rdBase
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from networkx.drawing.nx_pydot import read_dot
from networkx.algorithms import dag_longest_path_length
from networkx.algorithms.matching import is_perfect_matching
from networkx.algorithms.cycles import simple_cycles, find_cycle
from networkx.drawing.nx_pylab import *
from networkx.drawing.nx_agraph import *
from typing import List, Dict
from itertools import permutations
from references import amphetamine, cetirizine
from run import start_from_pickle, initial_computation, route_manager_reader

result_dir = "/home/shibukawa/exec_log/cetirizine_reaxys_6_10800000_-1/"
reactions = "/home/shibukawa/ShibuChem/retrosynthesis/reaxys_reactions"
threshold = 100
mode = "enumeration"
result_dir = "/home/shibukawa/exec_log/cetirizine_reaxys_6_10800000_-1/"
pt = os.path.join(result_dir, "pt.txt")
timePoint = os.path.join(result_dir, "timePoint.txt")
mol_file = open(os.path.join(result_dir,"finalResult.txt"))
pt_file = open(pt, 'r')
tp_file = open(timePoint, 'r')
time_points = tp_file.readlines()
time_points = [t.split(':')[1].split(',')[0] for t in time_points]
tp_file.close()
print(len(time_points))
id_to_string = {}
for l in mol_file:
    if l.startswith("digraph") or l.startswith("l") or l.startswith("in") or l.startswith("end") or l.startswith("$.digraph"):
        continue
    else:
        try:
            l = l.rstrip("\n").split(",")
            node_id = l[0]
            string = ''.join(l[1:])
            if len(string) == 0:
                continue
            id_to_string[node_id] = string
        except:
            pass
route_nums = []
total_times = []
for i in range(len(time_points)):
    G = None
    tree_size = i + 1
    true_route = cetirizine()
    route_manager = route_manager_reader(reactions, result_dir, tree_size, threshold, true_route, ["similarity", "scscore", "step_num"], id_to_string, mode)
    start = time.time()
    route_manager.enumeration()
    end = time.time()
    num = len(route_manager.route_list)
    route_nums.append(num)
    total_times.append(int(time_points[i])/1000+float(end-start))
    print(total_times[i], num)
    if num > 500000:
        break