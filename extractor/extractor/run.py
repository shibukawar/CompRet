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
from extractor.references import amphetamine, cetirizine

def start_from_pickle(result_dir, tree_size: int):
    pickle_path = os.path.join(result_dir, "route_manager_" + str(tree_size) + ".pickle")
    route_manager = pickle.load(open(pickle_path, 'rb'))
    return route_manager

def route_manager_reader(reaction, result_dir, tree_size, threshold, reference, methods, id_to_string,  mode='enumeration'):
    
    smarts_dict = {}
    G = None
    prev = ""
    f = open(os.path.join(result_dir, "pt.txt"),'r')
    for i in range(tree_size):
        l = f.readline()
        if l == '':
            break
        prev = l
    c = 0
    f.close()
    f = open("./tmp.dot",'w')
    f.write(prev)
    f.close()
    G = read_dot("./tmp.dot")
    linecache.clearcache()
    route_manager = RouteManager(G, id_to_string, result_dir, reference, methods)

    return route_manager

def initial_computation(reaction, result_dir, tree_size, threshold, reference, methods, mode='enumeration'):
    
    smarts_dict = {}
    G = None
    prev = ""
    f = open(os.path.join(result_dir, "pt.txt"),'r')
    for i in range(tree_size):
        l = f.readline()
        if l == '':
            break
        prev = l
    c = 0
    f.close()
    #l = linecache.getline(result_dir+"pt.txt", n)
    f = open("./tmp.dot",'w')
    f.write(prev)
    f.close()
    #os.system("dot -Tpng ./tmp.dot -o nako.png")
    G = read_dot("./tmp.dot")
    linecache.clearcache()
    mol_file = open(os.path.join(result_dir,"finalResult.txt"))
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
    route_manager = RouteManager(G, id_to_string, result_dir, reference, methods)
    print(route_manager.target)
    start = time.time()
    #route_manager.extract_route()
    if mode == 'enumeration':
        route_manager.enumeration()
    elif mode =='sampling':
        route_manager.sampling()
    end = time.time()
    print("Elapsed time for {}:{} sec.".format(mode, end-start))
    #print("the number of route:", len(route_manager.route_list))
    #start = time.time()
    #route_manager.add_name()
    with open (os.path.join(result_dir,"route_manager_" + str(tree_size) + ".pickle"), "wb") as f1:
        pickle.dump(route_manager, f1)
    
    return route_manager