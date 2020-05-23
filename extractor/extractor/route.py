import sys, os, copy, linecache, time, heapq
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import rdDepictor, Draw, MolFromSmiles, MolFromSmarts, AllChem
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from rdkit.Chem.rdChemReactions import ReactionFromSmarts, CreateDifferenceFingerprintForReaction, CreateStructuralFingerprintForReaction
from rdkit.Chem.Fingerprints import FingerprintMols
from networkx.drawing.nx_pydot import read_dot
from networkx.algorithms import dag_longest_path_length
from networkx.algorithms.matching import is_perfect_matching
from typing import List, Dict
from itertools import permutations
from tqdm import tqdm
from standalone_model_numpy import SCScorer


class Route(nx.DiGraph):
    
    def __init__(self, incoming_graph_data=None, **attr):
        super(Route, self).__init__(incoming_graph_data=None, **attr)
        self.scores = {}
        self.intermediate_products = []
        self.starting_materials = []
        self.name = "route"
        self.ignore = False
        self.sampled_count = 0

    def __lt__(self, other):
        return len(self.nodes) < len(other.nodes)
    
    def __le__(self, other):
        return len(self.nodes) <= len(other.nodes)

    def __eq__(self, other):
        return len(self.nodes) == len(other.nodes)

    def __ne__(self, other):
        return len(self.nodes) != len(other.nodes)
    
    def __gt__(self, other):
        return len(self.nodes) > len(other.nodes)
    
    def __ge__(self, other):
        return len(self.nodes) >= len(other.nodes)

    def is_included(self, query: str, option="id"):
        #mol = Chem.MolFromSmiles(query)
        # TODO:変な入力を弾く
        if option == "id":
            return query in self.nodes()
        elif option == "smiles":
            return query in self.intermediate_products
    
    def get_leaf_nodes(self):
        #tmp = [n for n in self.nodes() if self.out_degree(n) == 0]
        return [n for n in self.nodes() if self.out_degree(n) == 0]
    
    def store_intermediate_products(self, G):
        # After completing route extraction, store all intermediate_products (root is excluded)
        for n in self.nodes:
            if n == '0':
                continue
            else:
                try:
                    sml = G.nodes[n]['smiles']
                    mol = Chem.MolFromSmiles(sml)
                    self.intermediate_products.append((sml, int(mol.GetNumAtoms())))
                except:
                    continue
        self.intermediate_products.sort(key=lambda x: x[1])

    def _store_intermediate_products(self):
        # After completing route extraction, store all intermediate_products (root is excluded)
        for n in self.nodes:
            if self.nodes[n]['node_type'] == "OR":
                if n == '0':
                    continue
                sml = self.nodes[n]['smiles']
                #print(sml)
                #print("hoge")
                mol = Chem.MolFromSmiles(sml)
                self.intermediate_products.append((sml, int(mol.GetNumAtoms())))
        self.intermediate_products.sort(key=lambda x: x[1])
        
class UniqueHeap(object):

    def __init__(self, _max_size=1000, _method="similarity"):
        self.h = []
        self.step_dict = {}
        self.max_size = _max_size
        self.method = _method
    
    def push(self, r: Route):
        if r.ignore:
            print('ignore')
            return
        tpl = (r.scores[self.method], r)
        if self.method == 'step_num':
            if not tpl[0] in self.step_dict.keys():
                self.step_dict[tpl[0]] = 1
            else:
                if self.step_dict[tpl[0]] < 10:
                    self.step_dict[tpl[0]] += 1
                else:
                    return
        else:
            for i in range(len(self.h)):
                if (self.h[i][0]-tpl[0])**2 < 0.00001:
                    return
        heapq.heappush(self.h, tpl)
        if len(self.h) > self.max_size:
            heapq.heappop(self.h)
            assert len(self.h) == self.max_size

    def sort(self):
        if self.method == "step_num":
            self.h = sorted(self.h, key=lambda x: x[0], reverse=False)
        else:
            self.h = sorted(self.h, key=lambda x: x[0], reverse=True)

    def get(self, index):
        return self.h[index]

class RouteManager(object):
    
    def __init__(self, _G, 
                       _id_dict,
                       _result_dir: str,
                       _ref=None,
                       _methods=["scscore", "similarity", "step_num"],
                       _project_root="/Users/shibukawar/Documents/lab/project/synthesis/CompRet/scscore/"):
        self.G = _G
        self.project_root = _project_root
        min_separation = 0.25
        model = SCScorer()
        model.restore(os.path.join(self.project_root, 'models', 'full_reaxys_model_1024uint8', 'model.ckpt-10654.as_numpy.json.gz'))
        self.id_dict = _id_dict
        for k in list(self.G.nodes):
            if self.G.nodes[k]['shape'] == 'circle':
                self.G.nodes[k]['smiles'] = self.id_dict[k]
                try:
                    self.G.nodes[k]['scscore'] = model.get_score_from_smi(self.id_dict[k])[1]
                except:
                    self.G.nodes[k]['scscore'] = 0
            else:
                self.G.nodes[k]['smarts'] = self.id_dict[k]
        self.fp_dict = {}
        for n in list(self.G.nodes):
            if self.G.nodes[n]['shape'] == 'box':
                sml = self.G.nodes[n]['smarts']
                if not sml in self.fp_dict.keys():
                    rxn = ReactionFromSmarts(sml)
                    fp_vec = CreateStructuralFingerprintForReaction(rxn)
                    self.fp_dict[sml] = np.array(list(fp_vec))
        self.fp_len = 4096
        self.route_list = []
        self.route_lists = {}
        for m in _methods:
            self.route_lists[m] = UniqueHeap(_method=m)
        self.methods = _methods
        self.target = self.id_dict['0']
        self.route_name_to_index = {}
        self.result_dir = _result_dir
        self.route_suffix = ''
        self.route_num_threshold = 0
        self.route_variation_dict = {}
        self.route_vectors = []
        self.reference = _ref
        self.count = 0
        self.s_time = 0
        self.calc_scores_time = 0
        self.ratio = 0.5
        self.prohibited = list()
        self.start_time = 0
    
    def route_to_vec(self):
        for r in self.route_list:
            fp_route = np.zeros(self.fp_len)
            for n in list(r.nodes):
                if self.G.nodes[n]['shape'] == 'box':
                    sml = self.G.nodes[n]['smarts']
                    fp_route += self.fp_dict[sml]
            self.route_vectors.append(fp_route)
            
    def route_to_vec(self, idx):
        for i in idx:
            r = self.route_list[i]
            fp_route = np.zeros(self.fp_len)
            for n in list(r.nodes):
                if self.G.nodes[n]['shape'] == 'box':
                    sml = self.G.nodes[n]['smarts']
                    fp_route += self.fp_dict[sml]
            self.route_vectors.append(fp_route)



    def get_statistics(self, problem_name):
        or_node_num = 0
        and_node_num = 0
        for n in self.G.nodes:
            if self.G.nodes[n]['shape'] == 'circle':
                or_node_num += 1
            elif self.G.nodes[n]['shape'] == 'box':
                and_node_num += 1
        edge_num = len(list(self.G.edges))
        route_num = len(self.route_list)
        print("{} & {} & {} & {} & {} & {}".format(problem_name, self.target, or_node_num, and_node_num, edge_num, route_num))
    
    def sampling(self, sampling_size=500000):
        g = Route()
        g.add_node('0')
        self.route_list = []
        _neighs_num = len(list(self.G.neighbors('0')))
        self.route_num_threshold = max(10000, int(sampling_size / _neighs_num))
        print("route_num_threshold:", self.route_num_threshold)
        for _n in self.G.neighbors('0'):
            self.route_variation_dict[_n] = 0
        self.sampler(g, '0')

    def sampler(self, g, n):
        if len(list(g)) > 1 and self.route_variation_dict[list(g.nodes)[1]] > self.route_num_threshold:
            return
        neighbors = list(self.G.neighbors(n))
        neighbors_to_search = [x for x in neighbors if not x in self.prohibited]
        if len(neighbors_to_search) == 0:
            for node in neighbors:
                self.prohibited.remove(node)
            return
        g_ = copy.deepcopy(g)
        p = neighbors_to_search[0]
        self.prohibited.append(p)
        sampling_threshold = float(5 / len(neighbors))
        v = np.random.uniform(0,1)
        continue_searching = False
        if n != '0' and g.sampled_count == 0 and v < sampling_threshold:
            g.sampled_count += 1
            continue_searching = True
        if n == '0' or (n != '0' and (g.sampled_count == 1 or continue_searching)):
            g.add_node(p)
            g.add_edge(n, p)
            for c in list(self.G.neighbors(p)):
                g.add_node(c)
                g.add_edge(p, c)
            leaves = g.get_leaf_nodes()
            leaves = [x for x in leaves if len(list(self.G.neighbors(x))) != 0]
            if len(leaves) > 0:
                t = leaves[0]
                self.sampler(g, t)
            else:
                self.route_list.append(g)
                self.route_variation_dict[list(g.nodes)[1]] += 1
        self.sampler(g_, n)

    def enumeration(self):
        g = Route()
        g.add_node('0')
        self.start_time = time.time()
        edge_list = []
        self.partition('0', edge_list)
        print("num:", len(self.route_list))

    def _partition(self, g, n):
        if time.time() - self.start_time > 1800:
            return
        if len(self.route_list) % 10000 == 0:
            print(len(self.route_list))
        neighbors = list(self.G.neighbors(n))
        neighbors_to_search = [x for x in neighbors if not x in self.prohibited]
        if len(neighbors_to_search) == 0:
            for node in neighbors:
                self.prohibited.remove(node)
            return
        g_ = copy.deepcopy(g)
        p = neighbors_to_search[0]
        g.add_node(p)
        g.add_edge(n, p)
        for c in list(self.G.neighbors(p)):
            g.add_node(c)
            g.add_edge(p, c)
        self.prohibited.append(p)
        leaves = g.get_leaf_nodes()
        leaves = [x for x in leaves if len(list(self.G.neighbors(x))) != 0]
        if len(leaves) > 0:
            t = leaves[0]
            self._partition(g, t)
        else:
            self.route_list.append(g)
        self._partition(g_, n)

    def partition(self, n, edge_list):
        if time.time() - self.start_time > 1800:
            return
        if len(self.route_list) > 0 and len(self.route_list) % 10000 == 0:
            print(len(self.route_list))
        neighbors = list(self.G.neighbors(n))
        neighbors_to_search = [x for x in neighbors if not x in self.prohibited]
        if len(neighbors_to_search) == 0:
            for node in neighbors:
                self.prohibited.remove(node)
            return
        p = neighbors_to_search[0]
        current_edge_list = []
        for e in edge_list:
            current_edge_list.append(e)
        edge_list.append((n, p))
        for c in list(self.G.neighbors(p)):
            edge_list.append((p, c))
        g = Route()
        for e in edge_list:
            g.add_edge(e[0], e[1])
        self.prohibited.append(p)
        try:
            _ = dag_longest_path_length(g)
            all_nodes = set()
            sources = set()
            for e1 in edge_list:
                all_nodes.add(e1[0])
                sources.add(e1[0])
                all_nodes.add(e1[1])
            leaves = list(all_nodes.difference(sources))
            leaves = [x for x in leaves if len(list(self.G.neighbors(x))) != 0]
            if len(leaves) > 0:
                t = leaves[0]
                self.partition(t, edge_list)
            else:
                self.route_list.append(g)
        except:
            pass
        self.partition(n, current_edge_list)


    def is_all_leaf_nodes_terminal(self, leaf_nodes):
        for n in leaf_nodes:
            if len(list(self.G.neighbors(n))) != 0:
                return False
        return True
    
            
    def calc_score_(self, r, method):
        if method == "similarity":
            r.store_intermediate_products(self.G)
            trg = copy.deepcopy(r.intermediate_products)
            ref = copy.deepcopy(self.reference.intermediate_products)
            n_ref = len(ref)
            n_trg = len(trg)
            if n_ref > n_trg:
                tmp = n_trg
                n_trg = n_ref
                n_ref = tmp
                tmp = copy.deepcopy(trg)
                trg = copy.deepcopy(ref)
                ref = tmp
            score_tmp = 0
            for j in range(n_ref):
                m = Chem.MolFromSmiles(ref[j][0])
                m_trg = Chem.MolFromSmiles(trg[j][0])
                ms = [m, m_trg]
                fps = [FingerprintMols.FingerprintMol(x) for x in ms]
                s = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                score_tmp += s
            r.scores["similarity"] = score_tmp
        elif method == "scscore":
            r.scores["scscore"] = np.mean([5-self.G.nodes[x]['scscore'] for x in list(r.nodes) if 'scscore' in self.G.nodes[x].keys()])
        elif method == "step_num":
            try:
                r.scores["step_num"] = dag_longest_path_length(r) / 2
            except:
                print("{} has been flagged as ignore because dag_longest_path_length(r) could not be computed.".format(r.name))
                r.scores["step_num"] = 1000
                r.ignore = True
        self.route_lists[method].push(r)

    # calculation scores for step_num, scscore, similarity
    def calc_scores(self, r):
        start = time.time()
        for m in self.methods:
            self.calc_score_(r, m)
        t = time.time() - start
        self.calc_scores_time += t

    def add_name(self):
        if not os.path.exists(self.result_dir):
            os.mkdir(self.result_dir)
        count: int = 0
        for r in self.route_list:
            dir_name = str(count).zfill(9) + "/"
            r.name += str(count).zfill(9)
            count += 1
    
    def calc_all_scores(self):
        for r in tqdm(self.route_list):
            self.calc_scores(r)
    
    def create_image(self, base_dir, label):
        if 'smiles' in self.G.nodes[label].keys():
            m = MolFromSmiles(self.G.nodes[label]['smiles'])
            img = Draw.MolsToGridImage([m], subImgSize=(200, 200), molsPerRow=1)
            if len(self.G.nodes[label]['smiles']) < 10:
                img = Draw.MolToImage(m, size=(200, 200))
            img.save(os.path.join(base_dir, os.path.join('mols', label+'.png')))
        elif 'smarts' in self.G.nodes[label].keys():
            rxn = AllChem.ReactionFromSmarts(self.G.nodes[label]['smarts'])
            img = Draw.ReactionToImage(rxn)
            img.save(os.path.join(base_dir, os.path.join('rxns', label+'.png')))

    def route_to_dot(self, route_dir: str, threshold: int):
        dots_dir = os.path.join(route_dir, "dots")
        mol_img_dir = os.path.join(route_dir, "mols")
        rxn_img_dir = os.path.join(route_dir, "rxns")
        if not os.path.exists(dots_dir):
            os.mkdir(dots_dir)
        if not os.path.exists(mol_img_dir):
            os.mkdir(mol_img_dir)
        if not os.path.exists(rxn_img_dir):
            os.mkdir(rxn_img_dir)
        for method in self.methods:
            self.route_lists[method].sort()
            f = open(os.path.join(route_dir, method + ".txt"), 'w')
            loop_num = min(threshold, len(self.route_lists[method].h))
            for i in tqdm(range(loop_num)):
                (val, r) = self.route_lists[method].get(i)
                if r.ignore:
                    continue
                route_name = r.name
                f.write(str(r.scores[method]) + ' ' + str(route_name) + '\n')
                tmp = os.path.join(route_dir ,route_name)
                for n in list(r.nodes):
                    self.create_image(route_dir, n)
                    img_path = os.path.join(mol_img_dir, str(n) + ".png")
                    rxn_path = os.path.join(rxn_img_dir, str(n) + ".png")
                    if os.path.exists(img_path):
                        r.nodes[n]['image'] = img_path
                        r.nodes[n]['label'] = " "
                        r.nodes[n]['shape'] = 'circle'
                    elif os.path.exists(rxn_path):
                        r.nodes[n]['image'] = rxn_path
                        r.nodes[n]['label'] = " "
                        r.nodes[n]['shape'] = 'box'
                A = nx.nx_agraph.to_agraph(r)
                A.write(os.path.join(dots_dir, r.name + '.dot')) 

    def rank_and_convert(self, n=100, suffix=''):
        self.route_suffix = suffix
        route_dir = os.path.join(self.result_dir, 'route'+suffix)
        print(route_dir)
        if not os.path.exists(route_dir):
            #print(route_dir)
            os.mkdir(route_dir)
        # rank top n routes in terms of score
        start = time.time()
        self.add_name()
        print("Done: add_name()")
        self.calc_all_scores()
        print("Done: calc_all_scores()")
        self.route_to_dot(route_dir, n)
        print("Done: route_to_dot()")
        end = time.time()
        print("Elapsed time for convert to dot:{} sec.".format(end-start))

