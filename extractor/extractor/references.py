from extractor.route import Route, RouteManager

def amphetamine():
    true_route_for_amphetamine = Route()
    sml1 = "CC(N(=O)=O)=CC1=CC=CC=C1"
    sml2 = "[H]C(C1=CC=CC=C1)=O"
    sml3 = "CCN(=O)=O"
    true_route_for_amphetamine.add_node('0', name='0', smiles="CC(N)CC1=CC=CC=C1", node_type="OR")
    true_route_for_amphetamine.add_node('1', name="hoge", smarts="fuga", node_type="AND")
    true_route_for_amphetamine.add_node('2', name='0', smiles=sml1, node_type="OR")
    true_route_for_amphetamine.add_node('3', name="piyo", smarts="hera", node_type="AND")
    true_route_for_amphetamine.add_node('4', name='0', smiles=sml2, node_type="OR")
    true_route_for_amphetamine.add_node('5', name='0', smiles=sml3, node_type="OR")
    true_route_for_amphetamine.add_edge('0','1')
    true_route_for_amphetamine.add_edge('1','2')
    true_route_for_amphetamine.add_edge('2','3')
    true_route_for_amphetamine.add_edge('3','4')
    true_route_for_amphetamine.add_edge('3','5')
    true_route_for_amphetamine._store_intermediate_products()
    return true_route_for_amphetamine

def _cetirizine():
    true_route_for_cetirizine = Route()
    trg = "ClC(C=C1)=CC=C1C(N2CCN(CCOCC(O)=O)CC2)C3=CC=CC=C3"
    smiles1 = "COC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1"
    smiles2 = "COC(=O)COCCCl"
    smiles3 = "OC(=O)COCCCl"
    smiles4 = "Clc1ccc(cc1)C(N1CCNCC1)c1ccccc1"
    smiles5 = "C1CNCCN1"
    smiles6 = "ClC(c1ccccc1)c1ccc(Cl)cc1"
    true_route_for_cetirizine.add_node('0', name='0', smiles=trg, node_type="OR")
    #true_route_for_cetirizine.add_edge('0','1')
    true_route_for_cetirizine.add_node('1', name='1', smiles=smiles1, node_type="OR")
    true_route_for_cetirizine.add_node('2', name='2', smiles=smiles2, node_type="OR")
    true_route_for_cetirizine.add_node('3', name='3', smiles=smiles3, node_type="OR")
    true_route_for_cetirizine.add_node('4', name='4', smiles=smiles4, node_type="OR")
    true_route_for_cetirizine.add_edge('0','1')
    true_route_for_cetirizine.add_edge('1','2')
    true_route_for_cetirizine.add_edge('1','4')
    #true_route_for_cetirizine.add_edge('2','3')
    #true_route_for_cetirizine.add_edge('4','5')
    #true_route_for_cetirizine.add_edge('4','6')
    #true_route_for_cetirizine.add_edge('4','6')
    true_route_for_cetirizine._store_intermediate_products()
    return true_route_for_cetirizine

def cetirizine():
    true_route_for_cetirizine = Route()
    trg = "ClC(C=C1)=CC=C1C(N2CCN(CCOCC(O)=O)CC2)C3=CC=CC=C3"
    smiles1 = "COC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1"
    smiles2 = "COC(=O)COCCCl"
    smiles3 = "OC(=O)COCCCl"
    smiles4 = "Clc1ccc(cc1)C(N1CCNCC1)c1ccccc1"
    smiles5 = "C1CNCCN1"
    smiles6 = "ClC(c1ccccc1)c1ccc(Cl)cc1"
    smiles7 = "CO"
    true_route_for_cetirizine.add_node('0', name='0', smiles=trg, node_type="OR")
    #true_route_for_cetirizine.add_edge('0','1')
    true_route_for_cetirizine.add_node('1', name='1', smiles=smiles1, node_type="OR")
    #true_route_for_cetirizine.add_edge('1','2')
    true_route_for_cetirizine.add_node('2', name='2', smiles=smiles2, node_type="OR")
    true_route_for_cetirizine.add_node('3', name='3', smiles=smiles3, node_type="OR")
    #true_route_for_cetirizine.add_edge('1','3')
    true_route_for_cetirizine.add_node('4', name='4', smiles=smiles4, node_type="OR")
    #true_route_for_cetirizine.add_edge('2','4')
    true_route_for_cetirizine.add_node('5', name='5', smiles=smiles5, node_type="OR")
    true_route_for_cetirizine.add_node('6', name='6', smiles=smiles6, node_type="OR")
    true_route_for_cetirizine.add_node('7', name='7', smiles=smiles7, node_type="OR")
    true_route_for_cetirizine.add_edge('0','1')
    true_route_for_cetirizine.add_edge('1','2')
    true_route_for_cetirizine.add_edge('1','4')
    #true_route_for_cetirizine.add_edge('2','3')
    #true_route_for_cetirizine.add_edge('4','5')
    #true_route_for_cetirizine.add_edge('4','6')
    #true_route_for_cetirizine.add_edge('4','6')
    true_route_for_cetirizine._store_intermediate_products()
    return true_route_for_cetirizine

def zolpidem():
    true_route_for_zolpidem = Route()
    trg = "CC(C=C1)=CN2C1=NC(C3=CC=C(C)C=C3)=C2CC(N(C)C)=O"
    smiles1 = "CC(C=C1)=CN2C1=NC(C3=CC=C(C)C=C3)=C2C[N+]#[C-]"
    smiles2 = "CC(C=C1)=CN2C1=NC(C3=CC=C(C)C=C3)=C2CN(C)C"
    smiles3 = "CC(C=C1)=CN2C1=NC(C3=CC=C(C)C=C3)=C2"
    smiles4 = "CC1=CC=C(C(CBr)=O)C=C1"
    smiles5 = "CC1=CN=C(N)C=C1"
    true_route_for_zolpidem.add_node('0', name='0', smiles=trg, node_type="OR")
    true_route_for_zolpidem.add_node('1', name='1', smiles=smiles1, node_type="OR")
    true_route_for_zolpidem.add_node('2', name='2', smiles=smiles2, node_type="OR")
    true_route_for_zolpidem.add_node('3', name='3', smiles=smiles3, node_type="OR")
    #true_route_for_cetirizine.add_edge('1','3')
    true_route_for_zolpidem.add_node('4', name='4', smiles=smiles4, node_type="OR")
    #true_route_for_cetirizine.add_edge('2','4')
    true_route_for_zolpidem.add_node('5', name='5', smiles=smiles5, node_type="OR")
    #true_route_for_cetirizine.add_node('6', name='6', smiles=smiles1, node_type="OR")
    true_route_for_zolpidem.add_edge('0','1')
    true_route_for_zolpidem.add_edge('1','2')
    true_route_for_zolpidem.add_edge('2','3')
    true_route_for_zolpidem.add_edge('3','4')
    true_route_for_zolpidem.add_edge('3','5')
    #true_route_for_cetirizine.add_edge('2','3')
    #true_route_for_cetirizine.add_edge('4','5')
    #true_route_for_cetirizine.add_edge('4','6')
    #true_route_for_cetirizine.add_edge('4','6')
    true_route_for_zolpidem._store_intermediate_products()
    return true_route_for_zolpidem
