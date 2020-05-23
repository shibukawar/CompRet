package retrosynthesis;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.TimeUnit;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

import chemaxon.calculations.clean.fixers.length.BondLengthChecker;
import chemaxon.checkers.StructureChecker;
import chemaxon.checkers.ValenceErrorChecker;
import chemaxon.checkers.result.StructureCheckerResult;
import chemaxon.formats.*;
import chemaxon.marvin.io.MolExportException;
import chemaxon.reaction.*;
import com.chemaxon.mapper.AutoMapper;
import chemaxon.struc.*;
import chemaxon.struc.RxnMolecule;
import chemaxon.struc.interfaces.chemicalgraph.BasicStructure;

public class Tree implements Serializable{
  public Node root;
  public List<Molecule> reactionData = new ArrayList<Molecule>();
  public List<String> reactionDataByName = new ArrayList<String>();
  public Set<String> moleculeFilter = new HashSet<String>();
  public List<String> dictionary = new ArrayList<String>();
  public List<Route> synthesisRoutes = new ArrayList<Route>();
  public List<Node> pathFromRoot = new ArrayList<Node>();
  public Map<Integer, String> availableCompoundList = new HashMap<Integer, String>();
  public static Map<String,Node> nodeHashMap = new HashMap<String,Node>();
  public boolean availabilityCheck = false;
  public boolean foundLastNode = false;
  public int maxDepth = 100;
  public int id = 0;
  public int fileNameCounter = 0;
  public int smilesLengthThreshold = 0;
  public int phiChild = 0;
  public int deltaSecond = 0;
  public String reactionLibrary;
  public String buildingBlocks;
  public int inf = Integer.MAX_VALUE / 2;
  public long startTime = 0;
  public long endTime = 0;
  public boolean stopByTime = true;
  private long timeThreshold = 600000;
  private long routeNumThreshold = -1;
    
  public Tree(Molecule m,boolean check,int maxDepth, long timeThreshold, long routeNumThreshold, String reactionLibraryName, String buildingBlocksName) throws IOException, ReactionException {
    String sml = MolExporter.exportToFormat(m,"smiles:u");
    this.root = new Node(m,"molecule",sml,0,id);
    this.root.isRoot = true;
    this.nodeHashMap.put(sml, this.root);
    this.smilesLengthThreshold = sml.length();
    this.maxDepth = maxDepth;
    this.reactionLibrary = reactionLibraryName;
    this.buildingBlocks = buildingBlocksName;
    double falsePositiveProbability = 0.01;
    int expectedNumberOfElements = 1000000000;
    if (timeThreshold == -1) {
      this.stopByTime = false;
    }
    this.timeThreshold = timeThreshold;
    this.routeNumThreshold = routeNumThreshold;
    availabilityCheck = check;
    System.out.println("Initializing reaction database...");
    reactionDataInitializer();
    System.out.println("Done.");
    if (availabilityCheck) {
      System.out.println("Initializing compounds database...");
      compoundsDataInitializer();
      System.out.println("Done.");
    }
  }
    
  public void reactionDataInitializer() throws IOException, ReactionException {
    String prefix = this.reactionLibrary;
    System.out.printf("As reaction data, using %s\n", prefix);
    File dir = new File(prefix);
    String[] files = dir.list();
    for (int i = 0; i < files.length; ++i) {
      if (files[i].equals(".DS_Store")) {
        continue;
      }
      if (i % 1000 == 0) {
        System.out.println(i);
      }
      reactionDataByName.add(files[i]);
      byte[] sBytes = Files.readAllBytes(Paths.get(prefix+files[i]));
      String tmp = new String(sBytes,StandardCharsets.UTF_8);
      try {
        Molecule rxn = MolImporter.importMol(tmp, "smarts");
        AutoMapper mapper = new AutoMapper();
        mapper.mapReaction(rxn);
        reactionData.add(rxn);
      } finally {
        continue;
      }
    }
  }

  public void compoundsDataInitializer() {
    String filePath = this.buildingBlocks;
    FileReader fr = null;
    BufferedReader br = null;
    try {
      fr = new FileReader(filePath);
      br = new BufferedReader(fr);
      String line;
      System.out.println("StartBlocks");
      while ((line = br.readLine()) != null) {
        Molecule tmp = MolImporter.importMol(line,"smiles");
        String key = MolExporter.exportToFormat(tmp, "smiles:u");
        System.out.println(line + "," + key);
        this.moleculeFilter.add(key);
      }
      System.out.println("EndBlocks");
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    } finally {
      try {
        br.close();
        fr.close();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
    
    
  }
    
  public boolean commerciallyAvailableCompound(Molecule m) throws IOException {
    String sml = MolExporter.exportToFormat(m,"smiles:u");
    return this.moleculeFilter.contains(sml);
  }
  
  static Boolean isValidValence(Molecule mol) {
    ValenceErrorChecker checker = new ValenceErrorChecker();
    StructureCheckerResult result = checker.check(mol);
    return result == null;
  }

  static Boolean isValidReactionResultByAtom(Molecule product, Molecule[] reactants, int threshold) {
    int productAtomCount = product.getAtomCount();
    int reactantTotalAtomCount = 0;
    for (Molecule m : reactants) {
      reactantTotalAtomCount += m.getAtomCount();
    }
    return threshold + reactantTotalAtomCount > productAtomCount;
  }

  static Boolean isValidReactionResultByReaction(String sml, Molecule[] reactants, Molecule rxn) {
    try {
      Reactor reactor = new Reactor();
      reactor.setIgnoreRules(Reactor.IGNORE_REACTIVITY);
      reactor.setReaction(rxn);
      reactor.setReactants(reactants);
      Molecule[] p = null;
      while ((p = reactor.react()) != null) {
        try {
          if (MolExporter.exportToFormat(p[0], "smiles:u").equals(sml)) {
            return true;
          }
        } catch (IOException ioe) {
          return false;
        }
      }
      return false;
    } catch (ReactionException re) {
      return false;
    }
  }

  public void depthFirstProofNumberSearch() throws ReactionException, IOException, InterruptedException {
    System.out.println("Searching...");
    long s, e = 0;
    this.startTime = System.currentTimeMillis();
    ArrayList<Node> init = new ArrayList<Node>();
    multipleIterativeOrDeepening(this.root, 0, init);
    e = System.currentTimeMillis();
    System.out.println("end");
    Visualize v = new Visualize();
    v.printProofTree(this.root);
    v.printAll();
    for (Map.Entry me : this.availableCompoundList.entrySet()) {
      System.out.println(me.getKey() + "," + me.getValue());
    }
    int idx = selectNode(this.root);
    if (idx != -1) {
      System.out.println(this.root.children.get(idx).id);
    }
    System.out.println("Terminate!!!");
  }
  
  
  @SuppressWarnings("deprecation")
  public void expandNode(Node node,long depth) throws ReactionException, IOException, InterruptedException {
    int index = 0;
    System.out.println("expand...");
    if (node.isExpanded) {
      System.out.println("return because of expanded flag");
      return;
    }
    if (node.smiles.equals("Clc1ccc(cc1)C(N1CCNCC1)c1ccccc1")) {
      System.out.println("She left him just alive enough for me to see him die.");
      TimeUnit.SECONDS.sleep(1);
    }
    for (int i = 0; i < this.reactionData.size(); ++i) {
      // retrosynthesis
      Molecule rxn = this.reactionData.get(i);
      Molecule[] reactants = new Molecule[1];
      reactants[0] = node.mol;
      Molecule[] products = null;

      try {
        Reactor reactor = new Reactor();
        reactor.setIgnoreRules(Reactor.IGNORE_REACTIVITY);
        reactor.setReaction(rxn);
        reactor.setReverse(true);
        reactor.setReactants(reactants);
        while ((products = reactor.react()) != null) {
          int k = products.length;

          boolean isValidResult = true;
          if (!isValidReactionResultByReaction(node.smiles, products, rxn)) {
            System.err.println("Invalid reaction result:");
            System.err.println("target:" + node.smiles);
            System.err.println("reaction:" + MolExporter.exportToFormat(rxn, "smarts"));
            for (Molecule p : products) {
              System.err.println("reactant:" + MolExporter.exportToFormat(p, "smiles:u"));
            }
            isValidResult = false;
            continue;
          } else {
            System.out.printf("valid result %s\n", node.smiles);
            System.out.printf("valid reaction %s\n", MolExporter.exportToFormat(rxn, "smarts"));
            for (Molecule r : products) {
              System.out.printf("valid reactant %s\n", MolExporter.exportToFormat(r, "smiles:u"));
            }
          }
          boolean isValidCompoundSet = true;
          for (int j = 0; j < k; ++j) {
            if (!isValidValence(products[j])) {
              isValidCompoundSet = false;
              try {
                System.err.println(MolExporter.exportToFormat(products[j], "smiles:u") + " is invalid.");
              } catch (IOException e) {
                System.err.println("Molecule could not be exported into SMILES format.");
              }
              break;
            }
          }
          
          if (isValidResult && isValidCompoundSet && k >= 1) {
            Node reactionNode = new Node(rxn,"reaction",MolExporter.exportToFormat(rxn, "smarts"),node.depth+1,node);
            int prevId = this.id;
            this.id += 1;
            reactionNode.setId(this.id);
            boolean nodeDelete = false;
            for (int t = 0; t < k; ++t) {
              try {
                String key = MolExporter.exportToFormat(products[t], "smiles:u");
                long e, s = 0;
                s = System.currentTimeMillis();
                if (!nodeHashMap.containsKey(key)) {
                  Node moleculeNode = new Node(products[t],"molecule",key, node.depth+2,node);
                  nodeHashMap.put(key, moleculeNode);
                  this.id += 1;
                  moleculeNode.setId(this.id);
                  reactionNode.children.add(moleculeNode);
                } else {
                  // Node joined!!!
                  Node prev = nodeHashMap.get(key);
                  prev.depth = Math.min(node.depth+2, prev.depth);
                  reactionNode.children.add(prev);
                }
                e = System.currentTimeMillis();
              } catch (MolExportException e) {
                nodeDelete = true;
              }
            }
            if (!nodeDelete) {
              node.children.add(reactionNode);
              index += 1;
            } else {
              reactionNode = null;
            }
          }
        }
      } catch (ReactionException e) {
        System.err.println("This reaction is not suitable for retrosynthesis.");
        continue;
      }
    }
    node.setPhi(deltaMin(node, false));
    node.setDelta(phiSum(node, false));
    node.updateDeltaWithCost();
    node.updatePhiWithCost();
    System.out.println("Done.");
  }
    
  public void multipleIterativeOrDeepening(Node node,int depth, ArrayList<Node> path) throws ReactionException, IOException, InterruptedException {
    // ノードの展開
    System.out.printf("OR deepening for id:%d\n", node.id);
    if (node.nodeType.equals("reaction")) {
      System.out.println("Illegal node type for multipleIterativeOrDeepening.");
      return;
    }
    this.endTime = System.currentTimeMillis();
    System.out.printf("Current computation time:%d\n", this.endTime - this.startTime);
    if ((this.endTime - this.startTime > this.timeThreshold) && this.stopByTime) {
      System.out.printf("Reached time threshold in OR:endTime %d, startTime %d, computaionTime %d, timeThreshold %d\n", this.endTime, this.startTime, this.endTime - this.startTime, this.timeThreshold);
      return;
    }
    for (int i = 0; i < path.size(); ++i) {
      node.pathFromRoot.add(path.get(i));
    }
    node.pathFromRoot.add(node);
    boolean isAvailable = false;
    if(!node.isRoot && !node.isChecked) {
      isAvailable = commerciallyAvailableCompound(node.mol);
      node.isChecked = true;
    }
    Visualize v = new Visualize();
    System.out.println("checking...");
    
    if (isAvailable) {
      System.out.println(node.id + " " + node.smiles +  " is available.");
      // if commercially available, set proof number to 0
      this.availableCompoundList.put(node.id, node.smiles);
      node.setPhi(0);
      node.setDelta(Integer.MAX_VALUE / 2);
      node.isProven = true;
      node.isAvailable = true;
      //node.ignore = true;
      moleculeTableSave(node);
      node.pathFromRoot.clear();
      node.routeNum = 1;
      System.out.printf("route num of %d is updated to %d because it is available.\n", node.id, node.routeNum);
      return;
    } else {
      System.out.println(node.id + " " + node.smiles +  " is not available.");
    }
    if (this.maxDepth <= node.depth) {
      // if not available and exceed maxdepth, set disproof number to 0
      System.out.println(node.id + " " + node.smiles +  " is too deep.");
      if (node.getProofNumber() != 0) {
        node.setPhi(Integer.MAX_VALUE / 2);
        node.setDelta(0);
        node.isProven = false;
        node.isDisproven = true;
        //node.ignore = true;
        moleculeTableSave(node);
      }
      node.pathFromRoot.clear();
      return;
    }
    
    moleculeTableLookUp(node);
    
    
    if (!node.isExpanded) {
      // 子ノードのindex
      expandNode(node,depth);
      node.isExpanded = true;
    }
    if (node.getPhiWithCost() >= node.phiThreshold || node.getDeltaWithCost() >= node.deltaThreshold) {
      node.phiThreshold = node.getPhiWithCost();
      node.deltaThreshold = node.getDeltaWithCost();
      node.pathFromRoot.clear();
      return;
    }
    if (node.getReachableChildrenSize() == 0) {
      System.out.println("finalUpdate");
      node.finalUpdate();
      node.pathFromRoot.clear();
      printNodeInfo(node);
      System.out.println("End of OrDeepening.");
      return;
    }
    System.out.printf("id:%d, phiThreshold:%d, deltaThreshold:%d, deltaMin(node, true):%d phiSum(node, true):%d\n", node.id, node.phiThreshold, node.deltaThreshold, deltaMin(node, true), phiSum(node, true));
    // cost付きのdeltaMin, phiSumが閾値を下回っている限りwhile文をまわす
    while (node.phiThreshold > deltaMin(node, true) && node.deltaThreshold > phiSum(node, true)) {
      this.phiChild = 0;
      this.deltaSecond = 0;
      int bestIndex = selectNode(node);
      node.lastSelectNodeIndex = bestIndex;
      if (bestIndex != -1) {
        System.out.printf("id:%d, selected id:%d, phiChild:%d, deltaSecond:%d\n", node.id, node.children.get(bestIndex).id, this.phiChild, this.deltaSecond);
        printNodeInfo(node.children.get(bestIndex));
      } else {
        System.out.printf("children num without ignore in OR %d\n", node.getReachableChildrenSize());
        System.out.printf("%d will be ignored because all children are ignored", node.id);
        break;
      }
      
      // 各閾値を更新
      node.children.get(bestIndex).phiThreshold = node.deltaThreshold + this.phiChild - phiSum(node, true); // これでいいのか
      System.out.printf("phiThreshold:%d, prevPth:%d\n", node.children.get(bestIndex).phiThreshold, node.children.get(bestIndex).prevPth);
      // if the value is not changed, infinite loop occurs
      if (node.children.get(bestIndex).phiThreshold == node.children.get(bestIndex).prevPth) {
        node.children.get(bestIndex).phiThreshold += 1;
        node.children.get(bestIndex).prevPth = node.children.get(bestIndex).phiThreshold; 
      } else if (node.children.get(bestIndex).phiThreshold <= node.children.get(bestIndex).prevPth) {
        System.out.printf("phiThreshold:%d %d\n", node.children.get(bestIndex).phiThreshold, node.children.get(bestIndex).phiThreshold+1);        
        node.children.get(bestIndex).phiThreshold = node.children.get(bestIndex).prevPth + 1;
        node.children.get(bestIndex).prevPth = node.children.get(bestIndex).phiThreshold;
      } else {
        node.children.get(bestIndex).prevPth = node.children.get(bestIndex).phiThreshold;
      }
      node.children.get(bestIndex).deltaThreshold = Math.min(node.phiThreshold, this.deltaSecond+1);
      System.out.printf("deltaThreshold:%d, prevDth:%d\n", node.children.get(bestIndex).deltaThreshold, node.children.get(bestIndex).prevDth);
      if (node.children.get(bestIndex).deltaThreshold == node.children.get(bestIndex).prevDth) {
        node.children.get(bestIndex).deltaThreshold += 1;
        node.children.get(bestIndex).prevDth = node.children.get(bestIndex).deltaThreshold; 
      } else if (node.children.get(bestIndex).deltaThreshold <= node.children.get(bestIndex).prevDth) {
        System.out.printf("deltaThreshold:%d %d\n", node.children.get(bestIndex).deltaThreshold, node.children.get(bestIndex).deltaThreshold+1);
        node.children.get(bestIndex).deltaThreshold = node.children.get(bestIndex).prevDth + 1;
        node.children.get(bestIndex).prevDth = node.children.get(bestIndex).deltaThreshold;
      } else {
        node.children.get(bestIndex).prevDth = node.children.get(bestIndex).deltaThreshold; 
      }
      node.setPhi(deltaMin(node, false));
      // System.out.printf("Before And Deepening, phi of %d is updated to %d\n", node.id, node.getPhi());
      node.setDelta(phiSum(node, false));
      // System.out.printf("Before And Deepening, delta of %d is updated to %d\n", node.id, node.getDelta());
      node.updatePhiWithCost();
      node.updateDeltaWithCost();
      //v.printWholeTree(this.root);
      multipleIterativeAndDeepening(node.children.get(bestIndex),node.depth, node.pathFromRoot);
      //v.printWholeTree(this.root);
      if (deltaMin(this.root, false) == 0) {
        this.endTime = System.currentTimeMillis();
        System.out.printf("Current computation time:%d\n", this.endTime - this.startTime);
        if ((this.endTime - this.startTime > this.timeThreshold) && this.stopByTime) {
          System.out.printf("Reached time threshold in OR:endTime %d, startTime %d, computaionTime %d, timeThreshold %d\n", this.endTime, this.startTime, this.endTime - this.startTime, this.timeThreshold);
          return;
        }
        node.updateRouteNum();
        System.out.println("trace proof tree");
        traceProofTree(this.root);
        //v.printWholeTree(node);
        //v.printWholeTree(this.root);
        this.foundLastNode = false;
        this.pathFromRoot.clear();
        v.printProofTree(this.root);
        for (Map.Entry me : this.availableCompoundList.entrySet()) {
          System.out.println(me.getKey() + "," + me.getValue());
        }
        this.endTime = System.currentTimeMillis();
        System.out.printf("found proof tree:%d,%d\n", this.endTime - this.startTime, node.routeNum);
        if (this.routeNumThreshold > 0 && this.root.routeNum > this.routeNumThreshold) {
          System.out.printf("reached routeNumThreshold %d with %d\n", this.routeNumThreshold, this.root.routeNum);
          return;
        }
      }
      node.setPhi(deltaMin(node, false));
      // System.out.printf("phi of %d is updated to %d\n", node.id, node.getPhi());
      node.setDelta(phiSum(node, false));
      // System.out.printf("delta of %d is updated to %d\n", node.id, node.getDelta());
      node.updatePhiWithCost();
      node.updateDeltaWithCost();
      if (node.getProofNumber() == 0 || node.getDisProofNumber() == 0) {
        System.out.printf("In OrDeepening, break with %d\n", node.id);
        break;
      }
      this.endTime = System.currentTimeMillis();
      System.out.printf("Current computation time:%d\n", this.endTime - this.startTime);
      if ((this.endTime - this.startTime > this.timeThreshold) && this.stopByTime) {
        System.out.printf("Reached time threshold in OR:endTime %d, startTime %d, computaionTime %d, timeThreshold %d\n", this.endTime, this.startTime, this.endTime - this.startTime, this.timeThreshold);
        return;
      }
    }
    node.pathFromRoot.clear();
    // printNodeInfo(node);
    System.out.println("End of OrDeepening");
  }

  public void multipleIterativeAndDeepening(Node node,int depth, ArrayList<Node> path) throws ReactionException, IOException, InterruptedException {
    if (node.nodeType.equals("molecule")) {
      System.out.println("Illegal node type for multipleIterativeAndDeepening.");
      return;
    }
    this.endTime = System.currentTimeMillis();
    System.out.printf("Current computation time:%d\n", this.endTime - this.startTime);
    if ((this.endTime - this.startTime > this.timeThreshold) && this.stopByTime) {
      System.out.printf("Reached time threshold in AND:endTime %d, startTime %d, computaionTime %d, timeThreshold %d\n", this.endTime, this.startTime, this.endTime - this.startTime, this.timeThreshold);
      return;
    }
    System.out.printf("AND deepening for id:%d\n", node.id);
    assert node.nodeType.equals("reaction");
    
    if (!node.isExpanded) {
      node.isExpanded = true;
      node.setPhi(deltaMin(node, false));
      node.setDelta(phiSum(node, false));
      if (node.getDisProofNumber() == 0) {
        node.ignore = true;
        System.out.printf("$$$id:%d is flagged as ignored\n", node.id);
        return;
      }
      node.updatePhiWithCost();
      node.updateDeltaWithCost();
    }
    for (int i = 0; i < path.size(); ++i) {
      node.pathFromRoot.add(path.get(i));
    }
    node.pathFromRoot.add(node);
    if (node.getPhiWithCost() >= node.phiThreshold || node.getDeltaWithCost() >= node.deltaThreshold) {
      System.out.printf("id:%d, phiWithCost:%d, phiThreshold:%d, deltaWithCost:%d, deltaThreshold:%d\n", node.id, node.getPhiWithCost(), node.phiThreshold, node.getDeltaWithCost(), node.deltaThreshold);
      node.phiThreshold = node.getPhiWithCost();
      node.deltaThreshold = node.getDeltaWithCost();
      node.pathFromRoot.clear();
      return;
    }
    assert node.getPhi() >= 0 : "Phi below 0";
    if (node.getReachableChildrenSize() == 0) {
      node.finalUpdate();
      node.pathFromRoot.clear();
      // printNodeInfo(node);
      return;
    }
    // System.out.printf("id:%d, phiThreshold:%d, deltaThreshold:%d, deltaMin(node, true):%d phiSum(node, true):%d\n", node.id, node.phiThreshold, node.deltaThreshold, deltaMin(node, true), phiSum(node, true));
    while (node.phiThreshold > deltaMin(node, true) && node.deltaThreshold > phiSum(node, true)) {
      int bestIndex = selectNode(node);
      node.lastSelectNodeIndex = bestIndex;
      
      if (bestIndex != -1) {
        System.out.printf("id:%d, selected id:%d, phiChild:%d, deltaSecond:%d\n", node.id, node.children.get(bestIndex).id, this.phiChild, this.deltaSecond);
      } else {
        System.out.printf("%d will be ignored because all children are proven, disproven, or ignored.\n", node.id);
        break;
      }
      node.children.get(bestIndex).phiThreshold = node.deltaThreshold + this.phiChild - phiSum(node, true);
      System.out.printf("phiThreshold:%d, prevPth:%d\n", node.children.get(bestIndex).phiThreshold, node.children.get(bestIndex).prevPth);
      if (node.children.get(bestIndex).phiThreshold == node.children.get(bestIndex).prevPth) {
        node.children.get(bestIndex).phiThreshold += 1;
        node.children.get(bestIndex).prevPth = node.children.get(bestIndex).phiThreshold; 
      } else if (node.children.get(bestIndex).phiThreshold <= node.children.get(bestIndex).prevPth) {
        System.out.printf("phiThreshold:%d %d\n", node.children.get(bestIndex).phiThreshold, node.children.get(bestIndex).phiThreshold+1);        
        node.children.get(bestIndex).phiThreshold = node.children.get(bestIndex).prevPth + 1;
        node.children.get(bestIndex).prevPth = node.children.get(bestIndex).phiThreshold;
      } else {
        node.children.get(bestIndex).prevPth = node.children.get(bestIndex).phiThreshold;
      }
      node.children.get(bestIndex).deltaThreshold = Math.min(node.phiThreshold,this.deltaSecond+1);
      System.out.printf("deltaThreshold:%d, prevDth:%d\n", node.children.get(bestIndex).deltaThreshold, node.children.get(bestIndex).prevDth);
      if (node.children.get(bestIndex).deltaThreshold == node.children.get(bestIndex).prevDth) {
        node.children.get(bestIndex).deltaThreshold += 1;
        node.children.get(bestIndex).prevDth = node.children.get(bestIndex).deltaThreshold; 
      } else if (node.children.get(bestIndex).deltaThreshold <= node.children.get(bestIndex).prevDth) {
        System.out.printf("deltaThreshold:%d %d\n", node.children.get(bestIndex).deltaThreshold, node.children.get(bestIndex).deltaThreshold+1);
        node.children.get(bestIndex).deltaThreshold = node.children.get(bestIndex).prevDth + 1;
        node.children.get(bestIndex).prevDth = node.children.get(bestIndex).deltaThreshold;
      } else {
        node.children.get(bestIndex).prevDth = node.children.get(bestIndex).deltaThreshold; 
      }
      node.setPhi(deltaMin(node, false));
      // System.out.printf("Before Or Deepening, phi of %d is updated to %d\n", node.id, node.getPhi());
      node.setDelta(phiSum(node, false));
      // System.out.printf("Before Or Deepening, delta of %d is updated to %d\n", node.id, node.getDelta());
      node.updatePhiWithCost();
      node.updateDeltaWithCost();
      Visualize v = new Visualize();
      //v.printWholeTree(this.root);
      multipleIterativeOrDeepening(node.children.get(bestIndex),node.depth, node.pathFromRoot);
      //v.printWholeTree(this.root);
      System.out.println("ginyaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
      // OrDeepeningから返ってきたら証明数、反証数を更新
      node.setPhi(deltaMin(node, false));
      System.out.printf("phi of %d is updated to %d\n", node.id, node.getPhi());
      node.setDelta(phiSum(node, false));
      System.out.printf("delta of %d is updated to %d\n", node.id, node.getDelta());
      node.updatePhiWithCost();
      node.updateDeltaWithCost();
      //v.printWholeTree(this.root);
      if (node.getProofNumber() == 0) {
        System.out.printf("In And Deepening, break with %d because it is proven.\n", node.id);
        break;
      } else if (node.getDisProofNumber() == 0) {
        //反証されたANDノードは以降の探索では無視するのでignoreフラグを立てる
        System.out.printf("In And Deepening, break with %d because it is disproven.\n", node.id);
        node.ignore = true;
        break;
      }
      this.endTime = System.currentTimeMillis();
      System.out.printf("Current computation time:%d\n", this.endTime - this.startTime);
      if ((this.endTime - this.startTime > this.timeThreshold) && this.stopByTime) {
        System.out.printf("Reached time threshold in AND:endTime %d, startTime %d, computaionTime %d, timeThreshold %d\n", this.endTime, this.startTime, this.endTime - this.startTime, this.timeThreshold);
        return;
      }
    }
    node.pathFromRoot.clear();
    // printNodeInfo(node);
    // System.out.println("End of AndDeepening.");
  }
    
  public int selectNode(Node node) throws IOException {
    int bestIndex = -1;
    int deltaChild = Integer.MAX_VALUE / 2;
    this.deltaSecond = Integer.MAX_VALUE / 2;
    this.phiChild = Integer.MAX_VALUE / 2;
    //Visualize v = new Visualize();
    //v.treeToDot(node, false);
    //v.printDot();
    System.out.println("start selecting");
    for (int i = 0; i < node.children.size(); ++i) {
      Node n = node.children.get(i);
      if (node.nodeType.equals("molecule")) {
        if (deltaMin(n, false) == 0) {
          n.setPhi(deltaMin(n, false));
          n.setDelta(phiSum(n, false));
          n.ignore = true;
        }
      }
      //printNodeInfo(node.children.get(i));
      if (node.children.get(i).ignore) {
        //System.out.printf("In selectNode for %d, skip %d because it is flagged as ignore.\n", node.id, node.children.get(i).id);
        continue;
      }
      if (node.children.get(i).getProofNumber() == 0) {
        //System.out.printf("In selectNode for %d, skip %d because it is proven.\n", node.id, node.children.get(i).id);
        continue;
      }
      if (node.children.get(i).getDisProofNumber() == 0) {
        //System.out.printf("In selectNode for %d, skip %d because it is disproven.\n", node.id, node.children.get(i).id);
        continue;
      }
      if (node.nodeType.equals("molecule") && node.isMakingCycle(node.children.get(i))) {
        //System.out.printf("%d is making cycle! with %d\n", node.children.get(i).id, node.id);
        continue;
      }
      //System.out.println(i);
      if (node.nodeType.equals("molecule")) {
        moleculeTableLookUp(node);
      }
      int dc = node.children.get(i).getDeltaWithCost();
      int pc = node.children.get(i).getPhiWithCost();
      if (dc < deltaChild) {
        bestIndex = i;
        this.deltaSecond = deltaChild;
        deltaChild = dc;
        this.phiChild = pc;
      } else if (dc < this.deltaSecond) {
        this.deltaSecond = dc;
      }
    }
    return bestIndex;
  }
  
  public void traceProofTree(Node node) throws IOException {
    this.pathFromRoot.add(node);
    if (node.lastSelectNodeIndex != -1) {
      Node next = node.children.get(node.lastSelectNodeIndex);
      boolean nextRec = true;
      for (Node tmp : this.pathFromRoot) {
        if (tmp.id == next.id) {
          Node lastOR = this.pathFromRoot.get(this.pathFromRoot.size()-1);
          lastOR.children.remove(node);
          node = null;
          nextRec = false;
        }
      }
      if (nextRec) {
        traceProofTree(next);
      } else {
        return;
      }
    } else {
      // finally
      for (int i = 0; i < this.pathFromRoot.size(); ++i) {
        System.out.println(this.pathFromRoot.get(i).id);
      }
      int pathLen = this.pathFromRoot.size();
      Node lastFoundNode = this.pathFromRoot.get(pathLen-1);
      boolean flagProofTree = false;
      if (lastFoundNode.getProofNumber() == 0) {
        flagProofTree = true;
      }
      
      Node lastFoundAndNode = null;
      
      // when joined to proven node which is not flagged as part of proof tree
      for (int i = this.pathFromRoot.size() - 1; i >= 0; --i) {
        Node tmp = this.pathFromRoot.get(i);
        if (tmp.nodeType.equals("reaction") && this.foundLastNode && !tmp.isPartOfProofTree) {
          for (int j = 0; j < tmp.children.size(); ++j) {
            if (j == tmp.lastSelectNodeIndex) {
              continue;
            }
            List<Node> sideRoot = new ArrayList<Node>();
            sideRoot.add(tmp);
            sideWalk(tmp.children.get(j), sideRoot);
          }
        }
        if (tmp.nodeType.equals("reaction") && !this.foundLastNode) {
          lastFoundAndNode = tmp;
          this.foundLastNode = true;
        }
      }
      lastFoundAndNode.ignore = true;  
      
      for (int i = 0; i < this.pathFromRoot.size(); ++i) {
        // System.out.printf("$$$id:%d is flagged as part of proof tree\n", this.pathFromRoot.get(i).id);
        this.pathFromRoot.get(i).isPartOfProofTree = true;
      }
      
      if (flagProofTree) {
        // System.out.printf("$$$id:%d is flagged as part of proof tree\n", lastFoundAndNode.id);
        lastFoundAndNode.isPartOfProofTree = true;
        for (int i = 0; i < lastFoundAndNode.children.size(); ++i) {
          // System.out.printf("id:%d is flagged as part of proof tree\n", lastFoundAndNode.children.get(i).id);
          lastFoundAndNode.children.get(i).isPartOfProofTree = true;
        }
      } else if (!flagProofTree) {
        int index = lastFoundAndNode.lastSelectNodeIndex;
        if (index != -1) {
          System.out.printf("$$$id:%d is flagged as not  part of proof tree\n", lastFoundAndNode.children.get(index).id);
        }
      }
      System.out.printf("%d is flagged as ignore\n",lastFoundAndNode.id);
      System.out.println("Found last node");
    }
    if (node.children.size() > 0) {
        node.setPhi(deltaMin(node, false));
        node.setDelta(phiSum(node, false));
        if ((node.getPhi() == 0 || node.getDelta() == 0) && node.nodeType.equals("reaction")) {
          node.ignore = true;
          // System.out.printf("$$$id:%d is flagged as ignored\n", node.id);
        }
        node.updatePhiWithCost();
        node.updateDeltaWithCost();
    }
    node.updateRouteNum();
  }
  
  public void sideWalk(Node node, List<Node> sideRoot) throws IOException {
    // System.out.printf("side walk for id:%d.\n", node.id);
    // System.out.printf("$$$id:%d is flagged as part of proof tree\n", node.id);
    node.isPartOfProofTree = true;
    if (node.lastSelectNodeIndex != -1) {
      Node next = node.children.get(node.lastSelectNodeIndex);
      sideRoot.add(next);
      if (node.getProofNumber() == 0) {
        // System.out.printf("$$$id:%d is flagged as part of proof tree in sideWalk\n", node.id);
        next.isPartOfProofTree = true;
      }
      sideWalk(next, sideRoot);
    } else {
      Node lastFoundAndNode = null;
      for (int i = sideRoot.size() - 1; i >= 0; --i) {
        if (sideRoot.get(i).nodeType.equals("reaction")) {
          lastFoundAndNode = sideRoot.get(i);
          break;
        }
      }
      // lastFoundAndNode.ignore = true;
      System.out.printf("%d is flagged as ignore\n",lastFoundAndNode.id);
      //lastFoundAndNode.setPhi(0);
      //lastFoundAndNode.setDelta(Integer.MAX_VALUE / 2);
    }
    if (node.children.size() > 0) {
      node.setPhi(deltaMin(node, false));
      System.out.printf("phi of %d has been updated to %d in sideWalk\n", node.id, node.getPhi());
      node.setDelta(phiSum(node, false));
      System.out.printf("delta of %d has been updated to %d in sideWalk\n", node.id, node.getDelta());
      node.updatePhiWithCost();
      node.updateDeltaWithCost();
    }
    node.updateRouteNum();
  }

  public static void moleculeTableLookUp(Node node) {
    String key = node.smiles;
    if (nodeHashMap.containsKey(key)) {
      node.setPhi(nodeHashMap.get(key).getPhi());
      node.setDelta(nodeHashMap.get(key).getDelta());
    } else {
      node.setPhi(1);
      node.setDelta(1);
    }
  }
  
  public void moleculeTableSave(Node node) {
    String key = node.smiles;
    if (!nodeHashMap.containsKey(key)) {
      nodeHashMap.put(key, node);
    } else {
      Node tmp = nodeHashMap.get(key);
      tmp.setPhi(node.getPhi());
      tmp.setDelta(node.getDelta());
      tmp.isProven = node.isProven;
      tmp.ignore = node.ignore;
    }
  }
  
  
  public int deltaMin(Node node, boolean useCost) throws IOException {
    if (useCost) {
      assert node.getReachableChildrenSize() > 0;
    }
    int d = Integer.MAX_VALUE / 2;
    for (int i = 0; i < node.children.size(); ++i) {
      Node n = node.children.get(i);
      if (n.ignore) {
        //System.out.printf("In deltaMin for %d, %d is skipped because it is flagged as ignored.\n", node.id, n.id);
        continue;
      }
      if (useCost && node.nodeType.equals("molecule") && n.getDelta() == 0) {
        continue;
      }
      
      if (useCost && (n.getDelta() == 0 || n.getPhi() ==0)) {
        continue;
      }
      if (n.nodeType.equals("molecule")) {
        moleculeTableLookUp(n);
      }
      
      int tmp_d = 0;
      if (useCost) {
        tmp_d = n.getDeltaWithCost();
      } else {
        tmp_d = n.getDelta();
      }
      d = Math.min(d, tmp_d);
    }
    if (useCost) {
      System.out.printf("deltaMin with cost for %d is %d\n", node.id, d);
    } else {
      System.out.printf("deltaMin for %d is %d\n", node.id, d);
    }
    return d;
  }
  
  public int phiSum(Node node, boolean useCost) {
    if (useCost) {
      assert node.getReachableChildrenSize() > 0;
    }
    int p = 0;
    for (int i = 0; i < node.children.size(); ++i) {
      Node n = node.children.get(i);
      if (n.ignore) {
        continue;
      }
      if (useCost && node.nodeType.equals("molecule") && n.getPhi() > Integer.MAX_VALUE / 3) {
        continue;
      }
      if (useCost && (n.getDelta() == 0 || n.getPhi() ==0)) {
        continue;
      }
      if (useCost) {
        p += n.getPhiWithCost();
      } else {
        p += n.getPhi();
      }
      if (p > Integer.MAX_VALUE / 3) {
        break;
      }
    }
    if (useCost) {
      System.out.printf("phiSum with cost for %d is %d\n", node.id, p);
    } else {
      System.out.printf("phiSum for %d is %d\n", node.id, p);
    }
    return p;
  }
  
  public void printNodeInfo(Node node) throws IOException {
    System.out.println("Node information");
    System.out.println("Node id");
    System.out.println(node.id);
    System.out.println(node.nodeType);
    if (node.nodeType.equals("molecule")) {
      System.out.println(node.smiles);
    } else {
      System.out.println(node.reaction);
    }
    System.out.printf("depth:%d\n", node.depth);
    if (node.isExpanded) {
      System.out.println("Already expanded.");
    } else {
      System.out.println("Not expanded yet.");
    }
    if (node.ignore) {
      System.out.println("Must be ignored");
    } else {
      System.out.println("To be searched");
    }
    if (node.isProven) {
      System.out.println("Proven.");
    }
    if (node.isDisproven) {
      System.out.println("Disproven.");
    }
    if (node.isPartOfProofTree) {
      System.out.printf("%d is part of proof tree\n",node.id);
    }
    System.out.printf("children:%d\n",node.children.size());
    if (node.nodeType.equals("reaction")) {
      for (int i = 0; i < node.children.size(); ++i) {
        System.out.println(node.children.get(i).smiles);
      }
    }
    System.out.printf("phi:%d\n",node.getPhi());
    System.out.printf("phiWithCost:%d\n", node.getPhiWithCost());
    System.out.printf("deltaWithCost:%d\n", node.getDeltaWithCost());
    System.out.printf("phi th:%d\n",node.phiThreshold);
    System.out.printf("delta:%d\n",node.getDelta());
    System.out.printf("delta th:%d\n",node.deltaThreshold);
    //System.out.printf("cost:%d\n",node.getCost());
    //System.out.printf("depth:%d\n",node.depth);
    
    if (node.lastSelectNodeIndex != -1) {
      System.out.printf("last select node index:%d\n", node.lastSelectNodeIndex);
    }
    if (node.isRoot) {
      System.out.println("root");
    }
  }
}
