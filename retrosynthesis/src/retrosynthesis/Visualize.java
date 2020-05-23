package retrosynthesis;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class Visualize {
  public String g;
  public Set<Edge> edgeSet = new HashSet<Edge>();
  public long nodeCount;
  public Map<Integer, String> allCompounds = new HashMap<Integer, String>();
  public Map<Integer, String> allReactions = new HashMap<Integer, String>();
  Visualize() {}
  
  public void printProofTree(Node node) {
    this.g = "digraph g{";
    this.edgeSet = new HashSet<Edge>();
    this.nodeCount = 0;
    proofTreeToDot(node);
    printDot();
    printAll();
    System.out.printf("node count:%d\n", this.nodeCount);
  }
  
  public void printWholeTree(Node node, boolean out) {
    this.g = "$.digraph g{";
    this.edgeSet = new HashSet<Edge>();
    this.nodeCount = 0;
    treeToDot(node, false);
    if (out) {
      printDot();
    }
    System.out.printf("node count:%d\n", this.nodeCount);
  }
  
  // visualize by depth-first search
  public void treeToDot(Node n, boolean proofTreeMode) {
    this.nodeCount += 1;
    String parentId = Integer.toString(n.id);
    g += " \"" + parentId + "\"";
    String t = "";
    if (n.isPartOfProofTree) {
      t += "pt";
    }
    String t1 = "";
    if (n.ignore) {
      t1 = "ignored";
    } if (n.getProofNumber() != 0 && n.getDisProofNumber() != 0) {
      t1 = "unknown";
    }
    if (n.nodeType == "molecule") {
      String tmp = "";
      if (n.lastSelectNodeIndex != -1) {
        tmp = " [shape=circle, label=\"" + "id:" + parentId + "\\n"
                                         + "pn(phi):" + Integer.toString(n.getPhi()) + "\\n" + "dn(delta):" + Integer.toString(n.getDelta()) + "\\n"
                                         + "phi+cost:" + Long.toString(n.getPhiWithCost()) + "\\n" + "delta+cost:" + Long.toString(n.getDeltaWithCost()) + "\\n"
                                         + "pth:" + Integer.toString(n.phiThreshold) + "\\n" + "dth:" + Integer.toString(n.deltaThreshold) + "\\n"
                                         + "cost:" + Integer.toString(n.cost) + "\\n" + "depth:" + Integer.toString(n.depth) + "\\n"
                                         + "next:" + Integer.toString(n.children.get(n.lastSelectNodeIndex).id) + "\\n"
                                         + t + "\\n" + t1 +  "\"]";
      } else {
        tmp = " [shape=circle, label=\"" + "id:" + parentId + "\\n"
                                         + "pn(phi):" + Long.toString(n.getPhi()) + "\\n" + "dn(delta):" + Long.toString(n.getDelta()) + "\\n"
                                         + "phi+cost:" + Long.toString(n.getPhiWithCost()) + "\\n" + "delta+cost:" + Long.toString(n.getDeltaWithCost()) + "\\n"
                                         + "pth:" + Long.toString(n.phiThreshold) + "\\n" + "dth:" + Long.toString(n.deltaThreshold) + "\\n"
                                         + "cost:" + Integer.toString(n.cost) + "\\n" + "depth:" + Integer.toString(n.depth) + "\\n"
                                         + t + "\\n" + t1 + "\"]";
      }
      g += tmp;
    } else {
      String tmp = "";
      if (n.lastSelectNodeIndex != -1) {
        tmp = " [shape=box, label=\"" + "id:" + parentId + "\\n" 
                                      + "pn(delta):" + Long.toString(n.getDelta()) + "\\n" + "dn(phi):" + Long.toString(n.getPhi()) + "\\n"
                                      + "phi+cost:" + Long.toString(n.getPhiWithCost()) + "\\n" + "delta+cost:" + Long.toString(n.getDeltaWithCost()) + "\\n"
                                      + "pth:" + Long.toString(n.phiThreshold) + "\\n" + "dth:" + Long.toString(n.deltaThreshold) + "\\n" 
                                      + "cost:" + Integer.toString(n.cost) + "\\n" +  "depth:" + Integer.toString(n.depth) + "\\n" 
                                      + "next:" + Integer.toString(n.children.get(n.lastSelectNodeIndex).id) + "\\n"
                                      + t + "\\n" + t1 + "\"]";
      } else {
        tmp = " [shape=box, label=\"" + "id:" + parentId + "\\n"
                                      + "pn(delta):" + Long.toString(n.getDelta()) + "\\n" + "dn(phi):" + Long.toString(n.getPhi()) + "\\n"
                                      + "phi+cost:" + Long.toString(n.getPhiWithCost()) + "\\n" + "delta+cost:" + Long.toString(n.getDeltaWithCost()) + "\\n"
                                      + "pth:" + Long.toString(n.phiThreshold) + "\\n" + "dth:" + Long.toString(n.deltaThreshold) + "\\n"
                                      + "cost:" + Integer.toString(n.cost) + "\\n" + "depth:" + Integer.toString(n.depth) + "\\n"
                                      + t + "\\n" + t1 + "\"]";
      }
      g += tmp;
    }
    for (int i = 0; i < n.children.size(); ++i) {
      String id = Integer.toString(n.children.get(i).id);
      Edge e = new Edge(n.id,n.children.get(i).id);
      
      if (isExist(e)) {
        continue;
      } else {
        edgeSet.add(e);
      }
      //String tmp = "\"" + parentId + '"' + " -- " + "\"" + id + "\"" + " ";
      String tmp = "\"" + parentId + '"' + " -> " + "\"" + id + "\"" + " ";
      g += tmp;
      treeToDot(n.children.get(i), proofTreeMode);
    }
  }
  
  public void proofTreeToDot(Node n) {
    this.nodeCount += 1;
    String parentId = Integer.toString(n.id);
    String pathNum = Long.toString(n.routeNum);
    g += " \"" + parentId + "\"";
    if (n.nodeType.equals("molecule")) {
      this.allCompounds.put(n.id, n.smiles);
      String tmp = " [shape=circle, label=\"" + parentId + "\\n" 
                 //+ " smiles=\"" + n.smiles + "\\n" 
                 + pathNum + "\"]";
      //String tmp = " [shape=circle] "+ parentId;
      g += tmp;
    } else {
      this.allReactions.put(n.id, n.reaction);
      String tmp = " [shape=box, label=\"" + parentId + "\\n" 
                 //+ " smiles=\"" + n.reaction + "\\n" 
                 //+ " name=\"" + n.name + "\\n"
                 + pathNum + "\"]";
      //String tmp = " [shape=box] " + parentId;
      g += tmp;
    }
    for (int i = 0; i < n.children.size(); ++i) {
      if (!n.children.get(i).isPartOfProofTree) {
        continue;
      }
      String id = Integer.toString(n.children.get(i).id);
      Edge e = new Edge(n.id,n.children.get(i).id);
      
      if (isExist(e)) {
        continue;
      } else {
        edgeSet.add(e);
      }
      //String tmp = "\"" + parentId + '"' + " -- " + "\"" + id + "\"" + " ";
      String tmp = "\"" + parentId + '"' + " -> " + "\"" + id + "\"" + " ";
      g += tmp;
      proofTreeToDot(n.children.get(i));
    }
  }
  
  public void printDot() {
    g += "}";
    System.out.println(g);
  }
  
  public boolean isExist(Edge e) {
    for (Edge ee : edgeSet) {
      if (e.equals(ee)) {
        return true;
      } 
    }
    return false;
  }
  
  public void printAll() {
    for (Map.Entry me : this.allCompounds.entrySet()) {
      System.out.println(me.getKey() + "," + me.getValue());
    }
    for (Map.Entry me : this.allReactions.entrySet()) {
      System.out.println(me.getKey() + "," + me.getValue());
    }
  }
}
