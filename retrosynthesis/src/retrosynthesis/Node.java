package retrosynthesis;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import chemaxon.struc.*;

public class Node {
  public Molecule mol;
  public List<Node> parents = new ArrayList<Node>();
  public String name;
  public String smiles;
  public String reaction;
  private int phi = 1; //ORノードでは証明数
  private int phiWithCost = 1;
  private int delta = 1; //ORノードでは反証数
  private int deltaWithCost = 1;
  //public long phiThreshold = Long.MAX_VALUE - 1;
  public int phiThreshold = Integer.MAX_VALUE / 2 - 1;
  //public long deltaThreshold = Long.MAX_VALUE - 1;
  public int deltaThreshold = Integer.MAX_VALUE / 2 - 1;
  //public int prevPth = Integer.MAX_VALUE / 2 - 1;
  //public int prevDth = Integer.MAX_VALUE / 2 - 1;
  public int prevPth = -1;
  public int prevDth = -1;
  public int depth;
  public List<Node> children = new ArrayList<Node>();
  public ArrayList<Node> pathFromRoot = new ArrayList<Node>();
  public boolean isProven = false;
  public boolean isDisproven = false;
  public boolean isExpanded = false;
  public boolean isRoot = false;
  public boolean isChecked = false;
  public boolean isAvailable = false;
  public boolean ignore = false;
  public boolean isPartOfProofTree = false;
  public String nodeType;
  public int cost; 
  public int id;
  public int degree = 0;
  public int lastSelectNodeIndex = -1;
  public long routeNum = 0;
  private int proofNumberWithCost;
  
  public Node(Molecule m,String nodeType,String tmp,int depth,int id){
    this.mol = m;
    this.phi = 1;
    this.delta = 1;
    this.id = id;
    this.phiThreshold = Integer.MAX_VALUE / 2 - 1;
    //this.phiThreshold = 2000000 - 1;
    this.deltaThreshold = Integer.MAX_VALUE / 2 - 1;
    //this.deltaThreshold = 2000000 - 1;    
    this.depth = depth;
    this.nodeType = nodeType;
    if (nodeType.equals("molecule")) {
      this.smiles = tmp;
      this.cost = 0;
      this.phiWithCost = 1;
      this.deltaWithCost = 1;
      //this.routeNum = 0;
    } else {
      this.reaction = tmp;
      this.cost = 1;
      this.phiWithCost = 2;
      this.deltaWithCost = 2;
      //this.routeNum = 1;
    }
  }
  
  public Node(Molecule m,String nodeType,String tmp,int depth,Node p) {
    //this.parent = p;
    this.parents.add(p);
    this.mol = m;
    this.phi = 1;
    this.delta = 1;
    this.phiThreshold = Integer.MAX_VALUE / 2 - 1;
    this.deltaThreshold = Integer.MAX_VALUE / 2 - 1;
    //this.phiThreshold = 2000000 - 1;
    //this.deltaThreshold = 2000000 - 1;
    this.depth = depth;
    this.phiWithCost = 1;
    this.deltaWithCost = 1;
    this.nodeType = nodeType;
    if (nodeType.equals("molecule")) {
      this.smiles = tmp;
      //this.cost = 0;
      this.phiWithCost = 1;
      this.deltaWithCost = 1;
      //this.routeNum = 0;
    } else {
      //this.cost = 1;
      this.reaction = tmp;
      this.phiWithCost = 2;
      this.deltaWithCost = 2;
      //this.routeNum = 1;
    }
  }
  
  public Node(int id, String nodeType, int depth) {
    this.phi = 1;
    this.delta = 1;
    this.phiThreshold = Integer.MAX_VALUE / 2 - 1;
    this.deltaThreshold = Integer.MAX_VALUE / 2 - 1;
    this.depth = depth;
    this.phiWithCost = 1;
    this.deltaWithCost = 1;
    this.nodeType = nodeType;
    if (nodeType.equals("molecule")) {
      //this.cost = 0;
      this.phiWithCost = 1;
      this.deltaWithCost = 1;
    } else {
      //this.cost = 1;
      this.phiWithCost = 2;
      this.deltaWithCost = 2;
    }
  }
  
  //@Override
  public boolean equals(Node n1) {
    if (smiles.equals(n1.smiles) && depth == n1.depth) {
      return true;
    }
    return false;
  }

  public void setPhi(int p) {
    this.phi = p;
  }

  public void setDelta(int d) {
    this.delta = d;
  }

  public int getPhi() {
    return this.phi;
  }

  public int getDelta() {
    return this.delta;
  }
 
  public int getPhiWithCost() {
    return this.phiWithCost;
    /*
    if (this.nodeType.equals("molecule")) {
      return this.phi + this.cost;
    } else {
      return this.phi;
    }
    */
  }
  
  public int getDeltaWithCost() {
    return this.deltaWithCost;
    /*
    if (this.nodeType.equals("molecule")) {
      return this.delta;
    } else {
      return this.delta + this.cost;
    }
    //return this.delta;
    */
  }
  
  public void setId(int id) {
    this.id = id;
  }
  
  public int getProofNumber() {
    if (this.nodeType.equals("molecule")) {
      return getPhi();
    } else {
      return getDelta();
    }
  }
  /*
  public int getProofNumberWithCost() {
    return this.proofNumberWithCost;
  }
  */
  
  public int getDisProofNumber() {
    if (this.nodeType.equals("reaction")) {
      return getPhi();
    } else {
      return getDelta();
    }
  }
  
  public void updatePhiWithCost() {
    // cost付きphiを更新する
    // deltaMinのようなことをする(deltaはcost付きdelta)
    int tmp = 1000000;
    // 子ノードは存在するが、その全てがignore, 証明, 反証のいずれかだった場合、証明数・反証数の更新規則を変更
    if (this.children.size() > 0 && this.getReachableChildrenSize() == 0) {
      System.out.printf("final update for %d\n", this.id);
      //System.out.println(this.children.size());
      if (this.nodeType.equals("reaction") && this.getDisProofNumber() == 0) {
        this.ignore = true;
      }
      finalUpdate();
      return;
    } else if (this.children.size() == 0) {
      return;
    }
    for (int i = 0; i < this.children.size(); ++i) {
      Node n = this.children.get(i);
      // 証明・反証・ignoreのいずれかのノードは選択候補に入らないので無視
      if (n.getProofNumber() == 0 || n.getDisProofNumber() == 0 || n.ignore) {
        continue;
      }
      if (n.getDeltaWithCost() < tmp) {
        tmp = n.getDeltaWithCost();
        //this.phi = n.delta;
        this.cost = n.cost;
      }
    }
    if (this.nodeType.equals("reaction")) {
      this.cost += 1;
      this.phiWithCost = tmp + 1;
    } else {
      this.phiWithCost = tmp;
    }
  }
  
  public void updateDeltaWithCost() {
    int tmp = 0;
    if (this.children.size() > 0 && this.getReachableChildrenSize() == 0) {
      System.out.printf("final update for %d\n", this.id);
      //System.out.println(this.children.size());
      //this.ignore = true;
      finalUpdate();
      return;
    } else if (this.children.size() == 0) {
      return;
    }
    for (int i = 0; i < this.children.size(); ++i) {
      Node n = this.children.get(i);
      if (n.getProofNumber() == 0 || n.getDisProofNumber() == 0 || n.ignore) {
        continue;
      }
      tmp += n.getPhiWithCost();
    }
    if (this.nodeType.equals("reaction")) {
      this.cost += 1;
      System.out.printf("deltawithcost of id:%d is updated from %d to %d\n", this.id, this.deltaWithCost, tmp + 1);
      this.deltaWithCost = tmp + 1;
    } else {
      System.out.printf("deltawithcost of id:%d is updated from %d to %d\n", this.id, this.deltaWithCost, tmp);
      this.deltaWithCost = tmp;
    }
  }
  
  public void finalUpdate() {
    for (int i = 0; i < this.children.size(); ++i) {
      Node n = this.children.get(i);
      //if ((this.nodeType.equals("molecule") && (n.ignore || n.getProofNumber() == 0)) // 子ANDノードが全部反証されてた時も証明扱いになる
      if ((this.nodeType.equals("molecule") && n.getProofNumber() == 0)
          || (this.nodeType.equals("reaction") && n.getDisProofNumber() == 0)) {
        // ORノードに関して、子ノードが一つでも証明されていたら証明数を0に
        // ANDノードに関して、子ノードが一つでも反証されていたら反証数を0に
        this.phi = 0;
        this.delta = Integer.MAX_VALUE / 2;
        this.phiWithCost = 0;
        this.deltaWithCost = Integer.MAX_VALUE / 2;
        return;
      }
    }
    this.phi = Integer.MAX_VALUE / 2;
    this.delta = 0;
    this.phiWithCost = Integer.MAX_VALUE / 2;
    this.delta = 0;
    return;
  }
  
  public void appendParents(Node n) {
    for (Node p : this.parents) {
      if (n.equals(p)) {
        return;
      }
    }
    this.parents.add(n);
  }

  public int getReachableChildrenSize() {
    int k = this.children.size();
    for (int i = 0; i < this.children.size(); ++i) {
      Node n = this.children.get(i);
      if (this.id == 18) {
        //System.out.printf("%d %d %d %d %d\n", n.id, k, i, n.getProofNumber(), n.getDisProofNumber());
      }
      if (n.ignore || n.getProofNumber() == 0 || n.getDisProofNumber() == 0) {
        k = k - 1;
      }
    }
    if (this.id == 18) {
      //System.out.printf("Reachable children size of %d is %d.\n", this.id, k);
    }
    return k;
  }
  
  public void updateRouteNum() {
    System.out.printf("route path update for %d, its current routeNum is %d\n", this.id, this.routeNum);
    System.out.printf("children size:%d\n", this.children.size());
    long var = 0;
    if (this.children.size() == 0) {
      var = this.routeNum;
    }
    if (this.nodeType.equals("reaction")) {
      var = 1;
    }
    for (int i = 0; i < this.children.size(); ++i) {
      Node n = this.children.get(i);
      //System.out.printf("routeNum of %d is %d\n",this.children.get(i).id, this.children.get(i).routeNum);
      if (this.nodeType.equals("molecule")) {
        if (n.isPartOfProofTree) {
          var += n.routeNum;
        }
      } else {
        var *= n.routeNum;
        System.out.printf("var:%d, %d\n", var, n.routeNum);
      }
    }
    this.routeNum = var;
  }

  public boolean isMakingCycle(Node next) {
    boolean isCycle = this.detectCycleByPath(next);
    return isCycle;
    /*
    if (isCycle) {
      return isCycle;
    } else {
      //System.out.println("dfs");
      isCycle = this.detectCycleByDFS(next);
      if (isCycle) {
        return isCycle;
      }
    }
    return isCycle;
    */
  }
  
  public boolean detectCycleByPath(Node next) {
    //assert this.nodeType.equals("reaction");
    int n = this.pathFromRoot.size();
    for (int i = 0; i < this.pathFromRoot.size(); ++i) {
      //System.out.println(this.pathFromRoot.get(i).id);
    }
    //System.out.printf("next:%d\n", next.id);
    //System.out.println("Are you killing me!?");
    for (int i = 0; i < n; ++i) {
      //System.out.printf("next:%d history:%d\n", next, this.pathFromRoot.get(i));
      for (int j = 0; j < next.children.size(); ++j) {
        Node c = next.children.get(j);
        if (c.id == this.pathFromRoot.get(i).id) {
          for (int k = 0; k < this.children.size(); ++k) {
            if (this.children.get(k).id == next.id) {
              next.setPhi(0);
              next.setDelta(Integer.MAX_VALUE / 2);
              return true;
            }
          }
        }
      }
    }
    return false;
  }
  
  public boolean detectCycleByDFS(Node next) {
    Deque<Node> stack = new ArrayDeque<Node>();
    Set<Integer> nodeSet = new HashSet<Integer>();
    Node startNode = null;
    for (int i = 0; i < this.children.size(); ++i) {
      if (this.children.get(i).id == next.id) {
        startNode = this.children.get(i);
      }
    }
    stack.add(startNode);
    int count = 0;
    while (stack.size() != 0) {
      Node tmp = stack.pop();
      if (!nodeSet.contains(tmp.id)) {
        nodeSet.add(tmp.id);
        for (int i = 0; i < tmp.children.size(); ++i) {
          Node c = tmp.children.get(i);
          if (c.id == this.id) {
            Node p = tmp.parents.get(0);
            for (int j = 0; j < tmp.children.size(); ++j) {
              if (p.children.get(j).id == tmp.id) {
                System.out.printf("removed id:%d\n", tmp.id);
                p.children.remove(j);
                break;
              }
            }
            return true;
          }
          stack.push(c);
        }
      }
    }
    return false;
  }
}
