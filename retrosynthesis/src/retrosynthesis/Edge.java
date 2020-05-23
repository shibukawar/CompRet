package retrosynthesis;

public class Edge {
  public int v1;
  public int v2;
  
  Edge(int s, int t) {
    this.v1 = s;
    this.v2 = t;
  }
  
  
  boolean equals(Edge e) {
    if ((this.v1 == e.v1 && this.v2 == e.v2) || (this.v1 == e.v2 && this.v2 == e.v1)) {
      return true;
    } else {
      return false;
    }
  }
  
}
