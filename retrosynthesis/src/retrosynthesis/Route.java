package retrosynthesis;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class Route {
  public List<String> route = new ArrayList<String>();
  public Set<String> startingMaterials = new HashSet<String>();

  public Route() {}

  public void addRoute(String str) {
    route.add(str);
  }
  
  public void addStartingMaterials(String sml) {
    startingMaterials.add(sml);
  }

  public void printRoute() {
    for (int i = 0; i < route.size(); ++i) {
      System.out.println(route.get(i));
    }
  }
  
  public void printStartingMaterials() {
    for (String s : startingMaterials ) {
      System.out.println(s);
    }
  }
}
