package retrosynthesis;

import java.io.IOException;

import chemaxon.formats.*;
import chemaxon.struc.*;

public class Main {
  public static void main(String[] args) throws Exception { 
    System.out.println(args[0]);
    int d = Integer.parseInt(args[0]); 
    long t = Long.parseLong(args[1]);
    long rn = Long.parseLong(args[2]);
    String r = args[3];
    String b = args[4];
    String sml = args[5];
    Molecule a = MolImporter.importMol(sml,"smiles");
    Tree dfpn = new Tree(a,true, d, t, rn, r, b);
    dfpn.depthFirstProofNumberSearch();
  }
}
