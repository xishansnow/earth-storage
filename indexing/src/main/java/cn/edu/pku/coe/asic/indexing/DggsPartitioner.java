package cn.edu.pku.coe.asic.indexing;


import edu.ucr.cs.bdlab.beast.cg.SpatialPartitioner;
import edu.ucr.cs.bdlab.beast.common.BeastOptions;
import edu.ucr.cs.bdlab.beast.geolite.EnvelopeND;
import edu.ucr.cs.bdlab.beast.geolite.EnvelopeNDLite;
import edu.ucr.cs.bdlab.beast.synopses.AbstractHistogram;
import edu.ucr.cs.bdlab.beast.synopses.Summary;
import edu.ucr.cs.bdlab.beast.util.IntArray;


import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;
import java.io.PrintStream;
import java.util.Iterator;


public abstract class DggsPartitioner implements  SpatialPartitioner {

    @Override
   public void setup(BeastOptions conf, boolean disjoint) {


    }

    @Override
    public void construct(Summary summary, double[][] doubles, AbstractHistogram abstractHistogram, int i) {


    }
}
