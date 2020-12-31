package cn.edu.pku.coe.asic.indexing;

import edu.ucr.cs.bdlab.indexing.SpatialPartitioner;
import edu.ucr.cs.bdlab.stsynopses.AbstractHistogram;
import edu.ucr.cs.bdlab.stsynopses.Summary;
import org.apache.hadoop.conf.Configuration;
import edu.ucr.cs.bdlab.geolite.EnvelopeND;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;


public class DggsPartitioner implements  SpatialPartitioner {

    @Override
    public void setup(Configuration conf, boolean disjoint, PartitionCriterion pc, long value) {

    }

    @Override
    public void construct(Summary summary, double[][] doubles, AbstractHistogram abstractHistogram) {

    }

    @Override
    public void writeExternal(ObjectOutput out) throws IOException {

    }

    @Override
    public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {

    }

    @Override
    public int getPartitionCount() {
        return 0;
    }

    @Override
    public int getCoordinateDimension() {
        return 0;
    }

}
