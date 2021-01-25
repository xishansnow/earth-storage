package cn.edu.pku.asic.storage.indexing;

import edu.ucr.cs.bdlab.beast.cg.SpatialPartitioner;
import edu.ucr.cs.bdlab.beast.geolite.EnvelopeNDLite;

/**
 * A spatial partitioner that uses the RR*-Grove partitioning algorithm. Simply, it applies the RR*-tree node splitting
 * method iteratively until each partition contains betweem [m, M] records.
 */
@SpatialPartitioner.Metadata(
    disjointSupported = true,
    extension = "rrsgrove",
    description = "A partitioner that uses the RR*-tree node splitting algorithm on a sample of points to partition the space"
)
public class RRSGrovePartitioner extends edu.ucr.cs.bdlab.beast.indexing.RSGrovePartitioner {

  protected EnvelopeNDLite[] partitionPoints(double[][] coords, int max, int min) {
    return edu.ucr.cs.bdlab.beast.indexing.RRStarTree.partitionPoints(coords, min, max, true, fractionMinSplitSize, aux);
  }
}
