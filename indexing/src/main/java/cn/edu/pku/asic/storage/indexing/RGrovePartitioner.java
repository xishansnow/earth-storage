package cn.edu.pku.asic.storage.indexing;

import cn.edu.pku.asic.storage.common.cg.SpatialPartitioner;
import cn.edu.pku.asic.storage.common.cli.BeastOptions;
import cn.edu.pku.asic.storage.common.geolite.EnvelopeNDLite;
import cn.edu.pku.asic.storage.common.geolite.GeometryHelper;
import cn.edu.pku.asic.storage.common.synopses.AbstractHistogram;
import cn.edu.pku.asic.storage.common.synopses.Summary;
import cn.edu.pku.asic.storage.common.utils.IntArray;
import cn.edu.pku.asic.storage.common.utils.OperationParam;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * A partitioner that uses the R-Grove partitioning algorithm. It is based on the R-tree split node function but
 * it ensures that the configured minimum and maximum node sizes are met.
 */
@SpatialPartitioner.Metadata(
    disjointSupported = false,
    extension = "rgrove",
    description = "Recursively applies the R-tree node splitting algorithm to split a sample of points into partitions.")
public class RGrovePartitioner implements SpatialPartitioner {

  @OperationParam(
      description = "The desired ratio between the minimum and maximum partitions sizes ]0,1[",
      defaultValue = "0.95",
      required = false
  )
  public static final String MMRatio = "mmratio";

  @OperationParam(
          description = "The minimum fraction to use when applying the linear-time split algorithm",
          defaultValue = "0.95",
          required = false
  )
  public static final String MinSplitRatio = "RGrove.MinSplitRatio";

  /**The ratio between the minimum and maximum split sizes*/
  protected float mMRatio;

  /**The minimum fraction to use when applying the linear-time split algorithm*/
  protected float fractionMinSplitSize;

  /**The minimum and maximum coordinates for the partitions*/
  protected double[][] minCoord, maxCoord;

  /**MBR of the input space*/
  protected final EnvelopeNDLite inputMBR = new EnvelopeNDLite();

  @Override
  public void setup(BeastOptions conf, boolean disjoint) {
    if (disjoint)
      throw new RuntimeException("Disjoint partitions are not supported by the R-Grove partitioner");
    this.mMRatio = conf.getFloat(MMRatio, 0.95f);
    this.fractionMinSplitSize = conf.getFloat(MinSplitRatio, 0.0f);
  }

  @Override
  public void construct(Summary summary, @Required double[][] sample, AbstractHistogram histogram, int numPartitions) {
    inputMBR.set(summary);
    // Use the R-tree node splitting algorithm recursively on the sample points
    int numSamplePoints = sample[0].length;
    int numDimensions = sample.length;
    int M = (int) Math.ceil((double)numSamplePoints / numPartitions);
    int m = (int) Math.ceil(mMRatio * M);
    EnvelopeNDLite[] partitionMBRs = RTreeGuttman.partitionPoints(sample, m, M, fractionMinSplitSize);

    minCoord = new double[numDimensions][partitionMBRs.length];
    maxCoord = new double[numDimensions][partitionMBRs.length];
    for (int i = 0; i < partitionMBRs.length; i++) {
      for (int d = 0; d < numDimensions; d++) {
        minCoord[d][i] = (partitionMBRs[i]).getMinCoord(d);
        maxCoord[d][i] = (partitionMBRs[i]).getMaxCoord(d);
      }
    }
  }

  /**
   * Tests if a partition overlaps a given rectangle
   * @param partitionID the ID of the partition to compute its overlap
   * @param ienv the envelope that needs to be tested with the partition
   * @return {@code true} iff the given envelope overlaps the partition boundaries
   */
  protected boolean Partition_overlap(int partitionID, EnvelopeNDLite ienv) {
    for (int d = 0; d < getCoordinateDimension(); d++) {
      if (maxCoord[d][partitionID] <= ienv.getMinCoord(d) || ienv.getMaxCoord(d) <= minCoord[d][partitionID])
        return false;
    }
    return true;
  }

  /**
   * Computes the area of a partition.
   * @param partitionID the ID of the partition to compute its volume
   * @return the volume of the partition
   */
  protected double Partition_volume(int partitionID) {
    double vol = 1.0;
    for (int d = 0; d < getCoordinateDimension(); d++)
      vol *= maxCoord[d][partitionID] - minCoord[d][partitionID];
    return vol;
  }

  /**
   * Computes the expansion that will happen on an a partition when it is
   * enlarged to enclose a given rectangle.
   * @param iPartition the ID of the partition to check its expansion
   * @param env the MBR of the object to be added to the partition
   * @return the amount of expansion that will happen if the given envelope is added to the given partition.
   */
  protected double Partition_expansion(int iPartition, EnvelopeNDLite env) {
    double volBefore = 1.0, volAfter = 1.0;
    assert env.getCoordinateDimension() == this.getCoordinateDimension();
    for (int d = 0; d < getCoordinateDimension(); d++) {
      volBefore *= maxCoord[d][iPartition] - minCoord[d][iPartition];
      volAfter *= Math.max(maxCoord[d][iPartition], env.getMaxCoord(d)) -
          Math.min(minCoord[d][iPartition], env.getMinCoord(d));
    }
    return volAfter - volBefore;
  }

  @Override
  public void overlapPartitions(EnvelopeNDLite mbr, IntArray matchedPartitions) {
    matchedPartitions.clear();
    for (int i = 0; i < minCoord[0].length; i++) {
      if (Partition_overlap(i, mbr))
        matchedPartitions.add(i);
    }
  }

  @Override
  public int overlapPartition(EnvelopeNDLite mbr) {
    double minExpansion = Double.POSITIVE_INFINITY;
    int chosenPartition = -1;
    for (int iPartition = 0; iPartition < minCoord[0].length; iPartition++) {
      double expansion = Partition_expansion(iPartition, mbr);
      if (expansion < minExpansion) {
        minExpansion = expansion;
        chosenPartition = iPartition;
      } else if (expansion == minExpansion) {
        // Resolve ties by choosing the entry with the rectangle of smallest area
        if (Partition_volume(iPartition) < Partition_volume(chosenPartition))
          chosenPartition = iPartition;
      }
    }
    return chosenPartition;
  }

  @Override
  public void getPartitionMBR(int partitionID, EnvelopeNDLite mbr) {
    for (int d = 0; d < getCoordinateDimension(); d++) {
      mbr.setMinCoord(d, this.minCoord[d][partitionID]);
      mbr.setMaxCoord(d, this.maxCoord[d][partitionID]);
    }
  }

  @Override
  public int getPartitionCount() {
    return minCoord == null? 0 : minCoord[0].length;
  }

  @Override
  public boolean isDisjoint() {
    return false;
  }

  @Override
  public int getCoordinateDimension() {
    return minCoord == null? 0 : minCoord.length;
  }

  @Override
  public void writeExternal(ObjectOutput out) throws IOException {
    GeometryHelper.writeIEnvelope(inputMBR, out);
    out.writeInt(getCoordinateDimension());
    out.writeInt(getPartitionCount());
    for (int d = 0; d < getCoordinateDimension(); d++) {
      for (int i = 0; i < getPartitionCount(); i++) {
        out.writeDouble(minCoord[d][i]);
        out.writeDouble(maxCoord[d][i]);
      }
    }
  }

  @Override
  public void readExternal(ObjectInput in) throws IOException {
    GeometryHelper.readIEnvelope(this.inputMBR, in);
    int numDimensions = in.readInt();
    int numPartitions = in.readInt();
    if (minCoord == null || getPartitionCount() != numPartitions || getCoordinateDimension() != numDimensions) {
      minCoord = new double[numDimensions][numPartitions];
      maxCoord = new double[numDimensions][numPartitions];
    }
    for (int d = 0; d < getCoordinateDimension(); d++) {
      for (int i = 0; i < getPartitionCount(); i++) {
        minCoord[d][i] = in.readDouble();
        maxCoord[d][i] = in.readDouble();
      }
    }
  }

  @Override
  public EnvelopeNDLite getEnvelope() {
    return inputMBR;
  }
}
