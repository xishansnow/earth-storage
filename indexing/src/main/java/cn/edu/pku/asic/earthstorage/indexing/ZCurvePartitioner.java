/*
 * Copyright 2018 University of California, Riverside
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package cn.edu.pku.asic.earthstorage.indexing;

import cn.edu.pku.asic.earthstorage.common.cg.SpatialPartitioner;
import cn.edu.pku.asic.earthstorage.common.cli.AppOptions;
import cn.edu.pku.asic.earthstorage.common.geolite.EnvelopeND;
import cn.edu.pku.asic.earthstorage.common.geolite.EnvelopeNDLite;
import cn.edu.pku.asic.earthstorage.common.geolite.GeometryHelper;
import cn.edu.pku.asic.earthstorage.common.synopses.AbstractHistogram;
import cn.edu.pku.asic.earthstorage.common.synopses.Summary;
import cn.edu.pku.asic.earthstorage.common.utils.IntArray;
import org.apache.hadoop.util.IndexedSortable;
import org.apache.hadoop.util.QuickSort;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * A partitioner that partitions the space based on a the Morton space-filling curve, a.k.a., Z-Curve.
 */
@SpatialPartitioner.Metadata(
    disjointSupported = false,
    extension = "zcurve",
    description = "Sorts the sample points based on a Z-curve and partitions the Z-space into ranges with roughly equal number of points")
public class ZCurvePartitioner implements SpatialPartitioner {
  /**
   * The coordinates of the points that define the splits.
   * The first index is for the dimension an the second is for the coordinate.
   * The last point in this array is reserved for the use in the {@link #overlapPartition(EnvelopeND)} method to avoid
   * creating an array each time it is called.
   */
  protected double[][] splitPoints;

  /**
   * The MBR of the input domain.
   */
  final protected EnvelopeNDLite inputMBR = new EnvelopeNDLite();

  /**
   * A sorter that is used to initialize the partitions and to search for partitions as well.
   */
  protected ZCurveSortable zCurveSorter;

  /**
   * A class for comparing and sorting points based on the Z-Curve
   */
  protected static class ZCurveSortable implements IndexedSortable {

    /**A temporary array for holding the position of the points to compare*/
    int[] centeri;
    int[] centerj;
    /**The MBR of the entire input space*/
    EnvelopeNDLite mbr;
    /**The array of point coordinate to compare. The first index is for the dimensions and the second for the coordinate*/
    double[][] points;

    public ZCurveSortable(EnvelopeNDLite mbr) {
      this.mbr = mbr;
      centeri = new int[mbr.getCoordinateDimension()];
      centerj = new int[mbr.getCoordinateDimension()];
    }

    @Override
    public int compare(int i, int j) {
      // The mask that retains the bit to compare to. Initialized based on the most significant bit.
      int comparisonMask = 0;
      // Calculate the position in the grid for both points
      for (int d = 0; d < mbr.getCoordinateDimension(); d++) {
        centeri[d] = (int) (((points[d][i] - mbr.getMinCoord(d)) * Integer.MAX_VALUE) / (mbr.getSideLength(d)));
        centerj[d] = (int) (((points[d][j] - mbr.getMinCoord(d)) * Integer.MAX_VALUE) / (mbr.getSideLength(d)));
        comparisonMask = Math.max(comparisonMask, Integer.highestOneBit(centeri[d]));
        comparisonMask = Math.max(comparisonMask, Integer.highestOneBit(centerj[d]));
      }
      // Compare the two positions bit-by-bit starting with the most significant bit
      while (comparisonMask > 0) {
        // Compare the most significant bit in all the axes in order
        for (int d = 0; d < mbr.getCoordinateDimension(); d++) {
          int diff = (centeri[d] & comparisonMask) - (centerj[d] & comparisonMask);
          if (diff != 0)
            return diff;
        }
        // Try the less significant bit
        comparisonMask >>>= 1;
      }
      // Equal points according to the Z-curve, i.e., in the same cell
      return 0;
    }

    @Override
    public void swap(int i, int j) {
      double t;
      for (int d = 0; d < points.length; d++) {
        t = points[d][i];
        points[d][i] = points[d][j];
        points[d][j] = t;
      }
    }
  }

  @Override
  public void setup(AppOptions conf, boolean disjoint) {
    if (disjoint)
      throw new RuntimeException("Z-Cruve partitioner does not support disjoint partitions");
  }

  @Override
  public void construct(Summary summary, @Required double[][] sample, AbstractHistogram histogram, int numPartitions) {
    this.inputMBR.set(summary);
    final int numDimensions = sample.length;
    int sampleCount = sample[0].length;
    assert numPartitions > 0;
    if (numPartitions > 1) {
      // Sort all the points based on a Z-curve
      zCurveSorter = new ZCurveSortable(summary);
      zCurveSorter.points = sample;
      new QuickSort().sort(zCurveSorter, 0, sampleCount);
      // Create an array to hold the split points
      int numSplits = numPartitions - 1;
      // We reserve an additional coordinate to use with the overlapPartition method for efficiency
      splitPoints = new double[numDimensions][numSplits + 1];
      // Store the split points
      for (int i = 0; i < numSplits; i++) {
        int splitPoint = sampleCount * (i + 1) / numPartitions;
        for (int d = 0; d < numDimensions; d++)
          splitPoints[d][i] = sample[d][splitPoint];
      }
      // prepare the sorter to search in the partitions
      zCurveSorter.points = splitPoints;
    }
  }

  @Override
  public void overlapPartitions(EnvelopeNDLite mbr, IntArray matchedPartitions) {
    throw new RuntimeException("Z-Curve partitioner does not support disjoint partitions");
  }

  @Override
  synchronized public int overlapPartition(EnvelopeNDLite mbr) {
    // TODO can we change this method to be thread-safe?
    // To do so, we need to remove the additional entry in splitPoints and the temporary arrays in the comparator

    if (splitPoints == null || splitPoints.length == 0)
      return 0; // No splits. Just one big partition that contains all the data.
    assert mbr.getCoordinateDimension() == this.getCoordinateDimension();
    // Set the last point after the end to the center of the given MBR
    int searchPoint = splitPoints[0].length - 1;
    for (int d = 0; d < this.getCoordinateDimension(); d++) {
      splitPoints[d][searchPoint] = (mbr.getMinCoord(d) + mbr.getMaxCoord(d)) / 2.0;
    }
    // Run a binary search to find the index of the split directly after the given point
    int start = 0, end = splitPoints[0].length - 1;
    while (start < end) {
      int middle = (start + end) / 2;
      int c = zCurveSorter.compare(searchPoint, middle);
      if (c <= 0)
        end = middle;
      else
        start = middle + 1;
    }
    return start;
  }

  @Override
  public void getPartitionMBR(int partitionID, EnvelopeNDLite mbr) {
    // This method is hard to implement. For simplicity, we just return the largest possible MBR
    for (int d = 0; d < getCoordinateDimension(); d++) {
      mbr.setMinCoord(d, Double.NEGATIVE_INFINITY);
      mbr.setMaxCoord(d, Double.POSITIVE_INFINITY);
    }
  }

  @Override
  public int getPartitionCount() {
    return splitPoints == null? 1 : (splitPoints[0].length + 1);
  }

  @Override
  public boolean isDisjoint() {
    // This partitions does not support disjoint partitions
    return false;
  }

  @Override
  public int getCoordinateDimension() {
    return inputMBR.getCoordinateDimension();
  }

  @Override
  public void writeExternal(ObjectOutput out) throws IOException {
    GeometryHelper.writeIEnvelope(inputMBR, out);
    int numSplits = splitPoints == null? 0 : splitPoints[0].length - 1; // Last entry is only used temporarily
    out.writeInt(numSplits);
    for (int d = 0; d < getCoordinateDimension(); d++) {
      for (int i = 0; i < numSplits; i++) {
        out.writeDouble(splitPoints[d][i]);
      }
    }
  }

  @Override
  public void readExternal(ObjectInput in) throws IOException {
    GeometryHelper.readIEnvelope(inputMBR, in);
    int numSplits = in.readInt();
    if (numSplits == 0) {
      splitPoints = null;
      zCurveSorter = null;
      return;
    }
    splitPoints = new double[inputMBR.getCoordinateDimension()][numSplits + 1];
    if (numSplits > 0) {
      zCurveSorter = new ZCurveSortable(inputMBR);
      zCurveSorter.points = splitPoints;
      for (int d = 0; d < getCoordinateDimension(); d++) {
        for (int i = 0; i < numSplits; i++) {
          splitPoints[d][i] = in.readDouble();
        }
      }
    }
  }

  @Override
  public EnvelopeNDLite getEnvelope() {
    return inputMBR;
  }
}
