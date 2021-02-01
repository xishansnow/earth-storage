/*
 * Copyright 2020 University of California, Riverside
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
 */package cn.edu.pku.asic.storage.indexing;

import cn.edu.pku.asic.storage.common.cg.HCurveSortable;
import cn.edu.pku.asic.storage.common.cg.SpatialPartitioner;
import cn.edu.pku.asic.storage.common.cli.AppOptions;
import cn.edu.pku.asic.storage.common.geolite.EnvelopeND;
import cn.edu.pku.asic.storage.common.geolite.EnvelopeNDLite;
import cn.edu.pku.asic.storage.common.geolite.GeometryHelper;
import cn.edu.pku.asic.storage.common.synopses.AbstractHistogram;
import cn.edu.pku.asic.storage.common.synopses.Summary;
import cn.edu.pku.asic.storage.common.utils.IntArray;
import org.apache.hadoop.util.QuickSort;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * A partitioner that partitions the space based on a the Morton space-filling curve, a.k.a., Z-Curve.
 */
@SpatialPartitioner.Metadata(
    disjointSupported = false,
    extension = "hcurve",
    description = "Sorts the sample points based on a Hilbert-curve and partitions the Hilbert-space into ranges with roughly equal number of points")
public class HCurvePartitioner implements SpatialPartitioner {
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
  protected HCurveSortable hCurveSorter;

  @Override
  public void setup(AppOptions conf, boolean disjoint) {
    if (disjoint)
      throw new RuntimeException("H-Cruve partitioner does not support disjoint partitions");
  }

  @Override
  public void construct(Summary summary, @Required double[][] sample, AbstractHistogram histogram, int numPartitions) {
    this.inputMBR.set(summary);
    final int numDimensions = sample.length;
    int sampleCount = sample[0].length;
    assert numPartitions > 0;
    if (numPartitions > 1) {
      // Sort all the points based on a Z-curve
      hCurveSorter = new HCurveSortable(summary);
      hCurveSorter.setPoints(sample);
      new QuickSort().sort(hCurveSorter, 0, sampleCount);

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
      hCurveSorter.setPoints(splitPoints);
    }
  }

  @Override
  public void overlapPartitions(EnvelopeNDLite mbr, IntArray matchedPartitions) {
    throw new RuntimeException("H-Curve partitioner does not support disjoint partitions");
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
      int c = hCurveSorter.compare(searchPoint, middle);
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
      hCurveSorter = null;
      return;
    }
    splitPoints = new double[inputMBR.getCoordinateDimension()][numSplits + 1];
    if (numSplits > 0) {
      hCurveSorter = new HCurveSortable(inputMBR);
      hCurveSorter.setPoints(splitPoints);
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
