/*
 * Copyright 2021 University of California, Riverside
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
package cn.edu.pku.asic.storage.common.synopses

import cn.edu.pku.asic.storage.common.cg.SpatialDataTypes.SpatialRDD
import cn.edu.pku.asic.storage.common.cg.WrappedSpatialPartition
import cn.edu.pku.asic.storage.common.geolite.IFeature
import org.apache.spark.{Partition, Partitioner, TaskContext}
/**
 * A [[SpatialRDD]] that has a summary associated with each partition.
 */
class SummarizedSpatialRDD(@transient _parentRDD: SpatialRDD, @transient summaries: Array[Summary])
  extends SpatialRDD(_parentRDD) {

  require(_parentRDD.getNumPartitions == summaries.length, "Should have one summary for each partition")

  override protected def getPartitions: Array[Partition] = {
    val ps = new Array[Partition](firstParent.getNumPartitions)
    for (i <- ps.indices) {
      ps(i) = new WrappedSpatialPartition(firstParent.partitions(i), summaries(i))
    }
    ps
  }

  override def compute(split: Partition, context: TaskContext): Iterator[IFeature] = {
    firstParent.compute(split.asInstanceOf[WrappedSpatialPartition].partition, context)
  }

  override val partitioner: Option[Partitioner] = firstParent.partitioner
}
