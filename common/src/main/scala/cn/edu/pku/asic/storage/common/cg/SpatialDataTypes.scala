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
 */
package cn.edu.pku.asic.storage.common.cg

import org.apache.spark.api.java.{JavaPairRDD, JavaRDD}
import org.apache.spark.rdd.RDD
import cn.edu.pku.asic.storage.common.geolite.IFeature

trait SpatialDataTypesMixin {

  /**A type alias for spatial RDDs*/
  type SpatialRDD = RDD[IFeature]

  /**A type alias for Java spatial RDDs*/
  type JavaSpatialRDD = JavaRDD[IFeature]

  /**
   * A type alias for a partitioned spatial RDD. The key is the partition number and the SpatialPartitioner
   * defines the boundaries of each partition
   */
  type PartitionedSpatialRDD = RDD[(Int, IFeature)]

  /**
   * A type alias for a partitioned spatial RDD in Java
   */
  type JavaPartitionedSpatialRDD = JavaPairRDD[Integer, IFeature]
}

/**
 * Create an object out of the mixin to enable import of the data types
 */
object SpatialDataTypes extends SpatialDataTypesMixin
