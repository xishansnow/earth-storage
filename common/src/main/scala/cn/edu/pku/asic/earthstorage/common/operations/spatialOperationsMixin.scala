package cn.edu.pku.asic.earthstorage.common.operations

import cn.edu.pku.asic.earthstorage.common.cg.{CGOperationsMixin, SparkSpatialPartitioner}
import cn.edu.pku.asic.earthstorage.common.cg.SpatialDataTypes.{PartitionedSpatialRDD, SpatialRDD}
import cn.edu.pku.asic.earthstorage.common.cg.SpatialJoinAlgorithms.{ESJDistributedAlgorithm, ESJPredicate}
import cn.edu.pku.asic.earthstorage.common.geolite.IFeature
import org.apache.spark.rdd.RDD
import org.apache.spark.util.LongAccumulator

/**
 * Adds spatial operations to RDDs
 */
trait SpatialOperationsMixin extends CGOperationsMixin {

  /**
   * Additional functions for SpatialRDD
   *
   * @param rdd the underlying RDD
   */
  implicit class RDDSpatialFunctions2(rdd: SpatialRDD) {

    /**
     * Performs a spatial join between this RDD and another RDD
     *
     * @param rdd2          another RDD to be join with
     * @param joinPredicate the spatial join predicate
     * @param method        the spatial join algorithm
     * @param mbrCount      an optional counter to keep track of the number of MBR tests
     * @return a new RDD that contains pairs of matching features
     */
    /*def spatialJoin(rdd2: SpatialRDD, joinPredicate: ESJPredicate = ESJPredicate.Intersects,
                    method: ESJDistributedAlgorithm = null,
                    mbrCount: LongAccumulator = null): RDD[(IFeature, IFeature)] = {

            SpatiatlJoin.spatialJoin(rdd, rdd2, joinPredicate, method, mbrCount)
    }*/
  }

  /**
   * Shortcut functions for SpatialRDDs that are partitioned by any spatial partitioner
   *
   * @param partitionedRDD
   */
  implicit class RDDPartitionedSpatialFunctions2(partitionedRDD: PartitionedSpatialRDD) {
    require(partitionedRDD.partitioner.isDefined && partitionedRDD.partitioner.get.isInstanceOf[SparkSpatialPartitioner],
      "This function requires the RDD to be partitioned by a spatial partitioner")

    /**
     * Performs a spatial join between two spatially partitioned datasets. This method always uses the
     * distributed join algorithm since both input datasets are spatially partitioned.
     *
     * @param partitionedRDD2
     * @param joinPredicate
     * @param mbrCount
     * @return
     */
    /*def spatialJoin(partitionedRDD2: PartitionedSpatialRDD, joinPredicate: ESJPredicate = ESJPredicate.Intersects,
                    mbrCount: LongAccumulator = null): RDD[(IFeature, IFeature)] = {
      require(partitionedRDD2.partitioner.isDefined && partitionedRDD2.partitioner.get.isInstanceOf[SparkSpatialPartitioner],
        "This function requires the RDD to be partitioned by a spatial partitioner")
     // SpatialJoin.spatialJoinDJPartitionedRDDs(partitionedRDD, partitionedRDD2, joinPredicate, mbrCount)
    }*/
  }

}
