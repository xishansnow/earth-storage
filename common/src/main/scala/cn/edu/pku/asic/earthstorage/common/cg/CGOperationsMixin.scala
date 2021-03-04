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
package cn.edu.pku.asic.earthstorage.common.cg


import cn.edu.pku.asic.earthstorage.common.cg.Reprojector.TransformationInfo
import cn.edu.pku.asic.earthstorage.common.cg.SpatialDataTypes.SpatialRDD
import cn.edu.pku.asic.earthstorage.common.geolite.{EnvelopeNDLite, Feature, GeometryHelper, IFeature}
import cn.edu.pku.asic.earthstorage.common.synopses
import cn.edu.pku.asic.earthstorage.common.synopses._
import org.apache.spark.beast.CRSServer
import org.opengis.referencing.crs.CoordinateReferenceSystem

/**
 * A mixin that defines implicit functions to run spatial operations on Spatial RDDs
 */
trait CGOperationsMixin {

  /**
   * Additional functions for SpatialRDD
   * @param rdd the underlying RDD
   */
  implicit class RDDCGOperations(rdd: SpatialRDD) {
    /**
     * Tells whether a SpatialRDD is partitioned using any spatial partitioner or not
     * @return {@code true} if the RDD is partitioned using any spatial partitioner
     */
    def isSpatiallyPartitioned: Boolean =
      rdd.partitioner.isDefined && rdd.partitioner.get.isInstanceOf[SparkSpatialPartitioner]

    /**
     * Compute the geometric summary of a set of features which includes size (in bytes), number of records,
     * number of points, number of non empty geometries, average side length (width and height), and the geometry type.
     *
     * @return the computed summaries
     */
    def summary: Summary = Summary.computeForFeatures(rdd)

    /**
     * Computes a uniform histogram with the given size that counts number of features in each cell
     *
     * @param histogramSize the size of the histogram as the number of partitions along each dimension
     * @param prefixSum     compute the prefix sum on the result to speed up range tests
     * @return the created histogram
     */
    def uniformHistogramCount(histogramSize: Array[Int], prefixSum: Boolean = false): AbstractHistogram =
      rdd.uniformHistogramSize(histogramSize, prefixSum, sizeFunction = _ => 1)

    /**
     * Computes a uniform histogram with the given size that calculates the size of the data in each cell
     *
     * @param histogramSize the size of the histogram as the number of partitions along each dimension
     * @param prefixSum     compute the prefix sum on the result to speed up range tests
     * @param sizeFunction  an optional function that computes the size of a feature.
     * @return the created histogram
     */
    def uniformHistogramSize(histogramSize: Array[Int], prefixSum: Boolean = false,
                             sizeFunction: IFeature => Int = _.getStorageSize): AbstractHistogram = {
      val numBins = histogramSize.reduce(_*_)
      val computationMethod: HistogramOP.ComputationMethod = if (numBins < 100000) HistogramOP.TwoPass else HistogramOP.Sparse
      val histogram = HistogramOP.computeHistogram(rdd, sizeFunction, computationMethod,
        HistogramOP.PointHistogram, histogramSize:_*)
      if (prefixSum)
        new Prefix2DHistogram(histogram.asInstanceOf[UniformHistogram])
      else
        histogram
    }

    /**
     * Computes an Euler histogram that works better for geometries with extents, i.e., envelopes,
     * which calculates the number of records in each cell
     *
     * @param histogramSize the size of the histogram as the number of partitions along each dimension
     * @param prefixSum     compute the prefix sum on the result to speed up range tests
     * @return the created histogram
     */
    def eulerHistogramCount(histogramSize: Array[Int], prefixSum: Boolean = false): AbstractHistogram =
      rdd.eulerHistogramSize(histogramSize, prefixSum, _ => 1)

    /**
     * Computes an Euler histogram that works better for geometries with extents, i.e., envelopes,
     * which calculates the total size of features in each cell
     *
     * @param histogramSize the size of the histogram as the number of partitions along each dimension
     * @param prefixSum     compute the prefix sum on the result to speed up range tests
     * @return the created histogram
     */
    def eulerHistogramSize(histogramSize: Array[Int], prefixSum: Boolean = false,
                           sizeFunction: IFeature => Int = _.getStorageSize): AbstractHistogram = {
      val histogram = HistogramOP.computeHistogram(rdd, sizeFunction, HistogramOP.TwoPass,
        HistogramOP.EulerHistogram, histogramSize:_*)
      if (prefixSum)
        new synopses.PrefixEulerHistogram2D(histogram.asInstanceOf[EulerHistogram2D])
      else
        histogram
    }


    /**
     * Reproject the geometries in the given RDD to from the source CRS to the target CRS
     * @param sourceCRS the course coordinate reference system (CRS). The function will assume that all geometries
     *                  in this RDD are in the sourceCRS
     * @param targetCRS the target coordinate reference system (CRS). The returned RDD will have all geometries
     *                  in the targetCRS
     * @return the RDD that contains the converted geometries
     */
    def reproject(sourceCRS: CoordinateReferenceSystem, targetCRS: CoordinateReferenceSystem): SpatialRDD = {
      val transform: TransformationInfo = Reprojector.findTransformationInfo(sourceCRS, targetCRS, rdd.sparkContext.getConf)
      rdd.map(f => {
        val transformedGeometry = Reprojector.reprojectGeometry(f.getGeometry, transform)
        Feature.create(f, transformedGeometry)
      })
    }

    /**
     * Reprojects the geometries in this SpatialRDD to the target CRS.
     * @param targetCRS the target coordinate reference system
     * @return a new RDD after geometries are transformed
     */
    def reproject(targetCRS: CoordinateReferenceSystem): SpatialRDD = {
      // For efficiency, the MathTransform is created once here and serialized to all worker nodes
      // We assume that all geometries in the RDD have the same SRID which is generally true
      val sourceCRS = CRSServer.sridToCRS(rdd.first().getGeometry.getSRID, rdd.sparkContext.getConf)
      rdd.reproject(sourceCRS, targetCRS)
    }
  }

  /**
   * Additional for arrays of features
   * @param records an array of features
   */
  implicit class TraversableSpatialFunctions(records: Array[_ <: IFeature]) {
    def mbr : EnvelopeNDLite = {
      val iterator = records.iterator
      if (!iterator.hasNext)
        return new EnvelopeNDLite(2)
      val firstRecord = iterator.next()
      val allMBR = new EnvelopeNDLite(GeometryHelper.getCoordinateDimension(firstRecord.getGeometry))
      allMBR.merge(firstRecord.getGeometry)
      while (iterator.hasNext) {
        allMBR.merge(iterator.next().getGeometry)
      }
      allMBR
    }
  }
}

object CGOperationsMixin extends CGOperationsMixin