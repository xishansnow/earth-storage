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
package cn.edu.pku.asic.storage.common.cg

import cn.edu.pku.asic.storage.common.cg.SpatialDataTypes.SpatialRDD
import cn.edu.pku.asic.storage.common.geolite.Feature
import org.apache.spark.beast.sql.GeometryDataType
import org.apache.spark.sql.Row
import org.apache.spark.sql.types.{StructField, StructType}
import org.locationtech.jts.geom.{Envelope, Geometry}
import scala.+:
import scala.collection.mutable

/**
 * An iterator that takes a single geometry as input and produces a sequence of geometries by recursively breaking
 * the geometry into four pieces until a threshold is met in terms of number of points. In other words, all produced
 * geometries must have at most the threshold number of points.
 */
class GeometryQuadSplitter(geometry: Geometry, threshold: Int) extends Iterator[Geometry] {
  val toBreak: mutable.ArrayBuffer[Geometry] = new mutable.ArrayBuffer[Geometry]()
  toBreak.append(geometry)

  override def hasNext: Boolean = toBreak.nonEmpty

  override def next(): Geometry = {
    var g: Geometry = toBreak.remove(toBreak.length - 1)
    while (g != null && g.getNumPoints > threshold) {
      val e: Envelope = g.getEnvelopeInternal
      val centerX: Double = (e.getMinX + e.getMaxX) / 2
      val centerY: Double = (e.getMinY + e.getMaxY) / 2
      try {toBreak.append(g.intersection(g.getFactory.toGeometry(new Envelope(e.getMinX, centerX, e.getMinY, centerY))))}
      catch { case _: Exception => }
      try {toBreak.append(g.intersection(g.getFactory.toGeometry(new Envelope(centerX, e.getMaxX, e.getMinY, centerY))))}
      catch { case _: Exception => }
      try {toBreak.append(g.intersection(g.getFactory.toGeometry(new Envelope(e.getMinX, centerX, centerY, e.getMaxY))))}
      catch { case _: Exception => }
      try {g = g.intersection(g.getFactory.toGeometry(new Envelope(centerX, e.getMaxX, centerY, e.getMaxY)))}
      catch { case _: Exception => g = if (toBreak.isEmpty) null else toBreak.remove(toBreak.length - 1)}
    }
    g
  }
}

object GeometryQuadSplitter {
  /**
   * Splits all geometries in the given RDD[IFeature] into a new RDD[IFeature] where geometries are broken down
   * using the quad split partitioning approach. If [[keepBothGeometries]] is set to `true`, the resulting features
   * will contain both geometries where the broken geometry appears first. If [[keepBothGeometries]] is set to `false`,
   * on the simplified geometry is produced and the original geometry is removed.
   * @param spatialRDD the original rdd to split
   * @param threshold the quad split threshold
   * @param keepBothGeometries if `true`, both geometries will be kept in the output
   * @return
   */
  def splitRDD(spatialRDD: SpatialRDD, threshold: Int, keepBothGeometries: Boolean): SpatialRDD = {
    spatialRDD.mapPartitions(features => {
      features.flatMap(feature => {
        val smallGeometries: Iterator[Geometry] = new GeometryQuadSplitter(feature.getGeometry, threshold)
        if (keepBothGeometries) {
          smallGeometries.map(g => {
            // Keep both geometries
            val values: Seq[Any] = g +: Row.unapplySeq(feature).get
            val schema: Seq[StructField] = StructField("geometry", GeometryDataType) +: feature.schema
            new Feature(values.toArray, StructType(schema))
          })
        } else {
          smallGeometries.map(g => {
            // Replace the original geometry with the simplified one
            Feature.create(feature, g)
          })
        }
      })
    }, true)
  }
}