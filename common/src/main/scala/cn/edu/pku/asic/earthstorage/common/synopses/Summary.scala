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
package cn.edu.pku.asic.earthstorage.common.synopses

import cn.edu.pku.asic.earthstorage.common.cg.SpatialDataTypes.{JavaSpatialRDD, SpatialRDD}
import cn.edu.pku.asic.earthstorage.common.geolite.{EnvelopeND, EnvelopeNDLite, GeometryHelper, GeometryType, IFeature}
import com.fasterxml.jackson.core.JsonGenerator
import org.apache.spark.rdd.RDD
import org.apache.spark.{FutureAction, SparkContext}
import org.locationtech.jts.geom.Geometry

import java.io.{Externalizable, IOException, ObjectInput, ObjectOutput}
import scala.concurrent.duration.Duration
import scala.concurrent.{CanAwait, ExecutionContext, Future}
import scala.util.Try

/**
 * A fixed-size summary for a set of features. It includes the following aggregate values:
 *  - *size*: The estimated size of all features in bytes
 *  - *numFeatures*: Total number of features, i.e., records.
 *  - *numPoints*: The sum of number of points for all geometries in features.
 *  - *numNonEmptyGeometries*: Number of features with non-empty geometries
 *  - *sumSideLength*: The sum of all size lengths of all non-empty geometries. Can be combined with
 *    `numNonEmptyGeometries` to compute the average size per record.
 *  - *geometryType*: The geometry type of all non-empty geometries. This is the least inclusive type of all
 *    geometries. For example, if all geometries are points, the geometry type will be point. If it contains a mix
 *    of points and multi-points, a multipoint type is returned. If it contains a mix of points and polygons,
 *    a GeometryCollection type is used.
 */
class Summary extends EnvelopeNDLite with Externalizable {

  /** Total estimated storage size for the features */
  var size: Long = 0L

  /** Total number of features in the data file */
  var numFeatures: Long = 0L

  /** Total number of points in the geometries. e.g., linestrings and polygons can contribute more than one point */
  var numPoints: Long = 0L

  /** Count the number of non-empty geometries. This is used to calculate the average side length correctly */
  var numNonEmptyGeometries: Long = 0L

  /** The sum of side length along each dimension. Can be used to find the average side length. */
  var sumSideLength: Array[Double] = _

  /** The type of geometry stored in the dataset. */
  var geometryType = GeometryType.EMPTY

  /**
   * Copy constructor
   * @param other another summary to copy from
   */
  def this(other: Summary) {
    this()
    super.set(other)
    this.size = other.size
    this.numFeatures = other.numFeatures
    this.numPoints = other.numPoints
    this.numNonEmptyGeometries = other.numNonEmptyGeometries
    if (other.sumSideLength != null)
      this.sumSideLength = Array(other.sumSideLength: _*)
    this.geometryType = other.geometryType
  }

  protected def setNumPoints(p: Long): Unit = this.numPoints = p
  protected def setNumFeatures(f: Long): Unit = this.numFeatures = f
  protected def setSize(s: Long): Unit = this.size = s

  def incrementNumFeatures(inc: Long = 1): Unit = this.numFeatures += inc

  override def setCoordinateDimension(numDimensions: Int): Unit = {
    super.setCoordinateDimension(numDimensions)
    if (this.sumSideLength == null || this.sumSideLength.length != numDimensions)
      this.sumSideLength = new Array[Double](numDimensions)
  }

  /**
   * Expands the summary to enclose another partial summary. This can be used for incremental computation relying on
   * the associative and commutative properties on the attributes in the summary.
   * @param other the other summary to expand to
   * @return this summary so that the method can be called serially
   */
  def expandToSummary(other: Summary): Summary = {
    if (other.isEmpty())
      return this;
    assert(this.getCoordinateDimension == other.getCoordinateDimension,
      s"Incompatible number of dimensions ${this.getCoordinateDimension} != ${other.getCoordinateDimension}")
    merge(other)
    this.numFeatures += other.numFeatures
    this.size += other.size
    this.numPoints += other.numPoints
    for (d <- this.sumSideLength.indices)
      this.sumSideLength(d) += other.sumSideLength(d)
    this.numNonEmptyGeometries += other.numNonEmptyGeometries
    mergeWithGeometryType(other.geometryType)
    this
  }

  /**
   * Update the geometry type based on the given one.
   * @param otherGeometryType the geometry type to include in this summary
   */
  protected def mergeWithGeometryType(otherGeometryType: GeometryType): Unit =
    this.geometryType = this.geometryType.coerce(otherGeometryType)

  /**
   * Expands the summary to enclose the given geometry and update the size with the given size
   *
   * @param geometry a geometry to include in this summary
   * @param size     the size of the given geometry in bytes
   */
  def expandToGeometryWithSize(geometry: Geometry, size: Int): Unit = {
    merge(geometry)
    this.numFeatures += 1
    this.size += size
    this.mergeWithGeometryType(Summary.GeometryNameToType(geometry.getGeometryType))
    if (geometry.isEmpty) return
    this.numNonEmptyGeometries += 1
    this.numPoints += geometry.getNumPoints
    val e = new EnvelopeND(geometry.getFactory, this.getCoordinateDimension)
    e.merge(geometry)
    for (d <- sumSideLength.indices)
      sumSideLength(d) += e.getSideLength(d)
  }

  /**
   * Expands the summary to enclose the given geometry
   * @param geom the geometry to include in this summary
   */
  def expandToGeometry(geom: Geometry): Unit = expandToGeometryWithSize(geom, GeometryHelper.getGeometryStorageSize(geom))

  /**
   * Expands the summary to enclose the given feature. In addition to what {@link #expandToGeometry ( Geometry )} does,
   * it increase the total size to account for the non-geometric attributes in this feature.
   *
   * @param f the feature to include in this summary
   */
  def expandToFeature(f: IFeature): Unit = expandToGeometryWithSize(f.getGeometry, f.getStorageSize)

  override def setEmpty(): Unit = {
    super.setEmpty()
    this.numFeatures = 0
    this.size = 0
    this.numPoints = 0
    this.numNonEmptyGeometries = 0
    this.geometryType = GeometryType.EMPTY
  }

  @throws[IOException]
  override def writeExternal(out: ObjectOutput): Unit = {
    GeometryHelper.writeIEnvelope(this, out)
    out.writeLong(size)
    out.writeLong(numFeatures)
    out.writeLong(numPoints)
    out.writeLong(numNonEmptyGeometries)
    out.writeInt(geometryType.ordinal)
    for ($d <- 0 until getCoordinateDimension) {
      out.writeDouble(sumSideLength($d))
    }
  }

  @throws[IOException]
  override def readExternal(in: ObjectInput): Unit = {
    GeometryHelper.readIEnvelope(this, in)
    size = in.readLong
    numFeatures = in.readLong
    numPoints = in.readLong
    numNonEmptyGeometries = in.readLong
    geometryType = GeometryType.values()(in.readInt)
    for ($d <- 0 until getCoordinateDimension)
      sumSideLength($d) = in.readDouble
  }

  override def equals(other: Any): Boolean = {
    if (!super.equals(other)) return false
    val that = other.asInstanceOf[Summary]
    if (this.numFeatures != that.numFeatures) return false
    if (this.size != that.size) return false
    if (this.numPoints != that.numPoints) return false
    if (this.numNonEmptyGeometries != that.numNonEmptyGeometries) return false
    if (this.geometryType ne that.geometryType) return false
    for ($d <- this.sumSideLength.indices) {
      if (this.sumSideLength($d) != that.sumSideLength($d)) return false
    }
    true
  }

  def averageSideLength(i: Int): Double = sumSideLength(i) / numNonEmptyGeometries

  override def toString: String = {
    val str = new StringBuilder
    // Write MBR
    str.append("MBR: [(")
    for ($d <- 0 until getCoordinateDimension) {
      if ($d != 0) str.append(", ")
      str.append(getMinCoord($d))
    }
    str.append("), (")
    for ($d <- 0 until getCoordinateDimension) {
      if ($d != 0) str.append(", ")
      str.append(getMaxCoord($d))
    }
    str.append(")]")
    str.append(", size: ")
    str.append(size)
    str.append(", numFeatures: ")
    str.append(numFeatures)
    str.append(", numPoints: ")
    str.append(numPoints)
    str.append(", avgSideLength: [")
    for ($d <- 0 until getCoordinateDimension) {
      if ($d != 0) str.append(", ")
      str.append(averageSideLength($d))
    }
    str.append("]")
    str.toString
  }
}

object Summary {
  /**
   * A dictionary that maps each geometry name to a GeometryType instance
   */
  val GeometryNameToType: Map[String, GeometryType] = GeometryType.values.map(gt => gt.typename -> gt).toMap + ("LinearRing" -> GeometryType.LINESTRING)

  /**
   * Writes the summary of the dataset in JSON format
   *
   * @param jsonGenerator the generate to write the output to
   * @param summary       the summary of the dataset
   * @param feature       a sample feature to get the attribute names and types
   * @throws IOException if an error happens while writing the output
   */
  @throws[IOException]
  def writeSummaryWithSchema(jsonGenerator: JsonGenerator, summary: Summary, feature: IFeature): Unit = {
    jsonGenerator.writeStartObject()
    // Write MBR
    jsonGenerator.writeFieldName("extent")
    jsonGenerator.writeStartArray()
    jsonGenerator.writeNumber(summary.getMinCoord(0))
    jsonGenerator.writeNumber(summary.getMinCoord(1))
    jsonGenerator.writeNumber(summary.getMaxCoord(0))
    jsonGenerator.writeNumber(summary.getMaxCoord(1))
    jsonGenerator.writeEndArray()
    // Write sizes and number of records
    jsonGenerator.writeNumberField("size", summary.size)
    jsonGenerator.writeNumberField("num_features", summary.numFeatures)
    jsonGenerator.writeNumberField("num_non_empty_features", summary.numNonEmptyGeometries)
    jsonGenerator.writeNumberField("num_points", summary.numPoints)
    // Write average side length
    jsonGenerator.writeFieldName("avg_sidelength")
    jsonGenerator.writeStartArray()
    jsonGenerator.writeNumber(summary.averageSideLength(0))
    jsonGenerator.writeNumber(summary.averageSideLength(1))
    jsonGenerator.writeEndArray()
    // Write Coordinate Reference System (CRS)
    jsonGenerator.writeFieldName("srid")
    jsonGenerator.writeNumber(feature.getGeometry.getSRID)
    // Write attributes
    if (feature.getNumAttributes > 0) {
      jsonGenerator.writeFieldName("attributes")
      jsonGenerator.writeStartArray()
      for (i <- 0 until feature.getNumAttributes) {
        jsonGenerator.writeStartObject()
        var fieldName = feature.getAttributeName(i)
        if (fieldName == null || fieldName.length == 0) fieldName = "attribute#" + i
        jsonGenerator.writeStringField("name", fieldName)
        val attValue = feature.getAttributeValue(i)
        if (attValue != null) {
          val fieldType = attValue match {
            case _: String => "string"
            case _: Integer | _: java.lang.Long => "integer"
            case _: java.lang.Float | _: java.lang.Double => "number"
            case _: java.lang.Boolean => "boolean"
            case _ => "unknown"
          }
          jsonGenerator.writeStringField("type", fieldType)
        }
        jsonGenerator.writeEndObject() // End of attribute type

      }
      jsonGenerator.writeEndArray() // End of attributes

    }
    jsonGenerator.writeEndObject() // End of root object
  }

  /**
   * Compute partial summaries for the given [[SpatialRDD]].
   * It returns a set of Summaries in an RDD that can be combined to produce the final summary.
   * Typically, it returns one summary per partition in the input. If the number of partitions is
   * larger than 1,000, they are further coalesced into 1,000 partitions to reduce memory requirement
   * for a subsequent reduce action.
   * @param features the set of features to summarize
   * @param sizeFunction the function used to estimate the size of each feature
   * @return an RDD of summaries
   */
  private def computePartialSummaries(features: SpatialRDD, sizeFunction: IFeature => Int = f => f.getStorageSize): RDD[Summary] = {
    var partialSummaries = features.mapPartitions(iFeatures => {
      val summary = new Summary
      while (iFeatures.hasNext) {
        val feature = iFeatures.next()
        if (summary.getCoordinateDimension == 0)
          summary.setCoordinateDimension(GeometryHelper.getCoordinateDimension(feature.getGeometry))
        summary.expandToGeometryWithSize(feature.getGeometry, sizeFunction(feature))
      }
      Option(summary).iterator
    })
    // If the number of partitions is very large, Spark would fail because it will try to collect all
    // the partial summaries locally and would run out of memory in such case
    if (partialSummaries.getNumPartitions > 1000)
      partialSummaries = partialSummaries.coalesce(1000)
    partialSummaries
  }

  /**
   * Compute summary for features asynchronously
   * @param features the features to compute for
   * @param sizeFunction the size function to use
   * @return a future action for the summary
   */
  def computeForFeaturesAsync(features: SpatialRDD, sizeFunction: IFeature => Int = f => f.getStorageSize)
  : FutureAction[Summary] = {
    val partialSummaries: RDD[Summary] = computePartialSummaries(features, sizeFunction)
    val partialSummariesAsync: FutureAction[Seq[Summary]] = partialSummaries.collectAsync()
    new FutureAction[Summary] {
      var _finalResult: Summary = _
      private def getResult(theresult: Seq[Summary]): Summary = {
        if (_finalResult == null)
          _finalResult = theresult.reduce((mbr1, mbr2) => mbr1.expandToSummary(mbr2))
        _finalResult
      }

      override def cancel(): Unit = partialSummariesAsync.cancel()

      override def ready(atMost: Duration)(implicit permit: CanAwait): this.type = {
        partialSummariesAsync.ready(atMost)
        this
      }

      override def result(atMost: Duration)(implicit permit: CanAwait): Summary = {
        val theResult: Seq[Summary] = partialSummariesAsync.result(atMost)
        if (theResult == null)
          return null
        getResult(theResult)
      }

      override def onComplete[U](func: Try[Summary] => U)(implicit executor: ExecutionContext): Unit = {
        partialSummariesAsync.onComplete(tryTheResult => {
          val theResult = tryTheResult.get
          val finalResult = getResult(theResult)
          func.apply(Try.apply(finalResult))
        })
      }

      override def isCompleted: Boolean = partialSummariesAsync.isCompleted

      override def isCancelled: Boolean = partialSummariesAsync.isCompleted

      override def value: Option[Try[Summary]] = {
        if (_finalResult != null)
          Option(Try.apply(_finalResult))
        else if (!partialSummariesAsync.isCompleted)
          Option.empty[Try[Summary]]
        else {
          val theResult: Seq[Summary] = partialSummariesAsync.value.get.get
          Option(Try.apply(getResult(theResult)))
        }

      }

      override def jobIds: Seq[Int] = partialSummariesAsync.jobIds

      override def transform[S](f: Try[Summary] => Try[S])(implicit executor: ExecutionContext): Future[S] = {
        partialSummariesAsync.transform(tryTheResult => {
          val theResult = tryTheResult.get
          val finalResult = getResult(theResult)
          f.apply(Try.apply(finalResult))
        })
      }

      override def transformWith[S](f: Try[Summary] => Future[S])(implicit executor: ExecutionContext): Future[S] = {
        partialSummariesAsync.transformWith(tryTheResult => {
          val theResult = tryTheResult.get
          val finalResult = getResult(theResult)
          f.apply(Try.apply(finalResult))
        })
      }
    }
  }

  /**
   * Compute the summary of a set of features with a custom function for estimating the size of a features
   * 用指定的大小计算函数，，估算输入数据集的计算
   * @param features the set of features to summarize
   * @param sizeFunction a function that returns the size for each features
   * @return the summary of the input features
   */
  def computeForFeatures(features: SpatialRDD, sizeFunction: IFeature => Int = f => f.getStorageSize) : Summary = {
    val partialSummaries: RDD[Summary] = computePartialSummaries(features, sizeFunction)
    partialSummaries.reduce((mbr1, mbr2) => mbr1.expandToSummary(mbr2))
  }

  /**
   * Compute the MBR, count, and size of a JavaRDD of geometries
   * @param features the JavaRDD of features to compute the summary for
   * @return a Summary for the given features
   */
  def computeForFeatures(features : JavaSpatialRDD) : Summary = computeForFeatures(features.rdd)

  /**
   * Compute the MBR, count, and size of a JavaRDD of geometries
   * @param features the JavaRDD of features to compute the summary for
   * @return a Summary for the given features
   */
  def computeForFeaturesWithSize(features : JavaSpatialRDD, sizeFunction: IFeature => Int) : Summary =
    computeForFeatures(features.rdd, sizeFunction)

  /**
   * Create a summary accumulator that uses the method [[IFeature#getStorageSize]] to accumulate the sizes of
   * the features.
   *
   * @param sc the spark context to register the accumulator to
   * @return the initialized and registered accumulator
   */
  def createSummaryAccumulator(sc: SparkContext) : SummaryAccumulator = createSummaryAccumulator(sc, _.getStorageSize)

  /**
   * Create a summary accumulator that uses the method [[IFeature#getStorageSize]] to accumulate the sizes of
   * the features.
   *
   * @param sc the spark context to register the accumulator to
   * @return the initialized and registered accumulator
   */
  def createSummaryAccumulator(sc: SparkContext, sizeFunction: IFeature => Int) : SummaryAccumulator = {
    val accumulator = new SummaryAccumulator(sizeFunction)
    sc.register(accumulator, "SummaryAccumulator")
    accumulator
  }
}