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
package org.apache.spark.beast.sql

import cn.edu.pku.asic.storage.common.geolite.{GeometryReader, PointND}
import org.apache.spark.sql.catalyst.InternalRow
import org.apache.spark.sql.catalyst.analysis.TypeCheckResult
import org.apache.spark.sql.catalyst.expressions.Expression
import org.apache.spark.sql.catalyst.expressions.codegen.CodegenFallback
import org.apache.spark.sql.types.DataType
import org.locationtech.jts.geom._
import org.locationtech.jts.io.{WKBReader, WKTReader}

trait AbstractCreationFunction extends Expression with CodegenFallback {
  /**Geometry factory to create geometries*/
  val geometryFactory: GeometryFactory = GeometryReader.DefaultGeometryFactory

  val inputExpressions: Seq[Expression]

  val inputArity: Int

  override def checkInputDataTypes(): TypeCheckResult = {
    if (inputArity == -1 || inputExpressions.length == inputArity)
      TypeCheckResult.TypeCheckSuccess
    else
      TypeCheckResult.TypeCheckFailure(s"Function $prettyName expects $inputArity inputs but received ${inputExpressions.length}")
  }

  override def children: Seq[Expression] = inputExpressions

  override def nullable: Boolean = false

  override def foldable: Boolean = inputExpressions.forall(_.foldable)

  override def dataType: DataType = GeometryDataType

  override def eval(input: InternalRow): Any = {
    val inputValues: Seq[Any] = inputExpressions.map(e => e.eval(input))
    val result: Geometry = performFunction(inputValues)
    GeometryDataType.setGeometryInRow(result)
  }

  def performFunction(inputs: Seq[Any]): Geometry
}

/**
 * Creates a point from two or more coordinates.
 * @param inputExpressions
 */
case class ST_CreatePoint(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {
  override val inputArity: Int = -1

  override def performFunction(inputs: Seq[Any]): Geometry = {
    val coordinates: Array[Double] = inputs.map(x => x.asInstanceOf[Number].doubleValue()).toArray
    if (coordinates == 2)
      geometryFactory.createPoint(new CoordinateXY(coordinates(0), coordinates(1)))
    else if (coordinates == 3)
      geometryFactory.createPoint(new CoordinateXYZM(coordinates(0), coordinates(1), coordinates(2), Double.NaN))
    else
      new PointND(geometryFactory, coordinates.length, coordinates:_*)
  }
}

/**
 * Creates a line from multi coordinates.
 * @param inputExpressions
 */
case class ST_CreateLine(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {
  override val inputArity: Int = -1

  override def performFunction(inputs: Seq[Any]): Geometry = {
    val coordinates: Array[Double] = inputs.map(x => x.asInstanceOf[Number].doubleValue()).toArray
    if (coordinates.length % 2 != 0) {
      throw new IllegalArgumentException("The number of input values is illegal")
    }
    val coordArray: Array[Coordinate] = new Array[Coordinate](coordinates.length / 2)
    for (i <- coordinates.indices by 2) {
      coordArray(i) = new Coordinate(coordinates(i), coordinates(i + 1))
    }
    geometryFactory.createLineString(coordArray)
  }
}

/**
 * Creates a polygon from multi coordinates.
 * @param inputExpressions
 */
case class ST_CreatePolygon(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {
  override val inputArity: Int = -1

  override def performFunction(inputs: Seq[Any]): Geometry = {
    val coordinates: Array[Double] = inputs.map(x => x.asInstanceOf[Number].doubleValue()).toArray
    val len: Int = coordinates.length
    if (len < 8 || len % 2 != 0) {
      throw new IllegalArgumentException("The number of input values is illegal")
    }
    if (coordinates(0) != coordinates(len - 2) || coordinates(1) != coordinates(len - 1)) {
      throw new IllegalArgumentException("Points of LinearRing do not form a closed linestring")
    }
    val coordArray: Array[Coordinate] = new Array[Coordinate](coordinates.length / 2)
    for (i <- coordinates.indices by 2) {
      coordArray(i / 2) = new Coordinate(coordinates(i), coordinates(i + 1))
    }
    geometryFactory.createPolygon(coordArray)
  }
}

/**
 * Creates a box (rectangle) from the coordinates of the two corners
 * @param inputExpressions
 */
case class ST_CreateBox(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {
  override val inputArity: Int = 4

  override def performFunction(inputs: Seq[Any]): Geometry = {
    val coordinates = inputs.map(x => x.asInstanceOf[Number].doubleValue()).toArray
    val x1 = coordinates(0)
    val y1 = coordinates(1)
    val x2 = coordinates(2)
    val y2 = coordinates(3)
    val corners: Array[Coordinate] = new Array[Coordinate](5)
    corners(0) = new Coordinate(x1, y1)
    corners(1) = new Coordinate(x1, y2)
    corners(2) = new Coordinate(x2, y2)
    corners(3) = new Coordinate(x2, y1)
    corners(4) = corners(0)
    geometryFactory.createPolygon(corners)
  }
}

/**
 * Takes a list of point locations and IDs and creates either a linestring or
 * polygon based on whether the last point is the same as the first point or not
 * @param inputExpressions
 */
case class ST_CreateLinePolygon(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {
  override val inputArity: Int = -1

  override def performFunction(inputs: Seq[Any]): Geometry = {
    val length = inputs.length
    if (length < 3 || length % 3 != 0) {
      throw new IllegalArgumentException("The number of input values is illegal")
    }

    val pointIDs = new Array[Long](length / 3)
    val coordinates = new Array[Double](length / 3 * 2)
    for (i <- 0 to inputs.length) {
      if (i < inputs.length / 3) {
        pointIDs(i) = inputs(i).asInstanceOf[Number].longValue()
      } else {
        coordinates(i - length / 3) = inputs(i).asInstanceOf[Number].doubleValue()
      }
    }

    val is_polygon = pointIDs(0) == pointIDs(pointIDs.length - 1)
    val coords = new Array[Coordinate](coordinates.length / 2)
    for (i <- coordinates.indices by 2) {
      if (i == coordinates.length - 2 && is_polygon) {
        coords(i / 2) = coords(0)
      } else {
        coords(i / 2) = new Coordinate(coordinates(i), coordinates(i + 1))
      }
    }

    if (is_polygon) {
      if (coords.length == 1 || coords.length == 2) {
        geometryFactory.createPoint(coords(0))
      } else if (coords.length == 3) {
        geometryFactory.createLineString(coords.take(2))
      } else {
        geometryFactory.createPolygon(coords)
      }
    } else {
      geometryFactory.createLineString(coords)
    }
  }
}

/**
 * Creates a two-dimensional envelope from four coordinates (x1, y1, x2, y2)
 * @param inputExpressions
 */
case class ST_CreateEnvelope(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {

  override val inputArity: Int = 4

  override def performFunction(inputs: Seq[Any]): Geometry = {
    val coordinates: Array[Double] = inputs.map(x => x.asInstanceOf[Number].doubleValue()).toArray
    geometryFactory.toGeometry(new Envelope(coordinates(0), coordinates(2), coordinates(1), coordinates(3)))
  }
}

/**
 * Create a geometry from a WKT representation
 * @param inputExpressions
 */
case class ST_FromWKT(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {

  override val inputArity: Int = 1

  val wktParser = new WKTReader(geometryFactory)

  override def performFunction(inputs: Seq[Any]): Geometry = {
    wktParser.read(inputs.head.asInstanceOf[String])
  }
}

/**
 * Create a geometry from a WKB representation
 * @param inputExpressions
 */
case class ST_FromWKB(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {

  override val inputArity: Int = 1

  val wkbParser = new WKBReader(geometryFactory)

  override def performFunction(inputs: Seq[Any]): Geometry = {
    wkbParser.read(inputs.head.asInstanceOf[Array[Byte]])
  }
}

/**
 * Finds the minimal bounding rectangle (MBR) of a set of shapes
 * @param inputExpressions
 */
case class ST_Extent(override val inputExpressions: Seq[Expression])
  extends AbstractCreationFunction {

  override val inputArity: Int = -1

  override def performFunction(inputs: Seq[Any]): Geometry = {
    val geometries = inputs.map(x => x.asInstanceOf[Geometry]).toArray
    val geometryCollection = geometryFactory.createGeometryCollection(geometries)
    geometryCollection.getEnvelope
  }
}