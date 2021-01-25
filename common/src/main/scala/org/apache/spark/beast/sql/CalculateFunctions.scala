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

import org.apache.spark.sql.catalyst.InternalRow
import org.apache.spark.sql.catalyst.analysis.TypeCheckResult
import org.apache.spark.sql.catalyst.expressions.Expression
import org.apache.spark.sql.catalyst.expressions.codegen.CodegenFallback
import org.apache.spark.sql.catalyst.util.ArrayData
import org.apache.spark.sql.types.{DataType, DoubleType}
import org.locationtech.jts.geom.{Coordinate, Geometry}

/**
 * An abstract analysis function that takes one or more geometries and returns a number.
 */
trait ST_AbstractCalculate extends Expression with CodegenFallback {
  val inputExpressions: Seq[Expression]

  val inputArity: Int

  override def children: Seq[Expression] = inputExpressions

  override def checkInputDataTypes(): TypeCheckResult = {
    if (inputExpressions.length == inputArity && inputExpressions.forall(_.dataType == GeometryDataType))
      TypeCheckResult.TypeCheckSuccess
    else
      TypeCheckResult.TypeCheckFailure(s"Function $prettyName expects $inputArity geometry arguments")
  }

  override def nullable: Boolean = false

  override def foldable: Boolean = inputExpressions.forall(_.foldable)

  override def dataType: DataType = DoubleType

  override def eval(input: InternalRow): Double = {
    val geometries: Seq[Geometry] = inputExpressions.map(e => GeometryDataType.getGeometryFromArray(e.eval(input).asInstanceOf[ArrayData]))
    performFunction(geometries)
  }

  def performFunction(geometries: Seq[Geometry]): Double
}

/**
 * Computes the area of the input geometry
 * @param inputExpressions
 */
case class ST_Area(override val inputExpressions: Seq[Expression])
  extends ST_AbstractCalculate {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Double =
    geometries(0).getArea
}

/**
 * Computes the left bound of the input geometry
 * @param inputExpressions
 */
case class ST_XMin(override val inputExpressions: Seq[Expression])
  extends ST_AbstractCalculate {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Double = {
    val coords : Array[Coordinate] = geometries(0).getEnvelope.getCoordinates
    if (coords.length == 0) throw new IllegalArgumentException("ST_XMin cannot work on empty geometries")
    if (coords.length == 1) coords(0).x
    else if (coords.length == 2) Math.min(coords(0).x, coords(1).x)
    else Math.min(coords(0).x, coords(2).x)
  }
}

/**
 * Computes the right bound of the input geometry
 * @param inputExpressions
 */
case class ST_XMax(override val inputExpressions: Seq[Expression])
  extends ST_AbstractCalculate {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Double = {
    val coords : Array[Coordinate] = geometries(0).getEnvelope.getCoordinates
    if (coords.length == 0) throw new IllegalArgumentException("ST_XMax cannot work on empty geometries")
    if (coords.length == 1) coords(0).x
    else if (coords.length == 2) Math.max(coords(0).x, coords(1).x)
    else Math.max(coords(0).x, coords(2).x)
  }
}

/**
 * Computes the lower bound of the input geometry
 * @param inputExpressions
 */
case class ST_YMin(override val inputExpressions: Seq[Expression])
  extends ST_AbstractCalculate {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Double = {
    val coords : Array[Coordinate] = geometries(0).getEnvelope.getCoordinates
    if (coords.length == 0) throw new IllegalArgumentException("ST_YMin cannot work on empty geometries")
    if (coords.length == 1) coords(0).y
    else if (coords.length == 2) Math.min(coords(0).y, coords(1).y)
    else Math.min(coords(0).y, coords(2).y)
  }
}

/**
 * Computes the upper bound of the input geometry
 * @param inputExpressions
 */
case class ST_YMax(override val inputExpressions: Seq[Expression])
  extends ST_AbstractCalculate {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Double = {
    val coords : Array[Coordinate] = geometries(0).getEnvelope.getCoordinates
    if (coords.length == 0) throw new IllegalArgumentException("ST_YMax cannot work on empty geometries")
    if (coords.length == 1) coords(0).y
    else if (coords.length == 2) Math.max(coords(0).y, coords(1).y)
    else Math.max(coords(0).y, coords(2).y)
  }
}

/**
 * Computes the number of points of the input geometry
 * @param inputExpressions
 */
case class ST_NumPoints(override val inputExpressions: Seq[Expression])
  extends ST_AbstractCalculate {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Double = {
    geometries(0).getNumPoints
  }
}
