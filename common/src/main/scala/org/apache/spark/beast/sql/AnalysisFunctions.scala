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

import cn.edu.pku.asic.storage.common.geolite.GeometryReader
import org.apache.spark.sql.catalyst.InternalRow
import org.apache.spark.sql.catalyst.analysis.TypeCheckResult
import org.apache.spark.sql.catalyst.expressions.Expression
import org.apache.spark.sql.catalyst.expressions.codegen.CodegenFallback
import org.apache.spark.sql.catalyst.util.ArrayData
import org.apache.spark.sql.types.DataType
import org.locationtech.jts.geom.{Geometry, GeometryFactory}

/**
 * An abstract analysis function that takes one or more geometries and returns a geometry.
 */
trait ST_AbstractAnalysis extends Expression with CodegenFallback {
  /**Geometry factory to create geometries*/
  val geometryFactory: GeometryFactory = GeometryReader.DefaultGeometryFactory

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

  override def dataType: DataType = GeometryDataType

  override def eval(input: InternalRow): Any = {
    val geometries: Seq[Geometry] = inputExpressions.map(e => GeometryDataType.getGeometryFromArray(e.eval(input).asInstanceOf[ArrayData]))
    val result: Geometry = performFunction(geometries)
    GeometryDataType.setGeometryInRow(result)
  }

  def performFunction(geometries: Seq[Geometry]): Geometry
}

/**
 * Computes the smallest convex Polygon of the input geometry
 * @param inputExpressions
 */
case class ST_ConvexHull(override val inputExpressions: Seq[Expression])
  extends ST_AbstractAnalysis {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).convexHull()
}

/**
 * Computes the boundary of the input geometry
 * @param inputExpressions
 */
case class ST_Boundary(override val inputExpressions: Seq[Expression])
  extends ST_AbstractAnalysis {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).getBoundary
}

/**
 * Creates a new Geometry which is a normalized copy of the input one
 * @param inputExpressions
 */
case class ST_Norm(override val inputExpressions: Seq[Expression])
  extends ST_AbstractAnalysis {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).norm()
}

/**
 * Computes a new geometry which has all component coordinate sequences in reverse order (opposite orientation)
 * to the input one
 * @param inputExpressions
 */
case class ST_Reverse(override val inputExpressions: Seq[Expression])
  extends ST_AbstractAnalysis {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).reverse()
}

/**
 * Computes the centroid of the input geometry
 * @param inputExpressions
 */
case class ST_Centroid(override val inputExpressions: Seq[Expression])
  extends ST_AbstractAnalysis {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).getCentroid
}
/**
 * Computes the intersection of two input geometries
 * @param inputExpressions
 */
case class ST_Intersection(override val inputExpressions: Seq[Expression])
    extends ST_AbstractAnalysis {
  override val inputArity: Int = 2

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).intersection(geometries(1))
}

/**
 * Computes a Geometry representing the point-set which is contained in two input Geometries.
 * @param inputExpressions
 */
case class ST_Union(override val inputExpressions: Seq[Expression])
  extends ST_AbstractAnalysis {
  override val inputArity: Int = 2

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).union(geometries(1))
}

/**
 * Computes a Geometry representing the closure of the point-set of the points contained in the first input Geometry
 * that are not contained in the second one
 * @param inputExpressions
 */
case class ST_Difference(override val inputExpressions: Seq[Expression])
  extends ST_AbstractAnalysis {
  override val inputArity: Int = 2

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).difference(geometries(1))
}

/**
 * Computes the envelope (bounding box) of the input geometry
 * @param inputExpressions
 */
case class ST_Envelope(override val inputExpressions: Seq[Expression])
  extends ST_AbstractAnalysis {
  override val inputArity: Int = 1

  override def performFunction(geometries: Seq[Geometry]): Geometry =
    geometries(0).getEnvelope
}
