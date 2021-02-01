/*
 * Copyright 2018 University of California, Riverside
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
package cn.edu.pku.asic.storage.common.io

import cn.edu.pku.asic.storage.common.cli.{AppOptions, CLIOperation}
import cn.edu.pku.asic.storage.common.io.ReadWriteMixin._
import cn.edu.pku.asic.storage.common.utils.{OperationMetadata, OperationParam}
import org.apache.spark.SparkContext
import org.apache.spark.internal.Logging


@OperationMetadata(
  shortName =  "cat",
  description = "Writes a file to the output. Used to convert file format",
  inputArity = "1",
  outputArity = "1",
  inheritParams = Array(classOf[SpatialFileRDD], classOf[SpatialOutputFormat])
)
object Cat extends CLIOperation with Logging {

  @OperationParam(description = "Number of partitions in the output", defaultValue = "Same as input")
  val NumPartitions = "numpartitions"

  override def run(opts: AppOptions, inputs: Array[String], outputs: Array[String], sc: SparkContext): Any = {
    // Read the input features and create pairs to be able to use Hadoop output format
    var features = sc.spatialFile(inputs(0), opts)
    // Store to the output
    val oFormat = opts.getString(SpatialOutputFormat.OutputFormat, opts.getString(SpatialFileRDD.InputFormat))
    val outputPartitions = opts.get(NumPartitions)
    if (outputPartitions.isDefined)
      features = features.coalesce(outputPartitions.get.toInt)
    features.writeSpatialFile(outputs(0), oFormat, opts)
  }
}
