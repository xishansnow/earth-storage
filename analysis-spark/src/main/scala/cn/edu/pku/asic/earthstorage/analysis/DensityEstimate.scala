package cn.edu.pku.asic.earthstorage.analysis

import cn.edu.pku.asic.earthstorage.common.cli.{AppOptions, CLIOperation, OperationMetadata}
import cn.edu.pku.asic.earthstorage.common.io.{SpatialFileRDD, SpatialOutputFormat}
import org.apache.spark.SparkContext
import org.apache.spark.internal.Logging

@OperationMetadata(
  shortName = "sj",
  description = "Computes density that finds all overlapping features from two files.",
  inputArity = "2",
  outputArity = "?",
  inheritParams = Array(classOf[SpatialFileRDD], classOf[SpatialOutputFormat])
)
object DensityEstimate extends CLIOperation with Logging{
  /**
   * Run the main function using the given user command-line options and spark context
   *
   * @param opts    user options for configuring the operation
   * @param inputs  inputs provided by the user
   * @param outputs outputs provided by the user
   * @param sc      the Spark context used to run the operation
   * @return an optional result of this operation
   */
  override def run(opts: AppOptions, inputs: Array[String], outputs: Array[String], sc: SparkContext): Any = ???
}
