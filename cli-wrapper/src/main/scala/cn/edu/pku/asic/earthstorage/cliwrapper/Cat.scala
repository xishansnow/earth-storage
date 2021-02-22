package cn.edu.pku.asic.earthstorage.cliwrapper

import cn.edu.pku.asic.earthStorage._
import cn.edu.pku.asic.earthstorage.common.cli.{AppOptions, CLIOperation, OperationMetadata, OperationParam}
import cn.edu.pku.asic.earthstorage.common.io.{SpatialFileRDD, SpatialOutputFormat}
import org.apache.spark.SparkContext
import org.apache.spark.internal.Logging

@OperationMetadata(
  shortName = "cat",
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
