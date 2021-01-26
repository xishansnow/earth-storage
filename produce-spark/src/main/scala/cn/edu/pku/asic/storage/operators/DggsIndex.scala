package cn.edu.pku.asic.storage.operators

import cn.edu.pku.asic.storage.common.cg.SpatialPartitioner
import cn.edu.pku.asic.storage.common.cli.{AppOptions, CLIOperation}
import cn.edu.pku.asic.storage.indexing.indexing.{RGrovePartitioner, RSGrovePartitioner}
import cn.edu.pku.asic.storage.common.io.{SpatialFileRDD, SpatialOutputFormat}
import cn.edu.pku.asic.storage.common.utils.OperationMetadata
import cn.edu.pku.asic.storage.indexing.IndexHelper
import org.apache.spark.internal.Logging

import java.io.PrintStream
import java.util

/**
 * Builds an indexing over a set of features that can be either stored to disk or kept as an RDD.
 */
@OperationMetadata(
  shortName =  "indexing",
  description = "Builds a distributed spatial indexing",
  inputArity = "+",
  outputArity = "1",
  inheritParams = Array(classOf[SpatialFileRDD], classOf[SpatialOutputFormat], classOf[RGrovePartitioner], classOf[RSGrovePartitioner])
)
object DggsIndex extends CLIOperation with Logging{
  override def printUsage(out: PrintStream): Unit = {
    val partitioners: Map[String, Class[_ <: SpatialPartitioner]] = IndexHelper.partitioners
    out.println("The available indexes are:")
    partitioners.foreach(kv => {
      val indexerMetadata = kv._2.getAnnotation(classOf[SpatialPartitioner.Metadata])
      out.println(s"- ${kv._1}: ${indexerMetadata.description}")
    })
  }


  override def addDependentClasses(opts: AppOptions, classes: util.Stack[Class[_]]): Unit = {


  }

  override def setup(opts: AppOptions): Unit = {
    super.setup(opts)
  }
}
