package cn.edu.pku.coe.asic.earthstorage.operators

import edu.ucr.cs.bdlab.beast.cg.SpatialPartitioner
import edu.ucr.cs.bdlab.beast.common.{BeastOptions, CLIOperation}
import edu.ucr.cs.bdlab.beast.indexing.{IndexHelper, RGrovePartitioner, RSGrovePartitioner}
import edu.ucr.cs.bdlab.beast.io.{SpatialFileRDD, SpatialOutputFormat}
import edu.ucr.cs.bdlab.beast.util.OperationMetadata
import org.apache.spark.internal.Logging

import java.io.PrintStream
import java.util

/**
 * Builds an index over a set of features that can be either stored to disk or kept as an RDD.
 */
@OperationMetadata(
  shortName =  "index",
  description = "Builds a distributed spatial index",
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


  override def addDependentClasses(opts: BeastOptions, classes: util.Stack[Class[_]]): Unit = {


  }

  override def setup(opts: BeastOptions): Unit = {
    super.setup(opts)
  }
}
