package cn.edu.pku.asic.earthstorage.compute

import cn.edu.pku.asic.earthstorage.common.cli.{AppOptions, CLIOperation}
import org.apache.spark.SparkContext
import org.apache.spark.internal.Logging

import java.io.PrintStream

object add extends CLIOperation with Logging{
  /**
   * Print the usage of this class (if any)
   *
   * @param out
   */
  override def printUsage(out: PrintStream): Unit = super.printUsage(out)

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
