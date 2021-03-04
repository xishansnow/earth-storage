package cn.edu.pku.asic.earthstorage.cliwrapper

import cn.edu.pku.asic.earthStorage._
import cn.edu.pku.asic.earthstorage.common.cg.SpatialDataTypes.{JavaSpatialRDD, SpatialRDD}
import cn.edu.pku.asic.earthstorage.common.cli.{AppOptions, CLIOperation, OperationMetadata, OperationParam}
import cn.edu.pku.asic.earthstorage.common.geolite.{EnvelopeNDLite, IFeature}
import cn.edu.pku.asic.earthstorage.common.io.{SpatialFileRDD, SpatialOutputFormat}
import cn.edu.pku.asic.earthstorage.common.operations.FeatureWriterSizeFunction
import cn.edu.pku.asic.earthstorage.common.synopses.{AbstractHistogram, HistogramOP, UniformHistogram}
import org.apache.spark.SparkContext
import org.apache.spark.internal.Logging

import java.io.IOException
import scala.annotation.varargs

/**
 * Computes the histogram of an input file.
 */
@OperationMetadata(
  shortName = "histogram",
  description = "Computes a uniform histogram for a geometry file",
  inputArity = "1",
  outputArity = "0",
  inheritParams = Array(classOf[SpatialFileRDD], classOf[SpatialOutputFormat])
)
object Histogram extends CLIOperation with Logging {

  /**
   * Number of buckets in the histogram
   */
  @OperationParam(
    description = "Total number of buckets in the histogram",
    defaultValue = "1000"
  )
  val NumBuckets = "numbuckets"

  @OperationParam(
    description = "Type of histogram {simple, euler}",
    defaultValue = "simple"
  )
  val HistogramType = "histogramtype"

  @OperationParam(
    description = "The value to compute in the histogram {count, size, writesize}",
    defaultValue = "count"
  )
  val HistogramValue = "histogramvalue"

  @OperationParam(
    description = "Method to compute the histogram {onepass, onehalfpass, twopasses, sparse}",
    defaultValue = "twopasses"
  )
  val ComputationMethod = "method"

  @varargs def computePointWriteSizeHistogram(features: SpatialRDD, mbb: EnvelopeNDLite, opts: AppOptions, numBuckets: Int*):
  UniformHistogram = {
    val sizeFunction = new FeatureWriterSizeFunction(opts)
    HistogramOP.computeHistogram(features, sizeFunction, HistogramOP.TwoPass, HistogramOP.PointHistogram, numBuckets: _*).asInstanceOf[UniformHistogram]
  }

  @varargs def computePointWriteSizeHistogram(features: JavaSpatialRDD, mbb: EnvelopeNDLite, opts: AppOptions, numBuckets: Int*):
  UniformHistogram =
    computePointWriteSizeHistogram(features.rdd, mbb, opts, numBuckets: _*)

  @throws(classOf[IOException])
  override def run(opts: AppOptions, inputs: Array[String], outputs: Array[String], sc: SparkContext): AbstractHistogram = {

    // Extract user parameters
    val numBuckets = opts.getInt(NumBuckets, 1000)
    val bo: AppOptions = opts
    val features = sc.spatialFile(inputs(0), opts)

    //获取计算直方图计数的函数，三种：统计要素个数、统计要素内存容量、统计要素输出时的物理存储容量，默认是统计要素个数
    val sizeFunction: IFeature => Int = opts.getString(HistogramValue, "count").toLowerCase match {
      case "count" => _ => 1 // 如果采用count计数方式，则计数函数返回1
      case "size" => _.getStorageSize // 如果采用size技数方式，则使用getStorageSize函数，技术函数返回要素占用的容量大小
      case "writesize" => new FeatureWriterSizeFunction(opts) // 如果采用writesize计数方式，则使用FeatureWriterSizeFunction函数，计数函数返回根据输出格式自动计算的容量大小
    }

    //获取直方图计算的方法，四种: onepass，onehalfpass，twopass，sparse，默认是twopass
    val method: HistogramOP.ComputationMethod = opts.getString(ComputationMethod, "twopasses").toLowerCase() match {
      case "twopasses" => HistogramOP.TwoPass
      case "onehalfpass" => HistogramOP.OneHalfPass
      case "onepass" => HistogramOP.OnePass
      case "sparse" => HistogramOP.Sparse
    }

    //获取直方图类型，两种：simple，euler
    val htype: HistogramOP.HistogramType = opts.getString(HistogramType, "simple").toLowerCase match {
      case "simple" => HistogramOP.PointHistogram
      case "euler" => HistogramOP.EulerHistogram
    }

    //根据参数计算直方图，注意：这里是根据数据集最小外包举行的大小，自动生成的均匀网格，各维度网格数量为numBuckets（默认为1000个网格）
    HistogramOP.computeHistogram(features, sizeFunction, method, htype, numBuckets)
  }
}
