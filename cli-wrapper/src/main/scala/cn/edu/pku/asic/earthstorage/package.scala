package cn.edu.pku.asic

import cn.edu.pku.asic.earthstorage.common.cg.SpatialDataTypesMixin
import cn.edu.pku.asic.earthstorage.common.io.ReadWriteMixin
import cn.edu.pku.asic.earthstorage.common.operations.SpatialOperationsMixin
import cn.edu.pku.asic.earthstorage.indexing.IndexMixin
import org.apache.spark.earthStorage.SparkSQLRegistration
import org.apache.spark.sql.SparkSession

/**
 * Contains implicit conversions that simplify the access to spatial functions in Beast.
 * To use it, add the following command to your Scala program.
 *
 * import edu.ucr.cs.bdlab.beast._
 */
/*通过隐式类，为Spark增加SpatialRDD等类，并为SparkContext、RDD等类增加新功能*/

package object earthStorage extends ReadWriteMixin  //各种格式编解码
  with SpatialOperationsMixin                       //空间运算
  with SpatialDataTypesMixin                        //新类型定义
  with IndexMixin                                   //索引功能混入
  //  with VisualizationMixin
  //  with RaptorMixin
{
  if (SparkSession.getActiveSession.nonEmpty) {
    SparkSQLRegistration.registerUDT
    SparkSQLRegistration.registerUDF(SparkSession.getActiveSession.get)
  }
}
