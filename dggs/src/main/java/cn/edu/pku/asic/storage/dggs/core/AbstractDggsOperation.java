package cn.edu.pku.asic.storage.dggs.core;


import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

/**
 *=============================================================================================
 * @ProjectName:    earth-storage
 * @Package:        dggs
 * @Classname       AbstractDggsOperation
 * @Version         1.0.0
 * @Author          Guoliang PU
 * @Date            2021 1.29
 * @See
 *=============================================================================================
 * @Description     Define OGC DGGS AS operations that users can operate with DGGS.
 *                  The minimal set is in OGC Abstract Specification.
 *                  QuantisationOperations Set:
 *                      * toDataTiles()     --  将输入数据转换为数据瓦片（给cell赋值）
 *                      * toDataCells()     --  将输入数据转换为cell集合（给cell赋值）
 *                      * toCoordinates()   --  将输入数据转换为适当精度层级的cellID集合（presence）
 *                      * toTags()          --  将输入数据转换为地理标签（根据输入数据的空间占用范围生成地理标签）
 *                      * toGraphicCells()  --  将输入数据转换为图形cell
 *                      * toGraphicTiles()  --  将输入数据转换为图形瓦片
 *
 *                  AlgebraicOperations Set:
 *                      *  cellNavigation() -- 父子、邻居cell间的导航
 *                      *  spatialAnalysis()    --  利用九交模型计算cell之间的空间关系或空间查询对象间的空间关系
 *                  Interoperability Set:
 *                      *  query()      --  接受和解释外部的多种数据查询（SQL，WCPS，WCS，WPS等），并将其转换为内部DGGS的实现
 *                      *  broadcast()   --  按照外部的数据交付要求，将内部DGGS数据转换为交付格式（ASCII,GML,HDF,JSON,netCDF,XML等)输出
 *
 *                  generateDggs()  -- user can configure some DGGS parameters, and use this method to produce
 *
 */
public interface AbstractDggsOperation<TID,TCELL> {
    Log LOG = LogFactory.getLog(AbstractDggsOperation.class);



    @Target(ElementType.PARAMETER)
    @Retention(RetentionPolicy.RUNTIME)
    @interface Required {}

    @Target(ElementType.PARAMETER)
    @Retention(RetentionPolicy.RUNTIME)
    @interface Preferred {}

//    default void setup() {}
//    void construct();

    //AddressType:
    // GEO   -- longitude latitude
    // Q2DI  -- quad number and (i, j) coordinates on that quad
    // SEQNUM-- linear address (1 to size-of-DGG)
    // INTERLEAVE --   digit-interleaved form of Q2DI
    // PLANE --  (x, y) coordinates on unfolded ISEA plane
    // Q2DD -- quad number and (x, y) coordinates on that quad
    // PROJTRI -- triangle number and (x, y) coordinates within that triangle on the ISEA plane
    // VERTEX2DD -- vertex number, triangle number, and (x, y) coordinates on ISEA plane
    // AIGEN -- the polygon of the DGG cell in ARC/INFO Generate file format
//    void generateGrid();        // create all (or some portion) of the specified DGG
//    void gridStatistics();      // output a table of grid characteristics for the specified DGG
//    void transformPoints();     // perform address conversion for specified point set.
//    void binPointVals();        // bin a set of floating-point data values associated with geodetic coordinates into the cells of a DGG(s) specified
//    void binPointPresence();    // perform presence/absence binning into the DGG(s) specified
//


}
