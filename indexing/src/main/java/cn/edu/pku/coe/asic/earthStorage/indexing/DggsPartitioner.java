package cn.edu.pku.coe.asic.earthStorage.indexing;


import edu.ucr.cs.bdlab.beast.cg.SpatialPartitioner;
import edu.ucr.cs.bdlab.beast.common.BeastOptions;
import edu.ucr.cs.bdlab.beast.geolite.EnvelopeND;
import edu.ucr.cs.bdlab.beast.geolite.EnvelopeNDLite;
import edu.ucr.cs.bdlab.beast.synopses.AbstractHistogram;
import edu.ucr.cs.bdlab.beast.synopses.Summary;
import edu.ucr.cs.bdlab.beast.util.IntArray;


import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;
import java.io.PrintStream;
import java.util.Iterator;

/**
 *=============================================================================================
 * @ProjectName:    EarthStorage
 * @Package:        cn.edu.pku.coe.asic.earthStorage.indexing
 * @Classname       DggsPartitioner
 * @Version         v1.0.0
 * @Author          Guoliang PU
 * @Date            2021/1/1 17:08
 *=============================================================================================
 * @Description     TODO
 * @See
 */

public abstract class DggsPartitioner implements  SpatialPartitioner {
   /**-----------------------------------------------------------------------------
   * setup()        TODO
   *-------------------------------------------------------------------------------
   * @param		{conf}
   * @param		{disjoint}
   * @Return       None
   * @Exception    None
   * @Author       Guoliang PU {pgl@pku.edu.cn}
   * @Date         2021/1/1 17:46
   */
   @Override
   public void setup(BeastOptions conf, boolean disjoint) {
    //conf.

    }

    /**-----------------------------------------------------------------------------
    * construct()        TODO
    *-------------------------------------------------------------------------------
    * @param		{summary}
    * @param		{doubles}
    * @param		{abstractHistogram}
    * @param		{i}
    * @Return       None
    * @Exception    None
    * @Author       Guoliang PU {pgl@pku.edu.cn}
    * @Date         2021/1/1 17:46
    */
    @Override
    public void construct(Summary summary, double[][] doubles, AbstractHistogram abstractHistogram, int i) {


    }
}
