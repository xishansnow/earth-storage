package cn.edu.pku.asic.storage.partitioner;

import java.io.Externalizable;
import java.io.PrintStream;
import java.util.Iterator;
import cn.edu.pku.asic.storage.common.cg.SpatialPartitioner;

/**
 * =============================================================================================
 *
 * @ProjectName: EarthStorage
 * @Package: cn.edu.pku.coe.asic.earthStorage.partitioner
 * @Classname DggsPartitioner
 * @Version v1.0.0
 * @Author Guoliang PU
 * @Date 2021/1/1 17:08
 * =============================================================================================
 * @Description An abstract class for DGGS grid system. You can extended it by your own grid system.
 *
 * @See GeohashPartitioner, GeosotPartitioner, HexPartitioner, S2GeometryPartitioner
 */

public abstract class DggsPartitioner implements Externalizable {

    /**
     * The MBR of the input domain.
     */
    final protected EnvelopeNDLite inputMBR = new EnvelopeNDLite();

    protected GridSortable hGridSorter;


    /**
     * -----------------------------------------------------------------------------
     * setup()        TODO
     * -------------------------------------------------------------------------------
     *
     * @param        {conf}
     * @param        {disjoint}
     * @Return None
     * @Exception None
     * @Author Guoliang PU {pgl@pku.edu.cn}
     * @Date 2021/1/1 17:46
     */
    @Override
    public void setup(BeastOptions conf, boolean disjoint) {
        if (disjoint)
            //disjoint partitions not support.
            throw new RuntimeException("DGGS grid system does support disjoint partitions!");
    }


    /**
     * -----------------------------------------------------------------------------
     * construct()
     * Constructs the histogram using the given parameters. By default, only the summary parameter is provided
     * and the other two parameters {@code sample} and {@code histogram} are {@code null}. If any of these are needed,
     * the parameter should be annotated with {@link Required}. If any of them can be used to improve the results but
     * are not necessarily required, they can be annotated with {@link Preferred}.
     * -------------------------------------------------------------------------------
     *
     * @return void
     * @throws
     * @param        summary the summary of the input
     * @param        sample     a sample points from the input
     * @param        histogram    a histogram of the input
     * @param        numPartitions the desired number of partitions. This is treated as a loose hint and not a strict value
     * @author Guoliang PU {pgl@pku.edu.cn}
     * @date 2021/1/2 21:05
     */
    @Override
    public void construct(Summary summary, @Required double[][] sample, @Required AbstractHistogram histogram, int numPartitions) {
        //number of samples
        int numDimensions = sample.length;
        int sampleCount = sample[0].length;

    }

    /**
     * Returns the MBR of the underlying dataset on which the partitioner was created.
     * * @return the envelope of the entire partitioner
     */
    @Override
    public EnvelopeNDLite getEnvelope() {
        return null;
    }

    /**
     * Returns total number of partitions
     *
     * @return the total number of partitions
     */
    @Override
    public int getPartitionCount() {
        return 0;
    }

    @Override
    public void overlapPartitions(EnvelopeND mbr, IntArray matchedPartitions) {

    }

    /**
     * Returns the single partition that contains the center of the given envelope.
     * This method should return a value in the range [0, {@link #getPartitionCount()}[
     *
     * @param mbr the minimum bounding rectangle to find a partition
     * @return the ID of the best partition that the given MBR should be assigned to.
     */
    @Override
    public int overlapPartition(EnvelopeND mbr) {
        return 0;
    }

    @Override
    public EnvelopeNDLite getPartitionMBR(int partitionID) {
        return null;
    }

    @Override
    public Metadata getMetadata() {
        return null;
    }

    @Override
    public Iterator<EnvelopeNDLite> iterator() {
        return null;
    }

    /**
     * Print all the partition boundaries to a tab-separated file.
     *
     * @param out the output to write to
     */
    @Override
    public void toWKT(PrintStream out) {

    }

    /**
     * Find all the overlapping partitions in an envelope (MBR).
     *
     * @param mbr               the MBR to find overlapping partitions
     * @param matchedPartitions an array to use for matched partitions
     */
    @Override
    public void overlapPartitions(EnvelopeNDLite mbr, IntArray matchedPartitions) {

    }

    @Override
    public int overlapPartition(EnvelopeNDLite mbr) {
        return 0;
    }

    /**
     * Returns the minimum bounding rectangle (MBR) of the given partition.
     *
     * @param partitionID the ID of the partition to get its MBR.
     * @param mbr         output parameter for the MBR to avoid allocation of objects
     */
    @Override
    public void getPartitionMBR(int partitionID, EnvelopeNDLite mbr) {

    }

    /**
     * Whether this partitioner is configured to produce disjoint partitions
     *
     * @return {@code true} if this partitions will produce spatially disjoint partitions
     */
    @Override
    public boolean isDisjoint() {
        return false;
    }

    /**
     * Number of dimensions for this partitioner
     *
     * @return the dimensions of the partitioner
     */
    @Override
    public int getCoordinateDimension() {
        return 0;
    }

}
