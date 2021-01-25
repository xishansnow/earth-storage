package cn.edu.pku.asic.storage.partitioner;

import cn.edu.pku.asic.storage.common.cli.BeastOptions;
import cn.edu.pku.asic.storage.common.geolite.EnvelopeND;
import cn.edu.pku.asic.storage.common.geolite.EnvelopeNDLite;
import cn.edu.pku.asic.storage.common.synopses.AbstractHistogram;
import cn.edu.pku.asic.storage.common.synopses.Summary;
import cn.edu.pku.asic.storage.common.utils.IntArray;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;
import java.io.PrintStream;
import java.util.Iterator;

/**
 * @author Guoliang PU
 * @description Partitioner for S2 grid system
 * @see
 * @see
 */
public class S2GeometryPartitioner extends DggsPartitioner{

    @Override
    public void construct(Summary summary, double[][] doubles, AbstractHistogram abstractHistogram, int i) {
        super.construct(summary, doubles, abstractHistogram, i);
    }

    @Override
    public void setup(BeastOptions conf, boolean disjoint) {
        super.setup(conf, disjoint);
    }

    @Override
    public void overlapPartitions(EnvelopeNDLite envelopeNDLite, IntArray intArray) {

    }

    @Override
    public void overlapPartitions(EnvelopeND mbr, IntArray matchedPartitions) {

    }

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

    @Override
    public void toWKT(PrintStream out) {

    }

    @Override
    public int overlapPartition(EnvelopeNDLite envelopeNDLite) {
        return 0;
    }

    @Override
    public void getPartitionMBR(int i, EnvelopeNDLite envelopeNDLite) {

    }

    @Override
    public int getPartitionCount() {
        return 0;
    }

    @Override
    public boolean isDisjoint() {
        return false;
    }

    @Override
    public int getCoordinateDimension() {
        return 0;
    }

    @Override
    public EnvelopeNDLite getEnvelope() {
        return null;
    }

    /**
     * The object implements the writeExternal method to save its contents
     * by calling the methods of DataOutput for its primitive values or
     * calling the writeObject method of ObjectOutput for objects, strings,
     * and arrays.
     *
     * @param out the stream to write the object to
     * @throws IOException Includes any I/O exceptions that may occur
     * @serialData Overriding methods should use this tag to describe
     * the data layout of this Externalizable object.
     * List the sequence of element types and, if possible,
     * relate the element to a public/protected field and/or
     * method of this Externalizable class.
     */
    @Override
    public void writeExternal(ObjectOutput out) throws IOException {

    }

    /**
     * The object implements the readExternal method to restore its
     * contents by calling the methods of DataInput for primitive
     * types and readObject for objects, strings and arrays.  The
     * readExternal method must read the values in the same sequence
     * and with the same types as were written by writeExternal.
     *
     * @param in the stream to read data from in order to restore the object
     * @throws IOException            if I/O errors occur
     * @throws ClassNotFoundException If the class for an object being
     *                                restored cannot be found.
     */
    @Override
    public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {

    }
}
