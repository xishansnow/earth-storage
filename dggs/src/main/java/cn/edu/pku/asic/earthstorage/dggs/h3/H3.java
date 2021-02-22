package cn.edu.pku.asic.earthstorage.dggs.h3;

import cn.edu.pku.asic.earthstorage.dggs.core.AbstractDggs;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * =============================================================================================
 *
 * @ProjectName dggs
 * @Package cn.edu.pku.asic.earthstorage.dggs.h3
 * @ClassName H3
 * @Version 1.0.0
 * @Author Guoliang PU
 * @Date 2021/2/1 19:09
 * @See =============================================================================================
 * @Description
 */
public class H3 extends AbstractDggs<H3CellId,H3Cell,H3CellUnion> {
    private static H3 h3Instance;

    private H3() {
    }

    public static H3 getInstance() {
        if (h3Instance == null) {
            h3Instance = new H3();
        }
        return h3Instance;
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
