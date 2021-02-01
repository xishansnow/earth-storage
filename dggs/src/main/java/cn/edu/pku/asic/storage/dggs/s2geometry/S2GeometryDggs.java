package cn.edu.pku.asic.storage.dggs.s2geometry;

import cn.edu.pku.asic.storage.dggs.core.AbstractDggs;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 *=============================================================================================
 * @ProjectName    EarthStorage 
 * @Package        cn.edu.pku.asic.storage.dggs
 * @ClassName      S2GeometryDggs
 * @Version        1.0.0
 * @Author         Guoliang PU
 * @Date           2021/1/29 13:42
 * @See
 *=============================================================================================
 * @Description    DGG system for Google's S2.            
 */
//@AbstractDggs.Metadata(
//        dggsType = "",
//        dggsTopology = "abs",
//        dggsProj ="",
//        description = "S2Geometry DGG system wrapper."
//
//
//)
public class S2GeometryDggs extends AbstractDggs {


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
