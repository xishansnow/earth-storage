package cn.edu.pku.asic.earthstorage.dggs;

import com.google.common.primitives.UnsignedLong;


/**
 * @author Guoliang PU
 * @description A coding class for Z-Ordering curve.
 * @date 2021/1/1 10:33
 * @see
 */
public class ZOrderCodingUtil {
    /**
     * IJtoCellID:   TODO
     *
     * @param: I
     * @param: J
     * @param: level
     * @return： {@code java.lang.String}
     * @exception: None
     * @author： Guoliang PU {pgl@pku.edu.cn}
     * @date： 2021/1/1
     */
    public static String IJtoCellID(long I, long J, int level) {
        UnsignedLong uLong;
        long lLong = -1020;
        lLong = lLong << 3;

        return ("OK");
    }

    /**
     * cellIDtoIJ:   TODO
     *
     * @param: cellID
     * @param: I
     * @param: J
     * @param: level
     * @return： {@code boolean}
     * @exception: None
     * @author： Guoliang PU {pgl@pku.edu.cn}
     * @date： 2021/1/1
     */
    public static boolean cellIDtoIJ(long cellID, long I, long J, int level) {
        return (true);
    }

}
