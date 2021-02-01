package cn.edu.pku.asic.storage.dggs.s2geometry;

import cn.edu.pku.asic.storage.dggs.sphere.SpherePoint;
import cn.edu.pku.asic.storage.dggs.sphere.SpherePolyline;

import java.util.List;

/**
 * =============================================================================================
 *
 * @ProjectName dggs
 * @Package cn.edu.pku.asic.storage.dggs.s2geometry
 * @ClassName S2PolyLine
 * @Version 1.0.0
 * @Author Guoliang PU
 * @Date 2021/2/1 20:07
 * @See =============================================================================================
 * @Description
 */
public class S2PolyLine extends SpherePolyline implements S2Region {
    /**
     * Create a polyline that connects the given vertices. Empty polylines are
     * allowed. Adjacent vertices should not be identical or antipodal. All
     * vertices should be unit length.
     *
     * @param vertices
     */
    public S2PolyLine(List<SpherePoint> vertices) {
        super(vertices);
    }

    /**
     * Copy constructor.
     * <p>
     * TODO(dbeaumont): Now that S2Polyline is immutable, remove this.
     *
     * @param src
     */
    public S2PolyLine(SpherePolyline src) {
        super(src);
    }

    /**
     * If this method returns true, the region completely contains the given cell.
     * Otherwise, either the region does not contain the cell or the containment
     * relationship could not be determined.
     *
     * @param cell
     */
    @Override
    public boolean contains(S2Cell cell) {
        return false;
    }

    /**
     * If this method returns false, the region does not intersect the given cell.
     * Otherwise, either region intersects the cell, or the intersection
     * relationship could not be determined.
     *
     * @param cell
     */
    @Override
    public boolean mayIntersect(S2Cell cell) {
        return false;
    }
}
