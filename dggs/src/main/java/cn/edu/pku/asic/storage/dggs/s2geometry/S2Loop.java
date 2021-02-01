package cn.edu.pku.asic.storage.dggs.s2geometry;

import cn.edu.pku.asic.storage.dggs.sphere.SphereLoop;
import cn.edu.pku.asic.storage.dggs.sphere.SpherePoint;

import java.util.List;

/**
 * =============================================================================================
 *
 * @ProjectName dggs
 * @Package cn.edu.pku.asic.storage.dggs.s2geometry
 * @ClassName S2Loop
 * @Version 1.0.0
 * @Author Guoliang PU
 * @Date 2021/2/1 20:10
 * @See =============================================================================================
 * @Description
 */
public class S2Loop extends SphereLoop implements S2Region {
    /**
     * Initialize a loop connecting the given vertices. The last vertex is
     * implicitly connected to the first. All points should be unit length. Loops
     * must have at least 3 vertices.
     *
     * @param vertices
     */
    public S2Loop(List<SpherePoint> vertices) {
        super(vertices);
    }

    /**
     * Copy constructor.
     *
     * @param src
     */
    public S2Loop(SphereLoop src) {
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
