package cn.edu.pku.asic.storage.dggs.s2geometry;

import cn.edu.pku.asic.storage.dggs.sphere.SpherePolygon;

/**
 * =============================================================================================
 *
 * @ProjectName dggs
 * @Package cn.edu.pku.asic.storage.dggs.s2geometry
 * @ClassName S2Polygon
 * @Version 1.0.0
 * @Author Guoliang PU
 * @Date 2021/2/1 20:06
 * @See =============================================================================================
 * @Description
 */
public class S2Polygon extends SpherePolygon implements S2RegionInterface {
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
