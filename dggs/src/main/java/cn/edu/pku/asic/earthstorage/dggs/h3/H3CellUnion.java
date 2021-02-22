package cn.edu.pku.asic.earthstorage.dggs.h3;

import cn.edu.pku.asic.earthstorage.dggs.core.AbstractCellUnion;
import cn.edu.pku.asic.earthstorage.dggs.sphere.SphereCap;
import cn.edu.pku.asic.earthstorage.dggs.sphere.SphereLatLngRect;

/**
 * =============================================================================================
 *
 * @ProjectName dggs
 * @Package cn.edu.pku.asic.earthstorage.dggs.h3
 * @ClassName H3CellUnion
 * @Version 1.0.0
 * @Author Guoliang PU
 * @Date 2021/2/1 19:03
 * @See =============================================================================================
 * @Description
 */
public class H3CellUnion extends AbstractCellUnion<H3CellId,H3Cell> implements H3Region{
    /**
     * Return a bounding spherical cap.
     */
    @Override
    public SphereCap getCapBound() {
        return null;
    }

    /**
     * Return a bounding latitude-longitude rectangle.
     */
    @Override
    public SphereLatLngRect getRectBound() {
        return null;
    }

    /**
     * If this method returns true, the region completely contains the given cell.
     * Otherwise, either the region does not contain the cell or the containment
     * relationship could not be determined.
     *
     * @param cell
     */
    @Override
    public boolean contains(H3Cell cell) {
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
    public boolean mayIntersect(H3Cell cell) {
        return false;
    }
}
