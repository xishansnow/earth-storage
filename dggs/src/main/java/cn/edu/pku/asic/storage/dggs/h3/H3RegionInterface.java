package cn.edu.pku.asic.storage.dggs.h3;

import cn.edu.pku.asic.storage.dggs.s2geometry.S2Cell;

public interface H3RegionInterface {
    /**
     * If this method returns true, the region completely contains the given cell.
     * Otherwise, either the region does not contain the cell or the containment
     * relationship could not be determined.
     */
    public abstract boolean contains(H3Cell cell);

    /**
     * If this method returns false, the region does not intersect the given cell.
     * Otherwise, either region intersects the cell, or the intersection
     * relationship could not be determined.
     */
    public abstract boolean mayIntersect(H3Cell cell);
}
