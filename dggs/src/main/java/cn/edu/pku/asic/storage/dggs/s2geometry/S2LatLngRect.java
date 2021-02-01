package cn.edu.pku.asic.storage.dggs.s2geometry;

import cn.edu.pku.asic.storage.dggs.sphere.R1Interval;
import cn.edu.pku.asic.storage.dggs.sphere.S1Interval;
import cn.edu.pku.asic.storage.dggs.sphere.SphereLatLng;
import cn.edu.pku.asic.storage.dggs.sphere.SphereLatLngRect;

/**
 * =============================================================================================
 *
 * @ProjectName dggs
 * @Package cn.edu.pku.asic.storage.dggs.s2geometry
 * @ClassName S2LatLngRect
 * @Version 1.0.0
 * @Author Guoliang PU
 * @Date 2021/2/1 20:11
 * @See =============================================================================================
 * @Description
 */
public class S2LatLngRect extends SphereLatLngRect implements S2Region {
    /**
     * Construct a rectangle from minimum and maximum latitudes and longitudes. If
     * lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude line.
     *
     * @param lo
     * @param hi
     */
    public S2LatLngRect(SphereLatLng lo, SphereLatLng hi) {
        super(lo, hi);
    }

    /**
     * Construct a rectangle from latitude and longitude intervals.
     *
     * @param lat
     * @param lng
     */
    public S2LatLngRect(R1Interval lat, S1Interval lng) {
        super(lat, lng);
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
