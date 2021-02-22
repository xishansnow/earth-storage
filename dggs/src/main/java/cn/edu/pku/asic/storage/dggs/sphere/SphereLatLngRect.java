/*
 * Copyright 2005 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package cn.edu.pku.asic.storage.dggs.sphere;

import com.google.common.base.Preconditions;

/**
 * 球面经纬度矩形类
 * An S2LatLngRect represents a latitude-longitude rectangle. It is capable of
 * representing the empty and full rectangles as well as single points.
 *
 */

public strictfp class SphereLatLngRect implements SphereRegion {

  private final R1Interval lat;
  private final S1Interval lng;

  /**
   * Construct a rectangle from minimum and maximum latitudes and longitudes. If
   * lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude line.
   */
  public SphereLatLngRect(final SphereLatLng lo, final SphereLatLng hi) {
    lat = new R1Interval(lo.lat().radians(), hi.lat().radians());
    lng = new S1Interval(lo.lng().radians(), hi.lng().radians());
    // assert (isValid());
  }

  /** Construct a rectangle from latitude and longitude intervals. */
  public SphereLatLngRect(R1Interval lat, S1Interval lng) {
    this.lat = lat;
    this.lng = lng;
    // assert (isValid());
  }

  /** The canonical empty rectangle */
  public static SphereLatLngRect empty() {
    return new SphereLatLngRect(R1Interval.empty(), S1Interval.empty());
  }

  /** The canonical full rectangle. */
  public static SphereLatLngRect full() {
    return new SphereLatLngRect(fullLat(), fullLng());
  }

  /** The full allowable range of latitudes. */
  public static R1Interval fullLat() {
    return new R1Interval(-Sphere.M_PI_2, Sphere.M_PI_2);
  }

  /**
   * The full allowable range of longitudes.
   */
  public static S1Interval fullLng() {
    return S1Interval.full();
  }

  /**
   * Construct a rectangle from a center point (in lat-lng space) and size in
   * each dimension. If size.lng() is greater than 360 degrees it is clamped,
   * and latitudes greater than +/- 90 degrees are also clamped. So for example,
   * FromCenterSize((80,170),(20,20)) -> (lo=(60,150),hi=(90,-170)).
   */
  public static SphereLatLngRect fromCenterSize(SphereLatLng center, SphereLatLng size) {
    return fromPoint(center).expanded(size.mul(0.5));
  }

  /** Convenience method to construct a rectangle containing a single point. */
  public static SphereLatLngRect fromPoint(SphereLatLng p) {
    // assert (p.isValid());
    return new SphereLatLngRect(p, p);
  }

  /**
   * Convenience method to construct the minimal bounding rectangle containing
   * the two given points. This is equivalent to starting with an empty
   * rectangle and calling AddPoint() twice. Note that it is different than the
   * S2LatLngRect(lo, hi) constructor, where the first point is always used as
   * the lower-left corner of the resulting rectangle.
   */
  public static SphereLatLngRect fromPointPair(SphereLatLng p1, SphereLatLng p2) {
    // assert (p1.isValid() && p2.isValid());
    return new SphereLatLngRect(R1Interval.fromPointPair(p1.lat().radians(), p2
      .lat().radians()), S1Interval.fromPointPair(p1.lng().radians(), p2.lng()
      .radians()));
  }

  /**
   * Return a latitude-longitude rectangle that contains the edge from "a" to
   * "b". Both points must be unit-length. Note that the bounding rectangle of
   * an edge can be larger than the bounding rectangle of its endpoints.
   */
  public static SphereLatLngRect fromEdge(SpherePoint a, SpherePoint b) {
    // assert (S2.isUnitLength(a) && S2.isUnitLength(b));
    SphereLatLngRect r = fromPointPair(new SphereLatLng(a), new SphereLatLng(b));

    // Check whether the min/max latitude occurs in the edge interior.
    // We find the normal to the plane containing AB, and then a vector "dir" in
    // this plane that also passes through the equator. We use RobustCrossProd
    // to ensure that the edge normal is accurate even when the two points are
    // very close together.
    SpherePoint ab = Sphere.robustCrossProd(a, b);
    SpherePoint dir = SpherePoint.crossProd(ab, new SpherePoint(0, 0, 1));
    double da = dir.dotProd(a);
    double db = dir.dotProd(b);
    if (da * db >= 0) {
      // Minimum and maximum latitude are attained at the vertices.
      return r;
    }
    // Minimum/maximum latitude occurs in the edge interior. This affects the
    // latitude bounds but not the longitude bounds.
    double absLat = Math.acos(Math.abs(ab.z / ab.norm()));
    if (da < 0) {
      return new SphereLatLngRect(new R1Interval(r.lat().lo(), absLat), r.lng());
    } else {
      return new SphereLatLngRect(new R1Interval(-absLat, r.lat().hi()), r.lng());
    }
  }

  /**
   * Return true if the rectangle is valid, which essentially just means that
   * the latitude bounds do not exceed Pi/2 in absolute value and the longitude
   * bounds do not exceed Pi in absolute value.
   *
   */
  public boolean isValid() {
    // The lat/lng ranges must either be both empty or both non-empty.
    return (Math.abs(lat.lo()) <= Sphere.M_PI_2 && Math.abs(lat.hi()) <= Sphere.M_PI_2
      && lng.isValid() && lat.isEmpty() == lng.isEmpty());
  }

  // Accessor methods.
  public S1Angle latLo() {
    return S1Angle.radians(lat.lo());
  }

  public S1Angle latHi() {
    return S1Angle.radians(lat.hi());
  }

  public S1Angle lngLo() {
    return S1Angle.radians(lng.lo());
  }

  public S1Angle lngHi() {
    return S1Angle.radians(lng.hi());
  }

  public R1Interval lat() {
    return lat;
  }

  public S1Interval lng() {
    return lng;
  }

  public SphereLatLng lo() {
    return new SphereLatLng(latLo(), lngLo());
  }

  public SphereLatLng hi() {
    return new SphereLatLng(latHi(), lngHi());
  }

  /**
   * Return true if the rectangle is empty, i.e. it contains no points at all.
   */
  public boolean isEmpty() {
    return lat.isEmpty();
  }

  // Return true if the rectangle is full, i.e. it contains all points.
  public boolean isFull() {
    return lat.equals(fullLat()) && lng.isFull();
  }

  /**
   * Return true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses the 180
   * degree latitude line.
   */
  public boolean isInverted() {
    return lng.isInverted();
  }

  /** Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order. */
  public SphereLatLng getVertex(int k) {
    // Return the points in CCW order (SW, SE, NE, NW).
    switch (k) {
      case 0:
        return SphereLatLng.fromRadians(lat.lo(), lng.lo());
      case 1:
        return SphereLatLng.fromRadians(lat.lo(), lng.hi());
      case 2:
        return SphereLatLng.fromRadians(lat.hi(), lng.hi());
      case 3:
        return SphereLatLng.fromRadians(lat.hi(), lng.lo());
      default:
        throw new IllegalArgumentException("Invalid vertex index.");
    }
  }

  /**
   * Return the center of the rectangle in latitude-longitude space (in general
   * this is not the center of the region on the sphere).
   */
  public SphereLatLng getCenter() {
    return SphereLatLng.fromRadians(lat.getCenter(), lng.getCenter());
  }

  /**
   * Return the minimum distance (measured along the surface of the sphere)
   * from a given point to the rectangle (both its boundary and its interior).
   * The latLng must be valid.
   */
  public S1Angle getDistance(SphereLatLng p) {
    // The algorithm here is the same as in getDistance(S2LagLngRect), only
    // with simplified calculations.
    SphereLatLngRect a = this;

    Preconditions.checkState(!a.isEmpty());
    Preconditions.checkArgument(p.isValid());

    if (a.lng().contains(p.lng().radians())) {
      return S1Angle.radians(Math.max(0.0, Math.max(p.lat().radians() - a.lat().hi(),
                                                    a.lat().lo() - p.lat().radians())));
    }

    S1Interval interval = new S1Interval(a.lng().hi(), a.lng().complement().getCenter());
    double aLng = a.lng().lo();
    if (interval.contains(p.lng().radians())) {
      aLng = a.lng().hi();
    }

    SpherePoint lo = SphereLatLng.fromRadians(a.lat().lo(), aLng).toPoint();
    SpherePoint hi = SphereLatLng.fromRadians(a.lat().hi(), aLng).toPoint();
    SpherePoint loCrossHi =
        SphereLatLng.fromRadians(0, aLng - Sphere.M_PI_2).normalized().toPoint();
    return SphereEdgeUtil.getDistance(p.toPoint(), lo, hi, loCrossHi);
  }

  /**
   * Return the minimum distance (measured along the surface of the sphere) to
   * the given S2LatLngRect. Both S2LatLngRects must be non-empty.
   */
  public S1Angle getDistance(SphereLatLngRect other) {
    SphereLatLngRect a = this;
    SphereLatLngRect b = other;

    Preconditions.checkState(!a.isEmpty());
    Preconditions.checkArgument(!b.isEmpty());

    // First, handle the trivial cases where the longitude intervals overlap.
    if (a.lng().intersects(b.lng())) {
      if (a.lat().intersects(b.lat())) {
        return S1Angle.radians(0);  // Intersection between a and b.
      }

      // We found an overlap in the longitude interval, but not in the latitude
      // interval. This means the shortest path travels along some line of
      // longitude connecting the high-latitude of the lower rect with the
      // low-latitude of the higher rect.
      S1Angle lo, hi;
      if (a.lat().lo() > b.lat().hi()) {
        lo = b.latHi();
        hi = a.latLo();
      } else {
        lo = a.latHi();
        hi = b.latLo();
      }
      return S1Angle.radians(hi.radians() - lo.radians());
    }

    // The longitude intervals don't overlap. In this case, the closest points
    // occur somewhere on the pair of longitudinal edges which are nearest in
    // longitude-space.
    S1Angle aLng, bLng;
    S1Interval loHi = S1Interval.fromPointPair(a.lng().lo(), b.lng().hi());
    S1Interval hiLo = S1Interval.fromPointPair(a.lng().hi(), b.lng().lo());
    if (loHi.getLength() < hiLo.getLength()) {
      aLng = a.lngLo();
      bLng = b.lngHi();
    } else {
      aLng = a.lngHi();
      bLng = b.lngLo();
    }

    // The shortest distance between the two longitudinal segments will include
    // at least one segment endpoint. We could probably narrow this down further
    // to a single point-edge distance by comparing the relative latitudes of the
    // endpoints, but for the sake of clarity, we'll do all four point-edge
    // distance tests.
    SpherePoint aLo = new SphereLatLng(a.latLo(), aLng).toPoint();
    SpherePoint aHi = new SphereLatLng(a.latHi(), aLng).toPoint();
    SpherePoint aLoCrossHi =
        SphereLatLng.fromRadians(0, aLng.radians() - Sphere.M_PI_2).normalized().toPoint();
    SpherePoint bLo = new SphereLatLng(b.latLo(), bLng).toPoint();
    SpherePoint bHi = new SphereLatLng(b.latHi(), bLng).toPoint();
    SpherePoint bLoCrossHi =
        SphereLatLng.fromRadians(0, bLng.radians() - Sphere.M_PI_2).normalized().toPoint();

    return S1Angle.min(SphereEdgeUtil.getDistance(aLo, bLo, bHi, bLoCrossHi),
                       S1Angle.min(SphereEdgeUtil.getDistance(aHi, bLo, bHi, bLoCrossHi),
                                   S1Angle.min(SphereEdgeUtil.getDistance(bLo, aLo, aHi, aLoCrossHi),
                                               SphereEdgeUtil.getDistance(bHi, aLo, aHi, aLoCrossHi))));
  }

  /**
   * Return the width and height of this rectangle in latitude-longitude space.
   * Empty rectangles have a negative width and height.
   */
  public SphereLatLng getSize() {
    return SphereLatLng.fromRadians(lat.getLength(), lng.getLength());
  }

  /**
   * More efficient version of Contains() that accepts a S2LatLng rather than an
   * S2Point.
   */
  public boolean contains(SphereLatLng ll) {
    // assert (ll.isValid());
    return (lat.contains(ll.lat().radians()) && lng.contains(ll.lng()
      .radians()));

  }

  /**
   * Return true if and only if the given point is contained in the interior of
   * the region (i.e. the region excluding its boundary). The point 'p' does not
   * need to be normalized.
   */
  public boolean interiorContains(SpherePoint p) {
    return interiorContains(new SphereLatLng(p));
  }

  /**
   * More efficient version of InteriorContains() that accepts a S2LatLng rather
   * than an S2Point.
   */
  public boolean interiorContains(SphereLatLng ll) {
    // assert (ll.isValid());
    return (lat.interiorContains(ll.lat().radians()) && lng
      .interiorContains(ll.lng().radians()));
  }

  /**
   * Return true if and only if the rectangle contains the given other
   * rectangle.
   */
  public boolean contains(SphereLatLngRect other) {
    return lat.contains(other.lat) && lng.contains(other.lng);
  }

  /**
   * Return true if and only if the interior of this rectangle contains all
   * points of the given other rectangle (including its boundary).
   */
  public boolean interiorContains(SphereLatLngRect other) {
    return (lat.interiorContains(other.lat) && lng
      .interiorContains(other.lng));
  }

  /** Return true if this rectangle and the given other rectangle have any
  points in common. */
  public boolean intersects(SphereLatLngRect other) {
    return lat.intersects(other.lat) && lng.intersects(other.lng);
  }


  /**
   * Return true if and only if the interior of this rectangle intersects any
   * point (including the boundary) of the given other rectangle.
   */
  public boolean interiorIntersects(SphereLatLngRect other) {
    return (lat.interiorIntersects(other.lat) && lng
      .interiorIntersects(other.lng));
  }

  public SphereLatLngRect addPoint(SpherePoint p) {
    return addPoint(new SphereLatLng(p));
  }

  // Increase the size of the bounding rectangle to include the given point.
  // The rectangle is expanded by the minimum amount possible.
  public SphereLatLngRect addPoint(SphereLatLng ll) {
    // assert (ll.isValid());
    R1Interval newLat = lat.addPoint(ll.lat().radians());
    S1Interval newLng = lng.addPoint(ll.lng().radians());
    return new SphereLatLngRect(newLat, newLng);
  }

  /**
   * Return a rectangle that contains all points whose latitude distance from
   * this rectangle is at most margin.lat(), and whose longitude distance from
   * this rectangle is at most margin.lng(). In particular, latitudes are
   * clamped while longitudes are wrapped. Note that any expansion of an empty
   * interval remains empty, and both components of the given margin must be
   * non-negative.
   *
   * NOTE: If you are trying to grow a rectangle by a certain *distance* on the
   * sphere (e.g. 5km), use the ConvolveWithCap() method instead.
   */
  public SphereLatLngRect expanded(SphereLatLng margin) {
    // assert (margin.lat().radians() >= 0 && margin.lng().radians() >= 0);
    if (isEmpty()) {
      return this;
    }
    return new SphereLatLngRect(lat.expanded(margin.lat().radians()).intersection(
      fullLat()), lng.expanded(margin.lng().radians()));
  }

  /**
   * Return the smallest rectangle containing the union of this rectangle and
   * the given rectangle.
   */
  public SphereLatLngRect union(SphereLatLngRect other) {
    return new SphereLatLngRect(lat.union(other.lat), lng.union(other.lng));
  }

  /**
   * Return the smallest rectangle containing the intersection of this rectangle
   * and the given rectangle. Note that the region of intersection may consist
   * of two disjoint rectangles, in which case a single rectangle spanning both
   * of them is returned.
   */
  public SphereLatLngRect intersection(SphereLatLngRect other) {
    R1Interval intersectLat = lat.intersection(other.lat);
    S1Interval intersectLng = lng.intersection(other.lng);
    if (intersectLat.isEmpty() || intersectLng.isEmpty()) {
      // The lat/lng ranges must either be both empty or both non-empty.
      return empty();
    }
    return new SphereLatLngRect(intersectLat, intersectLng);
  }

  /**
   * Return a rectangle that contains the convolution of this rectangle with a
   * cap of the given angle. This expands the rectangle by a fixed distance (as
   * opposed to growing the rectangle in latitude-longitude space). The returned
   * rectangle includes all points whose minimum distance to the original
   * rectangle is at most the given angle.
   */
  public SphereLatLngRect convolveWithCap(S1Angle angle) {
    // The most straightforward approach is to build a cap centered on each
    // vertex and take the union of all the bounding rectangles (including the
    // original rectangle; this is necessary for very large rectangles).

    // Optimization: convert the angle to a height exactly once.
    SphereCap cap = SphereCap.fromAxisAngle(new SpherePoint(1, 0, 0), angle);

    SphereLatLngRect r = this;
    for (int k = 0; k < 4; ++k) {
      SphereCap vertexCap = SphereCap.fromAxisHeight(getVertex(k).toPoint(), cap
        .height());
      r = r.union(vertexCap.getRectBound());
    }
    return r;
  }

  /** Return the surface area of this rectangle on the unit sphere. */
  public double area() {
    if (isEmpty()) {
      return 0;
    }

    // This is the size difference of the two spherical caps, multiplied by
    // the longitude ratio.
    return lng().getLength() * Math.abs(Math.sin(latHi().radians()) - Math.sin(latLo().radians()));
  }

  /** Return true if two rectangles contains the same set of points. */
  @Override
  public boolean equals(Object that) {
    if (!(that instanceof SphereLatLngRect)) {
      return false;
    }
    SphereLatLngRect otherRect = (SphereLatLngRect) that;
    return lat().equals(otherRect.lat()) && lng().equals(otherRect.lng());
  }

  /**
   * Return true if the latitude and longitude intervals of the two rectangles
   * are the same up to the given tolerance (see r1interval.h and s1interval.h
   * for details).
   */
  public boolean approxEquals(SphereLatLngRect other, double maxError) {
    return (lat.approxEquals(other.lat, maxError) && lng.approxEquals(
      other.lng, maxError));
  }

  public boolean approxEquals(SphereLatLngRect other) {
    return approxEquals(other, 1e-15);
  }

  @Override
  public int hashCode() {
    int value = 17;
    value = 37 * value + lat.hashCode();
    return (37 * value + lng.hashCode());
  }

  // //////////////////////////////////////////////////////////////////////
  // S2Region interface (see {@code S2Region} for details):

  @Override
  public SphereRegion clone() {
    return new SphereLatLngRect(this.lo(), this.hi());
  }

  @Override
  public SphereCap getCapBound() {
    // We consider two possible bounding caps, one whose axis passes
    // through the center of the lat-long rectangle and one whose axis
    // is the north or south pole. We return the smaller of the two caps.

    if (isEmpty()) {
      return SphereCap.empty();
    }

    double poleZ, poleAngle;
    if (lat.lo() + lat.hi() < 0) {
      // South pole axis yields smaller cap.
      poleZ = -1;
      poleAngle = Sphere.M_PI_2 + lat.hi();
    } else {
      poleZ = 1;
      poleAngle = Sphere.M_PI_2 - lat.lo();
    }
    SphereCap poleCap = SphereCap.fromAxisAngle(new SpherePoint(0, 0, poleZ), S1Angle
      .radians(poleAngle));

    // For bounding rectangles that span 180 degrees or less in longitude, the
    // maximum cap size is achieved at one of the rectangle vertices. For
    // rectangles that are larger than 180 degrees, we punt and always return a
    // bounding cap centered at one of the two poles.
    double lngSpan = lng.hi() - lng.lo();
    if (Math.IEEEremainder(lngSpan, 2 * Sphere.M_PI) >= 0) {
      if (lngSpan < 2 * Sphere.M_PI) {
        SphereCap midCap = SphereCap.fromAxisAngle(getCenter().toPoint(), S1Angle
          .radians(0));
        for (int k = 0; k < 4; ++k) {
          midCap = midCap.addPoint(getVertex(k).toPoint());
        }
        if (midCap.height() < poleCap.height()) {
          return midCap;
        }
      }
    }
    return poleCap;
  }

  @Override
  public SphereLatLngRect getRectBound() {
    return this;
  }

  /** The point 'p' does not need to be normalized. */
  public boolean contains(SpherePoint p) {
    return contains(new SphereLatLng(p));
  }

  /**
   * Return true if the edge AB intersects the given edge of constant longitude.
   */
  private static boolean intersectsLngEdge(SpherePoint a, SpherePoint b,
                                           R1Interval lat, double lng) {
    // Return true if the segment AB intersects the given edge of constant
    // longitude. The nice thing about edges of constant longitude is that
    // they are straight lines on the sphere (geodesics).

    return Sphere.simpleCrossing(a, b, SphereLatLng.fromRadians(lat.lo(), lng)
      .toPoint(), SphereLatLng.fromRadians(lat.hi(), lng).toPoint());
  }

  /**
   * Return true if the edge AB intersects the given edge of constant latitude.
   */
  private static boolean intersectsLatEdge(SpherePoint a, SpherePoint b, double lat,
                                           S1Interval lng) {
    // Return true if the segment AB intersects the given edge of constant
    // latitude. Unfortunately, lines of constant latitude are curves on
    // the sphere. They can intersect a straight edge in 0, 1, or 2 points.
    // assert (S2.isUnitLength(a) && S2.isUnitLength(b));

    // First, compute the normal to the plane AB that points vaguely north.
    SpherePoint z = SpherePoint.normalize(Sphere.robustCrossProd(a, b));
    if (z.z < 0) {
      z = SpherePoint.neg(z);
    }

    // Extend this to an orthonormal frame (x,y,z) where x is the direction
    // where the great circle through AB achieves its maximium latitude.
    SpherePoint y = SpherePoint.normalize(Sphere.robustCrossProd(z, new SpherePoint(0, 0, 1)));
    SpherePoint x = SpherePoint.crossProd(y, z);
    // assert (S2.isUnitLength(x) && x.z >= 0);

    // Compute the angle "theta" from the x-axis (in the x-y plane defined
    // above) where the great circle intersects the given line of latitude.
    double sinLat = Math.sin(lat);
    if (Math.abs(sinLat) >= x.z) {
      return false; // The great circle does not reach the given latitude.
    }
    // assert (x.z > 0);
    double cosTheta = sinLat / x.z;
    double sinTheta = Math.sqrt(1 - cosTheta * cosTheta);
    double theta = Math.atan2(sinTheta, cosTheta);

    // The candidate intersection points are located +/- theta in the x-y
    // plane. For an intersection to be valid, we need to check that the
    // intersection point is contained in the interior of the edge AB and
    // also that it is contained within the given longitude interval "lng".

    // Compute the range of theta values spanned by the edge AB.
    S1Interval abTheta = S1Interval.fromPointPair(Math.atan2(
      a.dotProd(y), a.dotProd(x)), Math.atan2(b.dotProd(y), b.dotProd(x)));

    if (abTheta.contains(theta)) {
      // Check if the intersection point is also in the given "lng" interval.
      SpherePoint isect = SpherePoint.add(SpherePoint.mul(x, cosTheta), SpherePoint.mul(y,
        sinTheta));
      if (lng.contains(Math.atan2(isect.y, isect.x))) {
        return true;
      }
    }
    if (abTheta.contains(-theta)) {
      // Check if the intersection point is also in the given "lng" interval.
      SpherePoint intersection = SpherePoint.sub(SpherePoint.mul(x, cosTheta), SpherePoint.mul(y, sinTheta));
      if (lng.contains(Math.atan2(intersection.y, intersection.x))) {
        return true;
      }
    }
    return false;

  }

  @Override
  public String toString() {
    return "[Lo=" + lo() + ", Hi=" + hi() + "]";
  }
}
