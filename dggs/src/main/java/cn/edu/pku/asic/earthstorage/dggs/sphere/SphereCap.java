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
package cn.edu.pku.asic.earthstorage.dggs.sphere;


/**
 * 球冠类
 * This class represents a spherical cap, i.e. a portion of a sphere cut off by
 * a plane. The cap is defined by its axis and height. This representation has
 * good numerical accuracy for very small caps (unlike the (axis,
 * min-distance-from-origin) representation), and is also efficient for
 * containment tests (unlike the (axis, angle) representation).
 *
 * Here are some useful relationships between the cap height (h), the cap
 * opening angle (theta), the maximum chord length from the cap's center (d),
 * and the radius of cap's base (a). All formulas assume a unit radius.
 *
 * h = 1 - cos(theta) = 2 sin^2(theta/2) d^2 = 2 h = a^2 + h^2

 *
 */
public final strictfp class SphereCap implements SphereRegion {

  /**
   * Multiply a positive number by this constant to ensure that the result of a
   * floating point operation is at least as large as the true
   * infinite-precision result.
   */
  private static final double ROUND_UP = 1.0 + 1.0 / (1L << 52);

  private final SpherePoint axis;
  private final double height;

  // Caps may be constructed from either an axis and a height, or an axis and
  // an angle. To avoid ambiguity, there are no public constructors
  private SphereCap() {
    axis = new SpherePoint();
    height = 0;
  }

  private SphereCap(SpherePoint axis, double height) {
    this.axis = axis;
    this.height = height;
    // assert (isValid());
  }

  /**
   * Create a cap given its axis and the cap height, i.e. the maximum projected
   * distance along the cap axis from the cap center. 'axis' should be a
   * unit-length vector.
   */
  public static SphereCap fromAxisHeight(SpherePoint axis, double height) {
    // assert (S2.isUnitLength(axis));
    return new SphereCap(axis, height);
  }

  /**
   * Create a cap given its axis and the cap opening angle, i.e. maximum angle
   * between the axis and a point on the cap. 'axis' should be a unit-length
   * vector, and 'angle' should be between 0 and 180 degrees.
   */
  public static SphereCap fromAxisAngle(SpherePoint axis, S1Angle angle) {
    // The height of the cap can be computed as 1-cos(angle), but this isn't
    // very accurate for angles close to zero (where cos(angle) is almost 1).
    // Computing it as 2*(sin(angle/2)**2) gives much better precision.

    // assert (S2.isUnitLength(axis));
    double d = Math.sin(0.5 * angle.radians());
    return new SphereCap(axis, 2 * d * d);

  }

  /**
   * Create a cap given its axis and its area in steradians. 'axis' should be a
   * unit-length vector, and 'area' should be between 0 and 4 * M_PI.
   */
  public static SphereCap fromAxisArea(SpherePoint axis, double area) {
    // assert (S2.isUnitLength(axis));
    return new SphereCap(axis, area / (2 * Sphere.M_PI));
  }

  /** Return an empty cap, i.e. a cap that contains no points. */
  public static SphereCap empty() {
    return new SphereCap(new SpherePoint(1, 0, 0), -1);
  }

  /** Return a full cap, i.e. a cap that contains all points. */
  public static SphereCap full() {
    return new SphereCap(new SpherePoint(1, 0, 0), 2);
  }


  // Accessor methods.
  public SpherePoint axis() {
    return axis;
  }

  public double height() {
    return height;
  }

  public double area() {
    return 2 * Sphere.M_PI * Math.max(0.0, height);
  }

  /**
   * Return the cap opening angle in radians, or a negative number for empty
   * caps.
   */
  public S1Angle angle() {
    // This could also be computed as acos(1 - height_), but the following
    // formula is much more accurate when the cap height is small. It
    // follows from the relationship h = 1 - cos(theta) = 2 sin^2(theta/2).
    if (isEmpty()) {
      return S1Angle.radians(-1);
    }
    return S1Angle.radians(2 * Math.asin(Math.sqrt(0.5 * height)));
  }

  /**
   * We allow negative heights (to represent empty caps) but not heights greater
   * than 2.
   */
  public boolean isValid() {
    return Sphere.isUnitLength(axis) && height <= 2;
  }

  /** Return true if the cap is empty, i.e. it contains no points. */
  public boolean isEmpty() {
    return height < 0;
  }

  /** Return true if the cap is full, i.e. it contains all points. */
  public boolean isFull() {
    return height >= 2;
  }

  /**
   * Return the complement of the interior of the cap. A cap and its complement
   * have the same boundary but do not share any interior points. The complement
   * operator is not a bijection, since the complement of a singleton cap
   * (containing a single point) is the same as the complement of an empty cap.
   */
  public SphereCap complement() {
    // The complement of a full cap is an empty cap, not a singleton.
    // Also make sure that the complement of an empty cap has height 2.
    double cHeight = isFull() ? -1 : 2 - Math.max(height, 0.0);
    return SphereCap.fromAxisHeight(SpherePoint.neg(axis), cHeight);
  }

  /**
   * Return true if and only if this cap contains the given other cap (in a set
   * containment sense, e.g. every cap contains the empty cap).
   */
  public boolean contains(SphereCap other) {
    if (isFull() || other.isEmpty()) {
      return true;
    }
    return angle().radians() >= axis.angle(other.axis)
      + other.angle().radians();
  }

  /**
   * Return true if and only if the interior of this cap intersects the given
   * other cap. (This relationship is not symmetric, since only the interior of
   * this cap is used.)
   */
  public boolean interiorIntersects(SphereCap other) {
    // Interior(X) intersects Y if and only if Complement(Interior(X))
    // does not contain Y.
    return !complement().contains(other);
  }

  /**
   * Return true if and only if the given point is contained in the interior of
   * the region (i.e. the region excluding its boundary). 'p' should be a
   * unit-length vector.
   */
  public boolean interiorContains(SpherePoint p) {
    // assert (S2.isUnitLength(p));
    return isFull() || SpherePoint.sub(axis, p).norm2() < 2 * height;
  }

  /**
   * Increase the cap height if necessary to include the given point. If the cap
   * is empty the axis is set to the given point, but otherwise it is left
   * unchanged. 'p' should be a unit-length vector.
   */
  public SphereCap addPoint(SpherePoint p) {
    // Compute the squared chord length, then convert it into a height.
    // assert (S2.isUnitLength(p));
    if (isEmpty()) {
      return new SphereCap(p, 0);
    } else {
      // To make sure that the resulting cap actually includes this point,
      // we need to round up the distance calculation. That is, after
      // calling cap.AddPoint(p), cap.Contains(p) should be true.
      double dist2 = SpherePoint.sub(axis, p).norm2();
      double newHeight = Math.max(height, ROUND_UP * 0.5 * dist2);
      return new SphereCap(axis, newHeight);
    }
  }

  // Increase the cap height if necessary to include "other". If the current
  // cap is empty it is set to the given other cap.
  public SphereCap addCap(SphereCap other) {
    if (isEmpty()) {
      return new SphereCap(other.axis, other.height);
    } else {
      // See comments for FromAxisAngle() and AddPoint(). This could be
      // optimized by doing the calculation in terms of cap heights rather
      // than cap opening angles.
      double angle = axis.angle(other.axis) + other.angle().radians();
      if (angle >= Sphere.M_PI) {
        return new SphereCap(axis, 2); //Full cap
      } else {
        double d = Math.sin(0.5 * angle);
        double newHeight = Math.max(height, ROUND_UP * 2 * d * d);
        return new SphereCap(axis, newHeight);
      }
    }
  }

  // //////////////////////////////////////////////////////////////////////
  // S2Region interface (see {@code S2Region} for details):
  @Override
  public SphereCap getCapBound() {
    return this;
  }

  @Override
  public SphereLatLngRect getRectBound() {
    if (isEmpty()) {
      return SphereLatLngRect.empty();
    }

    // Convert the axis to a (lat,lng) pair, and compute the cap angle.
    SphereLatLng axisLatLng = new SphereLatLng(axis);
    double capAngle = angle().radians();

    boolean allLongitudes = false;
    double[] lat = new double[2], lng = new double[2];
    lng[0] = -Sphere.M_PI;
    lng[1] = Sphere.M_PI;

    // Check whether cap includes the south pole.
    lat[0] = axisLatLng.lat().radians() - capAngle;
    if (lat[0] <= -Sphere.M_PI_2) {
      lat[0] = -Sphere.M_PI_2;
      allLongitudes = true;
    }
    // Check whether cap includes the north pole.
    lat[1] = axisLatLng.lat().radians() + capAngle;
    if (lat[1] >= Sphere.M_PI_2) {
      lat[1] = Sphere.M_PI_2;
      allLongitudes = true;
    }
    if (!allLongitudes) {
      // Compute the range of longitudes covered by the cap. We use the law
      // of sines for spherical triangles. Consider the triangle ABC where
      // A is the north pole, B is the center of the cap, and C is the point
      // of tangency between the cap boundary and a line of longitude. Then
      // C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
      // we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
      // Here "a" is the cap angle, and "c" is the colatitude (90 degrees
      // minus the latitude). This formula also works for negative latitudes.
      //
      // The formula for sin(a) follows from the relationship h = 1 - cos(a).

      double sinA = Math.sqrt(height * (2 - height));
      double sinC = Math.cos(axisLatLng.lat().radians());
      if (sinA <= sinC) {
        double angleA = Math.asin(sinA / sinC);
        lng[0] = Math.IEEEremainder(axisLatLng.lng().radians() - angleA,
          2 * Sphere.M_PI);
        lng[1] = Math.IEEEremainder(axisLatLng.lng().radians() + angleA,
          2 * Sphere.M_PI);
      }
    }
    return new SphereLatLngRect(new R1Interval(lat[0], lat[1]), new S1Interval(lng[0], lng[1]));
  }


  public boolean contains(SpherePoint p) {
    // The point 'p' should be a unit-length vector.
    // assert (S2.isUnitLength(p));
    return SpherePoint.sub(axis, p).norm2() <= 2 * height;

  }


  /** Return true if two caps are identical. */
  @Override
  public boolean equals(Object that) {

    if (!(that instanceof SphereCap)) {
      return false;
    }

    SphereCap other = (SphereCap) that;
    return (axis.equals(other.axis) && height == other.height)
        || (isEmpty() && other.isEmpty()) || (isFull() && other.isFull());

  }

  @Override
  public int hashCode() {
    if (isFull()) {
      return 17;
    } else if (isEmpty()) {
      return 37;
    }
    int result = 17;
    result = 37 * result + axis.hashCode();
    long heightBits = Double.doubleToLongBits(height);
    result = 37 * result + (int) ((heightBits >>> 32) ^ heightBits);
    return result;
  }

  // /////////////////////////////////////////////////////////////////////
  // The following static methods are convenience functions for assertions
  // and testing purposes only.

  /**
   * Return true if the cap axis and height differ by at most "max_error" from
   * the given cap "other".
   */
  boolean approxEquals(SphereCap other, double maxError) {
    return (axis.aequal(other.axis, maxError) && Math.abs(height - other.height) <= maxError)
      || (isEmpty() && other.height <= maxError)
      || (other.isEmpty() && height <= maxError)
      || (isFull() && other.height >= 2 - maxError)
      || (other.isFull() && height >= 2 - maxError);
  }

  boolean approxEquals(SphereCap other) {
    return approxEquals(other, 1e-14);
  }

  @Override
  public String toString() {
    return "[Point = " + axis.toString() + " Height = " + height + "]";
  }
}
