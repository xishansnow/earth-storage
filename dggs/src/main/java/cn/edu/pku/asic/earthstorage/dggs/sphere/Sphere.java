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

/*  球面上的常数和球面几何计算函数 */

public strictfp class Sphere {

  // Declare some frequently used constants
  public static final double M_PI = Math.PI;
  public static final double M_1_PI = 1.0 / Math.PI;
  public static final double M_PI_2 = Math.PI / 2.0;
  public static final double M_PI_4 = Math.PI / 4.0;
  public static final double M_SQRT2 = Math.sqrt(2);
  public static final double M_E = Math.E;


  /**
   * Return a unique "origin" on the sphere for operations that need a fixed
   * reference point. It should *not* be a point that is commonly used in edge
   * tests in order to avoid triggering code to handle degenerate cases. (This
   * rules out the north and south poles.)
   */
  public static SpherePoint origin() {
    return new SpherePoint(0, 1, 0);
  }

  /**
   * Return true if the given point is approximately unit length (this is mainly
   * useful for assertions).
   */
  public static boolean isUnitLength(SpherePoint p) {
    return Math.abs(p.norm2() - 1) <= 1e-15;
  }

  /**
   * 判断球面上ab和cd是否相交于内部某点
   * Return true if edge AB crosses CD at a point that is interior to both
   * edges. Properties:
   *
   *  (1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d)
   *  (2) SimpleCrossing(c,d,a,b) == SimpleCrossing(a,b,c,d)
   */
  public static boolean simpleCrossing(SpherePoint a, SpherePoint b, SpherePoint c, SpherePoint d) {
    // We compute SimpleCCW() for triangles ACB, CBD, BDA, and DAC. All
    // of these triangles need to have the same orientation (CW or CCW)
    // for an intersection to exist. Note that this is slightly more
    // restrictive than the corresponding definition for planar edges,
    // since we need to exclude pairs of line segments that would
    // otherwise "intersect" by crossing two antipodal points.

    SpherePoint ab = SpherePoint.crossProd(a, b);
    SpherePoint cd = SpherePoint.crossProd(c, d);

    double acb = -ab.dotProd(c);
    double cbd = -cd.dotProd(b);
    double bda = ab.dotProd(d);
    double dac = cd.dotProd(a);

    return (acb * cbd > 0) && (cbd * bda > 0) && (bda * dac > 0);
  }

  /**
   * 更为鲁棒的外积计算函数
   * Return a vector "c" that is orthogonal to the given unit-length vectors "a"
   * and "b". This function is similar to a.CrossProd(b) except that it does a
   * better job of ensuring orthogonality when "a" is nearly parallel to "b",
   * and it returns a non-zero result even when a == b or a == -b.
   *
   *  It satisfies the following properties (RCP == RobustCrossProd):
   *
   *  (1) RCP(a,b) != 0 for all a, b (2) RCP(b,a) == -RCP(a,b) unless a == b or
   * a == -b (3) RCP(-a,b) == -RCP(a,b) unless a == b or a == -b (4) RCP(a,-b)
   * == -RCP(a,b) unless a == b or a == -b
   */
  public static SpherePoint robustCrossProd(SpherePoint a, SpherePoint b) {
    // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
    // approaches zero. This leads to situations where a.CrossProd(b) is not
    // very orthogonal to "a" and/or "b". We could fix this using Gram-Schmidt,
    // but we also want b.RobustCrossProd(a) == -b.RobustCrossProd(a).
    //
    // The easiest fix is to just compute the cross product of (b+a) and (b-a).
    // Given that "a" and "b" are unit-length, this has good orthogonality to
    // "a" and "b" even if they differ only in the lowest bit of one component.

    // assert (isUnitLength(a) && isUnitLength(b));
    SpherePoint x = SpherePoint.crossProd(SpherePoint.add(b, a), SpherePoint.sub(b, a));
    if (!x.equals(new SpherePoint(0, 0, 0))) {
      return x;
    }

    // The only result that makes sense mathematically is to return zero, but
    // we find it more convenient to return an arbitrary orthogonal vector.
    return ortho(a);
  }

  /**
   * 返回正交于
   * Return a unit-length vector that is orthogonal to "a". Satisfies Ortho(-a)
   * = -Ortho(a) for all a.
   */
  public static SpherePoint ortho(SpherePoint a) {
    // The current implementation in SpherePoint has the property we need,
    // i.e. Ortho(-a) = -Ortho(a) for all a.
    return a.ortho();
  }

  /**
   * 计算球面三角形的面积
   * Return the area of triangle ABC. The method used is about twice as
   * expensive as Girard's formula, but it is numerically stable for both large
   * and very small triangles. The points do not need to be normalized. The area
   * is always positive.
   *
   *  The triangle area is undefined if it contains two antipodal points, and
   * becomes numerically unstable as the length of any edge approaches 180
   * degrees.
   */
  static double area(SpherePoint a, SpherePoint b, SpherePoint c) {
    // This method is based on l'Huilier's theorem,
    //
    // tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
    //
    // where E is the spherical excess of the triangle (i.e. its area),
    // a, b, c, are the side lengths, and
    // s is the semiperimeter (a + b + c) / 2 .
    //
    // The only significant source of error using l'Huilier's method is the
    // cancellation error of the terms (s-a), (s-b), (s-c). This leads to a
    // *relative* error of about 1e-16 * s / min(s-a, s-b, s-c). This compares
    // to a relative error of about 1e-15 / E using Girard's formula, where E is
    // the true area of the triangle. Girard's formula can be even worse than
    // this for very small triangles, e.g. a triangle with a true area of 1e-30
    // might evaluate to 1e-5.
    //
    // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
    // dmin = min(s-a, s-b, s-c). This basically includes all triangles
    // except for extremely long and skinny ones.
    //
    // Since we don't know E, we would like a conservative upper bound on
    // the triangle area in terms of s and dmin. It's possible to show that
    // E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
    // Using this, it's easy to show that we should always use l'Huilier's
    // method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore,
    // if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
    // k3 is about 0.1. Since the best case error using Girard's formula
    // is about 1e-15, this means that we shouldn't even consider it unless
    // s >= 3e-4 or so.

    // We use volatile doubles to force the compiler to truncate all of these
    // quantities to 64 bits. Otherwise it may compute a value of dmin > 0
    // simply because it chose to spill one of the intermediate values to
    // memory but not one of the others.
    final double sa = b.angle(c);
    final double sb = c.angle(a);
    final double sc = a.angle(b);
    final double s = 0.5 * (sa + sb + sc);
    if (s >= 3e-4) {
      // Consider whether Girard's formula might be more accurate.
      double s2 = s * s;
      double dmin = s - Math.max(sa, Math.max(sb, sc));
      if (dmin < 1e-2 * s * s2 * s2) {
        // This triangle is skinny enough to consider Girard's formula.
        double area = girardArea(a, b, c);
        if (dmin < s * (0.1 * area)) {
          return area;
        }
      }
    }
    // Use l'Huilier's formula.
    return 4
        * Math.atan(
            Math.sqrt(
                Math.max(0.0,
                    Math.tan(0.5 * s) * Math.tan(0.5 * (s - sa)) * Math.tan(0.5 * (s - sb))
                        * Math.tan(0.5 * (s - sc)))));
  }

  /**
   * Return the area of the triangle computed using Girard's formula. This is
   * slightly faster than the Area() method above is not accurate for very small
   * triangles.
   */
  public static double girardArea(SpherePoint a, SpherePoint b, SpherePoint c) {
    // This is equivalent to the usual Girard's formula but is slightly
    // more accurate, faster to compute, and handles a == b == c without
    // a special case.

    SpherePoint ab = SpherePoint.crossProd(a, b);
    SpherePoint bc = SpherePoint.crossProd(b, c);
    SpherePoint ac = SpherePoint.crossProd(a, c);
    return Math.max(0.0, ab.angle(ac) - ab.angle(bc) + bc.angle(ac));
  }

  /**
   * Like Area(), but returns a positive value for counterclockwise triangles
   * and a negative value otherwise.
   */
  public static double signedArea(SpherePoint a, SpherePoint b, SpherePoint c) {
    return area(a, b, c) * robustCCW(a, b, c);
  }

  // About centroids:
  // ----------------
  //
  // There are several notions of the "centroid" of a triangle. First, there
  // // is the planar centroid, which is simply the centroid of the ordinary
  // (non-spherical) triangle defined by the three vertices. Second, there is
  // the surface centroid, which is defined as the intersection of the three
  // medians of the spherical triangle. It is possible to show that this
  // point is simply the planar centroid projected to the surface of the
  // sphere. Finally, there is the true centroid (mass centroid), which is
  // defined as the area integral over the spherical triangle of (x,y,z)
  // divided by the triangle area. This is the point that the triangle would
  // rotate around if it was spinning in empty space.
  //
  // The best centroid for most purposes is the true centroid. Unlike the
  // planar and surface centroids, the true centroid behaves linearly as
  // regions are added or subtracted. That is, if you split a triangle into
  // pieces and compute the average of their centroids (weighted by triangle
  // area), the result equals the centroid of the original triangle. This is
  // not true of the other centroids.
  //
  // Also note that the surface centroid may be nowhere near the intuitive
  // "center" of a spherical triangle. For example, consider the triangle
  // with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
  // The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
  // within a distance of 2*eps of the vertex B. Note that the median from A
  // (the segment connecting A to the midpoint of BC) passes through S, since
  // this is the shortest path connecting the two endpoints. On the other
  // hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
  // the surface is a much more reasonable interpretation of the "center" of
  // this triangle.

  /**
   * Return the centroid of the planar triangle ABC. This can be normalized to
   * unit length to obtain the "surface centroid" of the corresponding spherical
   * triangle, i.e. the intersection of the three medians. However, note that
   * for large spherical triangles the surface centroid may be nowhere near the
   * intuitive "center" (see example above).
   */
  public static SpherePoint planarCentroid(SpherePoint a, SpherePoint b, SpherePoint c) {
    return new SpherePoint((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0, (a.z + b.z + c.z) / 3.0);
  }

  /**
   * Returns the true centroid of the spherical triangle ABC multiplied by the
   * signed area of spherical triangle ABC. The reasons for multiplying by the
   * signed area are (1) this is the quantity that needs to be summed to compute
   * the centroid of a union or difference of triangles, and (2) it's actually
   * easier to calculate this way.
   */
  public static SpherePoint trueCentroid(SpherePoint a, SpherePoint b, SpherePoint c) {
    // I couldn't find any references for computing the true centroid of a
    // spherical triangle... I have a truly marvellous demonstration of this
    // formula which this margin is too narrow to contain :)

    // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));
    double sina = SpherePoint.crossProd(b, c).norm();
    double sinb = SpherePoint.crossProd(c, a).norm();
    double sinc = SpherePoint.crossProd(a, b).norm();
    double ra = (sina == 0) ? 1 : (Math.asin(sina) / sina);
    double rb = (sinb == 0) ? 1 : (Math.asin(sinb) / sinb);
    double rc = (sinc == 0) ? 1 : (Math.asin(sinc) / sinc);

    // Now compute a point M such that M.X = rX * det(ABC) / 2 for X in A,B,C.
    SpherePoint x = new SpherePoint(a.x, b.x, c.x);
    SpherePoint y = new SpherePoint(a.y, b.y, c.y);
    SpherePoint z = new SpherePoint(a.z, b.z, c.z);
    SpherePoint r = new SpherePoint(ra, rb, rc);
    return new SpherePoint(0.5 * SpherePoint.crossProd(y, z).dotProd(r),
        0.5 * SpherePoint.crossProd(z, x).dotProd(r), 0.5 * SpherePoint.crossProd(x, y).dotProd(r));
  }

  /**
   * Return true if the points A, B, C are strictly counterclockwise. Return
   * false if the points are clockwise or colinear (i.e. if they are all
   * contained on some great circle).
   *
   *  Due to numerical errors, situations may arise that are mathematically
   * impossible, e.g. ABC may be considered strictly CCW while BCA is not.
   * However, the implementation guarantees the following:
   *
   *  If SimpleCCW(a,b,c), then !SimpleCCW(c,b,a) for all a,b,c.
   *
   * In other words, ABC and CBA are guaranteed not to be both CCW
   */
  public static boolean simpleCCW(SpherePoint a, SpherePoint b, SpherePoint c) {
    // We compute the signed volume of the parallelepiped ABC. The usual
    // formula for this is (AxB).C, but we compute it here using (CxA).B
    // in order to ensure that ABC and CBA are not both CCW. This follows
    // from the following identities (which are true numerically, not just
    // mathematically):
    //
    // (1) x.CrossProd(y) == -(y.CrossProd(x))
    // (2) (-x).DotProd(y) == -(x.DotProd(y))

    return SpherePoint.crossProd(c, a).dotProd(b) > 0;
  }

  /**
   * WARNING! This requires arbitrary precision arithmetic to be truly robust.
   * This means that for nearly colinear AB and AC, this function may return the
   * wrong answer.
   *
   * <p>
   * Like SimpleCCW(), but returns +1 if the points are counterclockwise and -1
   * if the points are clockwise. It satisfies the following conditions:
   *
   *  (1) RobustCCW(a,b,c) == 0 if and only if a == b, b == c, or c == a (2)
   * RobustCCW(b,c,a) == RobustCCW(a,b,c) for all a,b,c (3) RobustCCW(c,b,a)
   * ==-RobustCCW(a,b,c) for all a,b,c
   *
   *  In other words:
   *
   *  (1) The result is zero if and only if two points are the same. (2)
   * Rotating the order of the arguments does not affect the result. (3)
   * Exchanging any two arguments inverts the result.
   *
   *  This function is essentially like taking the sign of the determinant of
   * a,b,c, except that it has additional logic to make sure that the above
   * properties hold even when the three points are coplanar, and to deal with
   * the limitations of floating-point arithmetic.
   *
   *  Note: a, b and c are expected to be of unit length. Otherwise, the results
   * are undefined.
   */
  public static int robustCCW(SpherePoint a, SpherePoint b, SpherePoint c) {
    return robustCCW(a, b, c, SpherePoint.crossProd(a, b));
  }

  /**
   * A more efficient version of RobustCCW that allows the precomputed
   * cross-product of A and B to be specified.
   *
   *  Note: a, b and c are expected to be of unit length. Otherwise, the results
   * are undefined
   */
  public static int robustCCW(SpherePoint a, SpherePoint b, SpherePoint c, SpherePoint aCrossB) {
    // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));

    // There are 14 multiplications and additions to compute the determinant
    // below. Since all three points are normalized, it is possible to show
    // that the average rounding error per operation does not exceed 2**-54,
    // the maximum rounding error for an operation whose result magnitude is in
    // the range [0.5,1). Therefore, if the absolute value of the determinant
    // is greater than 2*14*(2**-54), the determinant will have the same sign
    // even if the arguments are rotated (which produces a mathematically
    // equivalent result but with potentially different rounding errors).
    final double kMinAbsValue = 1.6e-15; // 2 * 14 * 2**-54

    double det = aCrossB.dotProd(c);

    // Double-check borderline cases in debug mode.
    // assert ((Math.abs(det) < kMinAbsValue) || (Math.abs(det) > 1000 * kMinAbsValue)
    //    || (det * expensiveCCW(a, b, c) > 0));

    if (det > kMinAbsValue) {
      return 1;
    }

    if (det < -kMinAbsValue) {
      return -1;
    }

    return expensiveCCW(a, b, c);
  }

  /**
   * A relatively expensive calculation invoked by RobustCCW() if the sign of
   * the determinant is uncertain.
   */
  private static int expensiveCCW(SpherePoint a, SpherePoint b, SpherePoint c) {
    // Return zero if and only if two points are the same. This ensures (1).
    if (a.equals(b) || b.equals(c) || c.equals(a)) {
      return 0;
    }

    // Now compute the determinant in a stable way. Since all three points are
    // unit length and we know that the determinant is very close to zero, this
    // means that points are very nearly colinear. Furthermore, the most common
    // situation is where two points are nearly identical or nearly antipodal.
    // To get the best accuracy in this situation, it is important to
    // immediately reduce the magnitude of the arguments by computing either
    // A+B or A-B for each pair of points. Note that even if A and B differ
    // only in their low bits, A-B can be computed very accurately. On the
    // other hand we can't accurately represent an arbitrary linear combination
    // of two vectors as would be required for Gaussian elimination. The code
    // below chooses the vertex opposite the longest edge as the "origin" for
    // the calculation, and computes the different vectors to the other two
    // vertices. This minimizes the sum of the lengths of these vectors.
    //
    // This implementation is very stable numerically, but it still does not
    // return consistent results in all cases. For example, if three points are
    // spaced far apart from each other along a great circle, the sign of the
    // result will basically be random (although it will still satisfy the
    // conditions documented in the header file). The only way to return
    // consistent results in all cases is to compute the result using
    // arbitrary-precision arithmetic. I considered using the Gnu MP library,
    // but this would be very expensive (up to 2000 bits of precision may be
    // needed to store the intermediate results) and seems like overkill for
    // this problem. The MP library is apparently also quite particular about
    // compilers and compilation options and would be a pain to maintain.

    // We want to handle the case of nearby points and nearly antipodal points
    // accurately, so determine whether A+B or A-B is smaller in each case.
    double sab = (a.dotProd(b) > 0) ? -1 : 1;
    double sbc = (b.dotProd(c) > 0) ? -1 : 1;
    double sca = (c.dotProd(a) > 0) ? -1 : 1;
    SpherePoint vab = SpherePoint.add(a, SpherePoint.mul(b, sab));
    SpherePoint vbc = SpherePoint.add(b, SpherePoint.mul(c, sbc));
    SpherePoint vca = SpherePoint.add(c, SpherePoint.mul(a, sca));
    double dab = vab.norm2();
    double dbc = vbc.norm2();
    double dca = vca.norm2();

    // Sort the difference vectors to find the longest edge, and use the
    // opposite vertex as the origin. If two difference vectors are the same
    // length, we break ties deterministically to ensure that the symmetry
    // properties guaranteed in the header file will be true.
    double sign;
    if (dca < dbc || (dca == dbc && a.lessThan(b))) {
      if (dab < dbc || (dab == dbc && a.lessThan(c))) {
        // The "sab" factor converts A +/- B into B +/- A.
        sign = SpherePoint.crossProd(vab, vca).dotProd(a) * sab; // BC is longest
                                                             // edge
      } else {
        sign = SpherePoint.crossProd(vca, vbc).dotProd(c) * sca; // AB is longest
                                                             // edge
      }
    } else {
      if (dab < dca || (dab == dca && b.lessThan(c))) {
        sign = SpherePoint.crossProd(vbc, vab).dotProd(b) * sbc; // CA is longest
                                                             // edge
      } else {
        sign = SpherePoint.crossProd(vca, vbc).dotProd(c) * sca; // AB is longest
                                                             // edge
      }
    }
    if (sign > 0) {
      return 1;
    }
    if (sign < 0) {
      return -1;
    }

    // The points A, B, and C are numerically indistinguishable from coplanar.
    // This may be due to roundoff error, or the points may in fact be exactly
    // coplanar. We handle this situation by perturbing all of the points by a
    // vector (eps, eps**2, eps**3) where "eps" is an infinitesmally small
    // positive number (e.g. 1 divided by a googolplex). The perturbation is
    // done symbolically, i.e. we compute what would happen if the points were
    // perturbed by this amount. It turns out that this is equivalent to
    // checking whether the points are ordered CCW around the origin first in
    // the Y-Z plane, then in the Z-X plane, and then in the X-Y plane.

    int ccw =
        planarOrderedCCW(new R2Vector(a.y, a.z), new R2Vector(b.y, b.z), new R2Vector(c.y, c.z));
    if (ccw == 0) {
      ccw =
          planarOrderedCCW(new R2Vector(a.z, a.x), new R2Vector(b.z, b.x), new R2Vector(c.z, c.x));
      if (ccw == 0) {
        ccw = planarOrderedCCW(
            new R2Vector(a.x, a.y), new R2Vector(b.x, b.y), new R2Vector(c.x, c.y));
        // assert (ccw != 0);
      }
    }
    return ccw;
  }


  public static int planarCCW(R2Vector a, R2Vector b) {
    // Return +1 if the edge AB is CCW around the origin, etc.
    double sab = (a.dotProd(b) > 0) ? -1 : 1;
    R2Vector vab = R2Vector.add(a, R2Vector.mul(b, sab));
    double da = a.norm2();
    double db = b.norm2();
    double sign;
    if (da < db || (da == db && a.lessThan(b))) {
      sign = a.crossProd(vab) * sab;
    } else {
      sign = vab.crossProd(b);
    }
    if (sign > 0) {
      return 1;
    }
    if (sign < 0) {
      return -1;
    }
    return 0;
  }

  public static int planarOrderedCCW(R2Vector a, R2Vector b, R2Vector c) {
    int sum = 0;
    sum += planarCCW(a, b);
    sum += planarCCW(b, c);
    sum += planarCCW(c, a);
    if (sum > 0) {
      return 1;
    }
    if (sum < 0) {
      return -1;
    }
    return 0;
  }

  /**
   * Return true if the edges OA, OB, and OC are encountered in that order while
   * sweeping CCW around the point O. You can think of this as testing whether
   * A <= B <= C with respect to a continuous CCW ordering around O.
   *
   * Properties:
   * <ol>
   *   <li>If orderedCCW(a,b,c,o) && orderedCCW(b,a,c,o), then a == b</li>
   *   <li>If orderedCCW(a,b,c,o) && orderedCCW(a,c,b,o), then b == c</li>
   *   <li>If orderedCCW(a,b,c,o) && orderedCCW(c,b,a,o), then a == b == c</li>
   *   <li>If a == b or b == c, then orderedCCW(a,b,c,o) is true</li>
   *   <li>Otherwise if a == c, then orderedCCW(a,b,c,o) is false</li>
   * </ol>
   */
  public static boolean orderedCCW(SpherePoint a, SpherePoint b, SpherePoint c, SpherePoint o) {
    // The last inequality below is ">" rather than ">=" so that we return true
    // if A == B or B == C, and otherwise false if A == C. Recall that
    // RobustCCW(x,y,z) == -RobustCCW(z,y,x) for all x,y,z.

    int sum = 0;
    if (robustCCW(b, o, a) >= 0) {
      ++sum;
    }
    if (robustCCW(c, o, b) >= 0) {
      ++sum;
    }
    if (robustCCW(a, o, c) > 0) {
      ++sum;
    }
    return sum >= 2;
  }

  /**
   * Return the angle at the vertex B in the triangle ABC. The return value is
   * always in the range [0, Pi]. The points do not need to be normalized.
   * Ensures that Angle(a,b,c) == Angle(c,b,a) for all a,b,c.
   *
   *  The angle is undefined if A or C is diametrically opposite from B, and
   * becomes numerically unstable as the length of edge AB or BC approaches 180
   * degrees.
   */
  public static double angle(SpherePoint a, SpherePoint b, SpherePoint c) {
    return SpherePoint.crossProd(a, b).angle(SpherePoint.crossProd(c, b));
  }

  /**
   * Return the exterior angle at the vertex B in the triangle ABC. The return
   * value is positive if ABC is counterclockwise and negative otherwise. If you
   * imagine an ant walking from A to B to C, this is the angle that the ant
   * turns at vertex B (positive = left, negative = right). Ensures that
   * TurnAngle(a,b,c) == -TurnAngle(c,b,a) for all a,b,c.
   *
   * @param a
   * @param b
   * @param c
   * @return the exterior angle at the vertex B in the triangle ABC
   */
  public static double turnAngle(SpherePoint a, SpherePoint b, SpherePoint c) {
    // This is a bit less efficient because we compute all 3 cross products, but
    // it ensures that turnAngle(a,b,c) == -turnAngle(c,b,a) for all a,b,c.
    double outAngle = SpherePoint.crossProd(b, a).angle(SpherePoint.crossProd(c, b));
    return (robustCCW(a, b, c) > 0) ? outAngle : -outAngle;
  }

  /**
   * Return true if two points are within the given distance of each other
   * (mainly useful for testing).
   */
  public static boolean approxEquals(SpherePoint a, SpherePoint b, double maxError) {
    return a.angle(b) <= maxError;
  }

  public static boolean approxEquals(SpherePoint a, SpherePoint b) {
    return approxEquals(a, b, 1e-15);
  }

  public static boolean approxEquals(double a, double b, double maxError) {
    return Math.abs(a - b) <= maxError;
  }

  public static boolean approxEquals(double a, double b) {
    return approxEquals(a, b, 1e-15);
  }
//
//  // Don't instantiate
//  private Sphere() {
//  }
}
