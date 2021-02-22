/*
 * Copyright 2006 Google Inc.
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

/**
 * （单位）球面上的点
 * An S2Point represents a point on the unit sphere as a 3D vector. Usually
 * points are normalized to be unit length, but some methods do not require
 * this.
 *
 */
public strictfp class SpherePoint implements Comparable<SpherePoint> {
  // coordinates of the points
  final public double x;
  final public double y;
  final public double z;

  public SpherePoint() {
    x = y = z = 0;
  }

  public SpherePoint(double x, double y, double z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public static SpherePoint minus(SpherePoint p1, SpherePoint p2) {
    return sub(p1, p2);
  }

  public static SpherePoint neg(SpherePoint p) {
    return new SpherePoint(-p.x, -p.y, -p.z);
  }

  public double norm2() {
    return x * x + y * y + z * z;
  }

  public double norm() {
    return Math.sqrt(norm2());
  }

  //Cross Product， 返回两个向量构成平面的的法向量
  public static SpherePoint crossProd(final SpherePoint p1, final SpherePoint p2) {
    return new SpherePoint(
        p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x);
  }

  public static SpherePoint add(final SpherePoint p1, final SpherePoint p2) {
    return new SpherePoint(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
  }

  public static SpherePoint sub(final SpherePoint p1, final SpherePoint p2) {
    return new SpherePoint(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
  }

  //Dot Product，表征或计算两个向量之间的夹角 b向量在a向量方向上的投影，当点积=0时表示两者正交，当点积=|a||b|时，两者共线
  public double dotProd(SpherePoint that) {
    return this.x * that.x + this.y * that.y + this.z * that.z;
  }

  public static SpherePoint mul(final SpherePoint p, double m) {
    return new SpherePoint(m * p.x, m * p.y, m * p.z);
  }

  public static SpherePoint div(final SpherePoint p, double m) {
    return new SpherePoint(p.x / m, p.y / m, p.z / m);
  }

  /** return a vector orthogonal to this one */
  public SpherePoint ortho() {
    int k = largestAbsComponent();
    SpherePoint temp;
    if (k == 1) {
      temp = new SpherePoint(1, 0, 0);
    } else if (k == 2) {
      temp = new SpherePoint(0, 1, 0);
    } else {
      temp = new SpherePoint(0, 0, 1);
    }
    return SpherePoint.normalize(crossProd(this, temp));
  }

  /** Return the index of the largest component fabs */
  public int largestAbsComponent() {
    SpherePoint temp = fabs(this);
    if (temp.x > temp.y) {
      if (temp.x > temp.z) {
        return 0;
      } else {
        return 2;
      }
    } else {
      if (temp.y > temp.z) {
        return 1;
      } else {
        return 2;
      }
    }
  }

  public static SpherePoint fabs(SpherePoint p) {
    return new SpherePoint(Math.abs(p.x), Math.abs(p.y), Math.abs(p.z));
  }

  public static SpherePoint normalize(SpherePoint p) {
    double norm = p.norm();
    if (norm != 0) {
      norm = 1.0 / norm;
    }
    return SpherePoint.mul(p, norm);
  }

  public double get(int axis) {
    return (axis == 0) ? x : (axis == 1) ? y : z;
  }

  /** Return the angle between two vectors in radians */
  public double angle(SpherePoint va) {
    return Math.atan2(crossProd(this, va).norm(), this.dotProd(va));
  }

  /**
   * Compare two vectors, return true if all their components are within a
   * difference of margin.
   */
  boolean aequal(SpherePoint that, double margin) {
    return (Math.abs(x - that.x) < margin) && (Math.abs(y - that.y) < margin)
        && (Math.abs(z - that.z) < margin);
  }

  @Override
  public boolean equals(Object that) {
    if (!(that instanceof SpherePoint)) {
      return false;
    }
    SpherePoint thatPoint = (SpherePoint) that;
    return this.x == thatPoint.x && this.y == thatPoint.y && this.z == thatPoint.z;
  }

  public boolean lessThan(SpherePoint vb) {
    if (x < vb.x) {
      return true;
    }
    if (vb.x < x) {
      return false;
    }
    if (y < vb.y) {
      return true;
    }
    if (vb.y < y) {
      return false;
    }
    if (z < vb.z) {
      return true;
    }
    return false;
  }

  // Required for Comparable
  @Override
  public int compareTo(SpherePoint other) {
    return (lessThan(other) ? -1 : (equals(other) ? 0 : 1));
  }

  @Override
  public String toString() {
    return "(" + x + ", " + y + ", " + z + ")";
  }

  public String toDegreesString() {
    SphereLatLng sphereLatLng = new SphereLatLng(this);
    return "(" + Double.toString(sphereLatLng.latDegrees()) + ", "
        + Double.toString(sphereLatLng.lngDegrees()) + ")";
  }

  /**
   * Calcualates hashcode based on stored coordinates. Since we want +0.0 and
   * -0.0 to be treated the same, we ignore the sign of the coordinates.
   */
  @Override
  public int hashCode() {
    long value = 17;
    value += 37 * value + Double.doubleToLongBits(Math.abs(x));
    value += 37 * value + Double.doubleToLongBits(Math.abs(y));
    value += 37 * value + Double.doubleToLongBits(Math.abs(z));
    return (int) (value ^ (value >>> 32));
  }
}
