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

import com.google.common.base.Objects;
import com.google.common.base.Preconditions;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * 球面上的线类
 * An S2Polyline represents a sequence of zero or more vertices connected by
 * straight edges (geodesics). Edges of length 0 and 180 degrees are not
 * allowed, i.e. adjacent vertices should not be identical or antipodal.
 *
 * <p>Note: Polylines do not have a Contains(S2Point) method, because
 * "containment" is not numerically well-defined except at the polyline
 * vertices.
 *
 */
public strictfp class SpherePolyline implements SphereRegion {
  private static final Logger log = Logger.getLogger(SpherePolyline.class.getCanonicalName());

  private final int numVertices;
  private final SpherePoint[] vertices;

  /**
   * Create a polyline that connects the given vertices. Empty polylines are
   * allowed. Adjacent vertices should not be identical or antipodal. All
   * vertices should be unit length.
   */
  public SpherePolyline(List<SpherePoint> vertices) {
    // assert isValid(vertices);
    this.numVertices = vertices.size();
    this.vertices = vertices.toArray(new SpherePoint[numVertices]);
  }

  /**
   * Copy constructor.
   *
   * TODO(dbeaumont): Now that S2Polyline is immutable, remove this.
   */
  public SpherePolyline(SpherePolyline src) {
    this.numVertices = src.numVertices();
    this.vertices = src.vertices.clone();
  }

  /**
   * Return true if the given vertices form a valid polyline.
   */
  public boolean isValid(List<SpherePoint> vertices) {
    // All vertices must be unit length.
    int n = vertices.size();
    for (int i = 0; i < n; ++i) {
      if (!Sphere.isUnitLength(vertices.get(i))) {
        log.info("Vertex " + i + " is not unit length");
        return false;
      }
    }

    // Adjacent vertices must not be identical or antipodal.
    for (int i = 1; i < n; ++i) {
      if (vertices.get(i - 1).equals(vertices.get(i))
          || vertices.get(i - 1).equals(SpherePoint.neg(vertices.get(i)))) {
        log.info("Vertices " + (i - 1) + " and " + i + " are identical or antipodal");
        return false;
      }
    }

    return true;
  }

  public int numVertices() {
    return numVertices;
  }

  public SpherePoint vertex(int k) {
    // assert (k >= 0 && k < numVertices);
    return vertices[k];
  }

  /**
   * Return the angle corresponding to the total arclength of the polyline on a
   * unit sphere.
   */
  public S1Angle getArclengthAngle() {
    double lengthSum = 0;
    for (int i = 1; i < numVertices(); ++i) {
      lengthSum += vertex(i - 1).angle(vertex(i));
    }
    return S1Angle.radians(lengthSum);
  }

  /**
   * Return the point whose distance from vertex 0 along the polyline is the
   * given fraction of the polyline's total length. Fractions less than zero or
   * greater than one are clamped. The return value is unit length. This cost of
   * this function is currently linear in the number of vertices.
   */
  public SpherePoint interpolate(double fraction) {
    // We intentionally let the (fraction >= 1) case fall through, since
    // we need to handle it in the loop below in any case because of
    // possible roundoff errors.
    if (fraction <= 0) {
      return vertex(0);
    }

    double lengthSum = 0;
    for (int i = 1; i < numVertices(); ++i) {
      lengthSum += vertex(i - 1).angle(vertex(i));
    }
    double target = fraction * lengthSum;
    for (int i = 1; i < numVertices(); ++i) {
      double length = vertex(i - 1).angle(vertex(i));
      if (target < length) {
        // This code interpolates with respect to arc length rather than
        // straight-line distance, and produces a unit-length result.
        double f = Math.sin(target) / Math.sin(length);
        return SpherePoint.add(SpherePoint.mul(vertex(i - 1), (Math.cos(target) - f * Math.cos(length))),
            SpherePoint.mul(vertex(i), f));
      }
      target -= length;
    }
    return vertex(numVertices() - 1);
  }

  // S2Region interface (see {@code S2Region} for details):

  /** Return a bounding spherical cap. */
  @Override
  public SphereCap getCapBound() {
    return getRectBound().getCapBound();
  }


  /** Return a bounding latitude-longitude rectangle. */
  @Override
  public SphereLatLngRect getRectBound() {
    SphereEdgeUtil.RectBounder bounder = new SphereEdgeUtil.RectBounder();
    for (int i = 0; i < numVertices(); ++i) {
      bounder.addPoint(vertex(i));
    }
    return bounder.getBound();
  }


  /**
   * Given a point, returns the index of the start point of the (first) edge on
   * the polyline that is closest to the given point. The polyline must have at
   * least one vertex. Throws IllegalStateException if this is not the case.
   */
  public int getNearestEdgeIndex(SpherePoint point) {
    Preconditions.checkState(numVertices() > 0, "Empty polyline");

    if (numVertices() == 1) {
      // If there is only one vertex, the "edge" is trivial, and it's the only one
      return 0;
    }

    // Initial value larger than any possible distance on the unit sphere.
    S1Angle minDistance = S1Angle.radians(10);
    int minIndex = -1;

    // Find the line segment in the polyline that is closest to the point given.
    for (int i = 0; i < numVertices() - 1; ++i) {
      S1Angle distanceToSegment = SphereEdgeUtil.getDistance(point, vertex(i), vertex(i + 1));
      if (distanceToSegment.lessThan(minDistance)) {
        minDistance = distanceToSegment;
        minIndex = i;
      }
    }
    return minIndex;
  }

  /**
   * Given a point p and the index of the start point of an edge of this polyline,
   * returns the point on that edge that is closest to p.
   */
  public SpherePoint projectToEdge(SpherePoint point, int index) {
    Preconditions.checkState(numVertices() > 0, "Empty polyline");
    Preconditions.checkState(numVertices() == 1 || index < numVertices() - 1, "Invalid edge index");
    if (numVertices() == 1) {
      // If there is only one vertex, it is always closest to any given point.
      return vertex(0);
    }
    return SphereEdgeUtil.getClosestPoint(point, vertex(index), vertex(index + 1));
  }

  @Override
  public boolean equals(Object that) {
    if (!(that instanceof SpherePolyline)) {
      return false;
    }

    SpherePolyline thatPolygon = (SpherePolyline) that;
    if (numVertices != thatPolygon.numVertices) {
      return false;
    }

    for (int i = 0; i < vertices.length; i++) {
      if (!vertices[i].equals(thatPolygon.vertices[i])) {
        return false;
      }
    }
    return true;
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(numVertices, Arrays.deepHashCode(vertices));
  }
}
