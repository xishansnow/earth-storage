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

import com.google.common.collect.*;

import java.util.*;
import java.util.logging.Logger;

/**
 * 球面多边形类
 * An S2Polygon is an S2Region object that represents a polygon. A polygon
 * consists of zero or more {@link SphereLoop loops} representing "shells" and
 * "holes". All loops should be oriented CCW, i.e. the shell or hole is on the
 * left side of the loop. Loops may be specified in any order. A point is
 * defined to be inside the polygon if it is contained by an odd number of
 * loops.
 *
 *  Polygons have the following restrictions:
 *
 *  - Loops may not cross, i.e. the boundary of a loop may not intersect both
 * the interior and exterior of any other loop.
 *
 *  - Loops may not share edges, i.e. if a loop contains an edge AB, then no
 * other loop may contain AB or BA.
 *
 *  - No loop may cover more than half the area of the sphere. This ensures that
 * no loop properly contains the complement of any other loop, even if the loops
 * are from different polygons. (Loops that represent exact hemispheres are
 * allowed.)
 *
 *  Loops may share vertices, however no vertex may appear twice in a single
 * loop.
 *
 */
public final strictfp class SpherePolygon implements SphereRegion, Comparable<SpherePolygon> {
  private static final Logger log = Logger.getLogger(SpherePolygon.class.getCanonicalName());

  private List<SphereLoop> loops;

  private SphereLatLngRect bound;
  private boolean hasHoles;
  private int numVertices;

  /**
   * Creates an empty polygon that should be initialized by calling Init().
   */
  public SpherePolygon() {
    this.loops = Lists.newArrayList();
    this.bound = SphereLatLngRect.empty();
    this.hasHoles = false;
    this.numVertices = 0;
  }

  /**
   * Convenience constructor that calls Init() with the given loops. Clears the
   * given list.
   */
  public SpherePolygon(List<SphereLoop> loops) {
    this.loops = Lists.newArrayList();
    this.bound = SphereLatLngRect.empty();

    init(loops);
  }

  /**
   * Copy constructor.
   */
  public SpherePolygon(SphereLoop loop) {
    this.loops = Lists.newArrayList();
    this.bound = loop.getRectBound();
    this.hasHoles = false;
    this.numVertices = loop.numVertices();

    loops.add(loop);
  }

  /**
   * Copy constructor.
   */
  public SpherePolygon(SpherePolygon src) {
    this.loops = Lists.newArrayList();
    this.bound = src.getRectBound();
    this.hasHoles = src.hasHoles;
    this.numVertices = src.numVertices;

    for (int i = 0; i < src.numLoops(); ++i) {
      loops.add(new SphereLoop(src.loop(i)));
    }
  }

  /**
   * Comparator (needed by Comparable interface). For two polygons to be
   * compared as equal: - the must have the same number of loops; - the loops
   * must be ordered in the same way (this is guaranteed by the total ordering
   * imposed by sortValueLoops). - loops must be logically equivalent (even if
   * ordered with a different starting point, e.g. ABCD and BCDA).
   */
  @Override
  public int compareTo(SpherePolygon other) {
    // If number of loops differ, use that.
    if (this.numLoops() != other.numLoops()) {
      return this.numLoops() - other.numLoops();
    }
    for (int i = 0; i < this.numLoops(); ++i) {
      int compare = this.loops.get(i).compareTo(other.loops.get(i));
      if (compare != 0) {
        return compare;
      }
    }
    return 0;
  }

  /**
   * Initialize a polygon by taking ownership of the given loops and clearing
   * the given list. This method figures out the loop nesting hierarchy and then
   * reorders the loops by following a preorder traversal. This implies that
   * each loop is immediately followed by its descendants in the nesting
   * hierarchy. (See also getParent and getLastDescendant.)
   */
  public void init(List<SphereLoop> loops) {
    // assert isValid(loops);
    // assert (this.loops.isEmpty());

    Map<SphereLoop, List<SphereLoop>> loopMap = Maps.newHashMap();
    // Yes, a null key is valid. It is used here to refer to the root of the
    // loopMap
    loopMap.put(null, Lists.<SphereLoop>newArrayList());

    for (SphereLoop loop : loops) {
      insertLoop(loop, null, loopMap);
      this.numVertices += loop.numVertices();
    }
    loops.clear();

    // Sort all of the lists of loops; in this way we guarantee a total ordering
    // on loops in the polygon. Loops will be sorted by their natural ordering,
    // while also preserving the requirement that each loop is immediately
    // followed by its descendants in the nesting hierarchy.
    //
    // TODO(andriy): as per kirilll in CL 18750833 code review comments:
    // This should work for now, but I think it's possible to guarantee the
    // correct order inside insertLoop by searching for the correct position in
    // the children list before inserting.
    sortValueLoops(loopMap);

    // Reorder the loops in depth-first traversal order.
    // Starting at null == starting at the root
    initLoop(null, -1, loopMap);

    // TODO(dbeaumont): Add tests or preconditions for these asserts (here and elesewhere).
    // forall i != j : containsChild(loop(i), loop(j), loopMap) == loop(i).containsNested(loop(j)));

    // Compute the bounding rectangle of the entire polygon.
    hasHoles = false;
    bound = SphereLatLngRect.empty();
    for (int i = 0; i < numLoops(); ++i) {
      if (loop(i).sign() < 0) {
        hasHoles = true;
      } else {
        bound = bound.union(loop(i).getRectBound());
      }
    }
  }

  /**
   * Release ownership of the loops of this polygon by appending them to the
   * given list. Resets the polygon to be empty.
   */
  public void release(List<SphereLoop> loops) {
    loops.addAll(this.loops);
    this.loops.clear();
    bound = SphereLatLngRect.empty();
    hasHoles = false;
    numVertices = 0;
  }

  /**
   * Return true if the given loops form a valid polygon. Assumes that that all
   * of the given loops have already been validated.
   */
  public static boolean isValid(final List<SphereLoop> loops) {
    // If a loop contains an edge AB, then no other loop may contain AB or BA.
    // We only need this test if there are at least two loops, assuming that
    // each loop has already been validated.
    if (loops.size() > 1) {
      Map<UndirectedEdge, LoopVertexIndexPair> edges = Maps.newHashMap();
      for (int i = 0; i < loops.size(); ++i) {
        SphereLoop lp = loops.get(i);
        for (int j = 0; j < lp.numVertices(); ++j) {
          UndirectedEdge key = new UndirectedEdge(lp.vertex(j), lp.vertex(j + 1));
          LoopVertexIndexPair value = new LoopVertexIndexPair(i, j);
          if (edges.containsKey(key)) {
            LoopVertexIndexPair other = edges.get(key);
            log.info(
                "Duplicate edge: loop " + i + ", edge " + j + " and loop " + other.getLoopIndex()
                    + ", edge " + other.getVertexIndex());
            return false;
          } else {
            edges.put(key, value);
          }
        }
      }
    }

    // Verify that no loop covers more than half of the sphere, and that no
    // two loops cross.
    for (int i = 0; i < loops.size(); ++i) {
      if (!loops.get(i).isNormalized()) {
        log.info("Loop " + i + " encloses more than half the sphere");
        return false;
      }
      for (int j = i + 1; j < loops.size(); ++j) {
        // This test not only checks for edge crossings, it also detects
        // cases where the two boundaries cross at a shared vertex.
        if (loops.get(i).containsOrCrosses(loops.get(j)) < 0) {
          log.info("Loop " + i + " crosses loop " + j);
          return false;
        }
      }
    }
    return true;
  }

  public int numLoops() {
    return loops.size();
  }

  public SphereLoop loop(int k) {
    return loops.get(k);
  }

  /**
   * Return the index of the parent of loop k, or -1 if it has no parent.
   */
  public int getParent(int k) {
    int depth = loop(k).depth();
    if (depth == 0) {
      return -1; // Optimization.
    }
    while (--k >= 0 && loop(k).depth() >= depth) {
      // spin
    }
    return k;
  }

  /**
   * Return the index of the last loop that is contained within loop k. Returns
   * num_loops() - 1 if k < 0. Note that loops are indexed according to a
   * preorder traversal of the nesting hierarchy, so the immediate children of
   * loop k can be found by iterating over loops (k+1)..getLastDescendant(k) and
   * selecting those whose depth is equal to (loop(k).depth() + 1).
   */
  public int getLastDescendant(int k) {
    if (k < 0) {
      return numLoops() - 1;
    }
    int depth = loop(k).depth();
    while (++k < numLoops() && loop(k).depth() > depth) {
      // spin
    }
    return k - 1;
  }

  private SphereAreaCentroid getAreaCentroid(boolean doCentroid) {
    double areaSum = 0;
    SpherePoint centroidSum = new SpherePoint(0, 0, 0);
    for (int i = 0; i < numLoops(); ++i) {
      SphereAreaCentroid areaCentroid = doCentroid ? loop(i).getAreaAndCentroid() : null;
      double loopArea = doCentroid ? areaCentroid.getArea() : loop(i).getArea();

      int loopSign = loop(i).sign();
      areaSum += loopSign * loopArea;
      if (doCentroid) {
        SpherePoint currentCentroid = areaCentroid.getCentroid();
        centroidSum =
            new SpherePoint(centroidSum.x + loopSign * currentCentroid.x,
                centroidSum.y + loopSign * currentCentroid.y,
                centroidSum.z + loopSign * currentCentroid.z);
      }
    }

    return new SphereAreaCentroid(areaSum, doCentroid ? centroidSum : null);
  }

  /**
   * Return the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops (this value return value is between 0 and 4*Pi)
   * and the true centroid of the polygon multiplied by the area of the polygon
   * (see s2geometry.h for details on centroids). Note that the centroid may not be
   * contained by the polygon.
   */
  public SphereAreaCentroid getAreaAndCentroid() {
    return getAreaCentroid(true);
  }

  /**
   * Return the area of the polygon interior, i.e. the region on the left side
   * of an odd number of loops. The return value is between 0 and 4*Pi.
   */
  public double getArea() {
    return getAreaCentroid(false).getArea();
  }

  /**
   * Return the true centroid of the polygon multiplied by the area of the
   * polygon (see s2geometry.h for details on centroids). Note that the centroid may not
   * be contained by the polygon.
   */
  public SpherePoint getCentroid() {
    return getAreaCentroid(true).getCentroid();
  }

  /**
   * Returns the shortest distance from a point P to this polygon, given as the
   * angle formed between P, the origin and the nearest point on the polygon to
   * P. This angle in radians is equivalent to the arclength along the unit
   * sphere.
   *
   * If the point is contained inside the polygon, the distance returned is 0.
   */
  public S1Angle getDistance(SpherePoint p) {
    if (contains(p)) {
      return S1Angle.radians(0);
    }

    // The furthest point from p on the sphere is its antipode, which is an
    // angle of PI radians. This is an upper bound on the angle.
    S1Angle minDistance = S1Angle.radians(Math.PI);
    for (int i = 0; i < numLoops(); i++) {
      minDistance = S1Angle.min(minDistance, loop(i).getDistance(p));
    }

    return minDistance;
  }


  /**
   * Return true if this polygon contains the given other polygon, i.e. if
   * polygon A contains all points contained by polygon B.
   */
  public boolean contains(SpherePolygon b) {
    // If both polygons have one loop, use the more efficient S2Loop method.
    // Note that S2Loop.contains does its own bounding rectangle check.
    if (numLoops() == 1 && b.numLoops() == 1) {
      return loop(0).contains(b.loop(0));
    }

    // Otherwise if neither polygon has holes, we can still use the more
    // efficient S2Loop::Contains method (rather than ContainsOrCrosses),
    // but it's worthwhile to do our own bounds check first.
    if (!bound.contains(b.getRectBound())) {
      // If the union of the bounding boxes spans the full longitude range,
      // it is still possible that polygon A contains B. (This is only
      // possible if at least one polygon has multiple shells.)
      if (!bound.lng().union(b.getRectBound().lng()).isFull()) {
        return false;
      }
    }
    if (!hasHoles && !b.hasHoles) {
      for (int j = 0; j < b.numLoops(); ++j) {
        if (!anyLoopContains(b.loop(j))) {
          return false;
        }
      }
      return true;
    }

    // This could be implemented more efficiently for polygons with lots of
    // holes by keeping a copy of the LoopMap computed during initialization.
    // However, in practice most polygons are one loop, and multiloop polygons
    // tend to consist of many shells rather than holes. In any case, the real
    // way to get more efficiency is to implement a sub-quadratic algorithm
    // such as building a trapezoidal map.

    // Every shell of B must be contained by an odd number of loops of A,
    // and every hole of A must be contained by an even number of loops of B.
    return containsAllShells(b) && b.excludesAllHoles(this);
  }

  /**
   * Return true if this polygon intersects the given other polygon, i.e. if
   * there is a point that is contained by both polygons.
   */
  public boolean intersects(SpherePolygon b) {
    // A.intersects(B) if and only if !complement(A).contains(B). However,
    // implementing a complement() operation is trickier than it sounds,
    // and in any case it's more efficient to test for intersection directly.

    // If both polygons have one loop, use the more efficient S2Loop method.
    // Note that S2Loop.intersects does its own bounding rectangle check.
    if (numLoops() == 1 && b.numLoops() == 1) {
      return loop(0).intersects(b.loop(0));
    }

    // Otherwise if neither polygon has holes, we can still use the more
    // efficient S2Loop.intersects method. The polygons intersect if and
    // only if some pair of loop regions intersect.
    if (!bound.intersects(b.getRectBound())) {
      return false;
    }
    if (!hasHoles && !b.hasHoles) {
      for (int i = 0; i < numLoops(); ++i) {
        for (int j = 0; j < b.numLoops(); ++j) {
          if (loop(i).intersects(b.loop(j))) {
            return true;
          }
        }
      }
      return false;
    }

    // Otherwise if any shell of B is contained by an odd number of loops of A,
    // or any shell of A is contained by an odd number of loops of B, there is
    // an intersection.
    return intersectsAnyShell(b) || b.intersectsAnyShell(this);
  }

  private static void addIntersection(SpherePoint a0,
                                      SpherePoint a1,
                                      SpherePoint b0,
                                      SpherePoint b1,
                                      boolean addSharedEdges,
                                      int crossing,
                                      List<ParametrizedS2Point> intersections) {
    if (crossing > 0) {
      // There is a proper edge crossing.
      SpherePoint x = SphereEdgeUtil.getIntersection(a0, a1, b0, b1);
      double t = SphereEdgeUtil.getDistanceFraction(x, a0, a1);
      intersections.add(new ParametrizedS2Point(t, x));
    } else if (SphereEdgeUtil.vertexCrossing(a0, a1, b0, b1)) {
      // There is a crossing at one of the vertices. The basic rule is simple:
      // if a0 equals one of the "b" vertices, the crossing occurs at t=0;
      // otherwise, it occurs at t=1.
      //
      // This has the effect that when two symmetric edges are encountered (an
      // edge an its reverse), neither one is included in the output. When two
      // duplicate edges are encountered, both are included in the output. The
      // "addSharedEdges" flag allows one of these two copies to be removed by
      // changing its intersection parameter from 0 to 1.
      double t = (a0 == b0 || a0 == b1) ? 0 : 1;
      if (!addSharedEdges && a1 == b1) {
        t = 1;
      }
      intersections.add(new ParametrizedS2Point(t, t == 0 ? a0 : a1));
    }
  }

  /**
   * Returns total number of vertices in all loops.
   */
  public int getNumVertices() {
    return this.numVertices;
  }


  public boolean isNormalized() {
    Multiset<SpherePoint> vertices = HashMultiset.<SpherePoint>create();
    SphereLoop lastParent = null;
    for (int i = 0; i < numLoops(); ++i) {
      SphereLoop child = loop(i);
      if (child.depth() == 0) {
        continue;
      }
      SphereLoop parent = loop(getParent(i));
      if (parent != lastParent) {
        vertices.clear();
        for (int j = 0; j < parent.numVertices(); ++j) {
          vertices.add(parent.vertex(j));
        }
        lastParent = parent;
      }
      int count = 0;
      for (int j = 0; j < child.numVertices(); ++j) {
        if (vertices.count(child.vertex(j)) > 0) {
          ++count;
        }
      }
      if (count > 1) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return true if two polygons have the same boundary except for vertex
   * perturbations. Both polygons must have loops with the same cyclic vertex
   * order and the same nesting hierarchy, but the vertex locations are allowed
   * to differ by up to "max_error". Note: This method mostly useful only for
   * testing purposes.
   */
  boolean boundaryApproxEquals(SpherePolygon b, double maxError) {
    if (numLoops() != b.numLoops()) {
      log.severe(
          "!= loops: " + Integer.toString(numLoops()) + " vs. " + Integer.toString(b.numLoops()));
      return false;
    }

    // For now, we assume that there is at most one candidate match for each
    // loop. (So far this method is just used for testing.)
    for (int i = 0; i < numLoops(); ++i) {
      SphereLoop aLoop = loop(i);
      boolean success = false;
      for (int j = 0; j < numLoops(); ++j) {
        SphereLoop bLoop = b.loop(j);
        if (bLoop.depth() == aLoop.depth() && bLoop.boundaryApproxEquals(aLoop, maxError)) {
          success = true;
          break;
        }
      }
      if (!success) {
        return false;
      }
    }
    return true;
  }

  // S2Region interface (see S2Region.java for details):

  /** Return a bounding spherical cap. */
  @Override
  public SphereCap getCapBound() {
    return bound.getCapBound();
  }


  /** Return a bounding latitude-longitude rectangle. */
  @Override
  public SphereLatLngRect getRectBound() {
    return bound;
  }


  /**
   * The point 'p' does not need to be normalized.
   */
  public boolean contains(SpherePoint p) {
    if (numLoops() == 1) {
      return loop(0).contains(p); // Optimization.
    }
    if (!bound.contains(p)) {
      return false;
    }
    boolean inside = false;
    for (int i = 0; i < numLoops(); ++i) {
      inside ^= loop(i).contains(p);
      if (inside && !hasHoles) {
        break; // Shells are disjoint.
      }
    }
    return inside;
  }

  // For each map entry, sorts the value list.
  private static void sortValueLoops(Map<SphereLoop, List<SphereLoop>> loopMap) {
    for (SphereLoop key : loopMap.keySet()) {
      Collections.sort(loopMap.get(key));
    }
  }

  private static void insertLoop(SphereLoop newLoop, SphereLoop parent, Map<SphereLoop, List<SphereLoop>> loopMap) {
    List<SphereLoop> children = loopMap.get(parent);

    if (children == null) {
      children = Lists.newArrayList();
      loopMap.put(parent, children);
    }

    for (SphereLoop child : children) {
      if (child.containsNested(newLoop)) {
        insertLoop(newLoop, child, loopMap);
        return;
      }
    }

    // No loop may contain the complement of another loop. (Handling this case
    // is significantly more complicated.)
    // assert (parent == null || !newLoop.containsNested(parent));

    // Some of the children of the parent loop may now be children of
    // the new loop.
    List<SphereLoop> newChildren = loopMap.get(newLoop);
    for (int i = 0; i < children.size();) {
      SphereLoop child = children.get(i);
      if (newLoop.containsNested(child)) {
        if (newChildren == null) {
          newChildren = Lists.newArrayList();
          loopMap.put(newLoop, newChildren);
        }
        newChildren.add(child);
        children.remove(i);
      } else {
        ++i;
      }
    }
    children.add(newLoop);
  }

  private void initLoop(SphereLoop loop, int depth, Map<SphereLoop, List<SphereLoop>> loopMap) {
    if (loop != null) {
      loop.setDepth(depth);
      loops.add(loop);
    }
    List<SphereLoop> children = loopMap.get(loop);
    if (children != null) {
      for (SphereLoop child : children) {
        initLoop(child, depth + 1, loopMap);
      }
    }
  }

  private int containsOrCrosses(SphereLoop b) {
    boolean inside = false;
    for (int i = 0; i < numLoops(); ++i) {
      int result = loop(i).containsOrCrosses(b);
      if (result < 0) {
        return -1; // The loop boundaries intersect.
      }
      if (result > 0) {
        inside ^= true;
      }
    }
    return inside ? 1 : 0; // True if loop B is contained by the polygon.
  }

  /** Return true if any loop contains the given loop. */
  private boolean anyLoopContains(SphereLoop b) {
    for (int i = 0; i < numLoops(); ++i) {
      if (loop(i).contains(b)) {
        return true;
      }
    }
    return false;
  }

  /** Return true if this polygon (A) contains all the shells of B. */
  private boolean containsAllShells(SpherePolygon b) {
    for (int j = 0; j < b.numLoops(); ++j) {
      if (b.loop(j).sign() < 0) {
        continue;
      }
      if (containsOrCrosses(b.loop(j)) <= 0) {
        // Shell of B is not contained by A, or the boundaries intersect.
        return false;
      }
    }
    return true;
  }

  /**
   * Return true if this polygon (A) excludes (i.e. does not intersect) all
   * holes of B.
   */
  private boolean excludesAllHoles(SpherePolygon b) {
    for (int j = 0; j < b.numLoops(); ++j) {
      if (b.loop(j).sign() > 0) {
        continue;
      }
      if (containsOrCrosses(b.loop(j)) != 0) {
        // Hole of B is contained by A, or the boundaries intersect.
        return false;
      }
    }
    return true;
  }

  /** Return true if this polygon (A) intersects any shell of B. */
  private boolean intersectsAnyShell(SpherePolygon b) {
    for (int j = 0; j < b.numLoops(); ++j) {
      if (b.loop(j).sign() < 0) {
        continue;
      }
      if (containsOrCrosses(b.loop(j)) != 0) {
        // Shell of B is contained by A, or the boundaries intersect.
        return true;
      }
    }
    return false;
  }

  /**
   * A human readable representation of the polygon
   */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("Polygon: (").append(numLoops()).append(") loops:\n");
    for (int i = 0; i < numLoops(); ++i) {
      SphereLoop s2Loop = loop(i);
      sb.append("loop <\n");
      for (int v = 0; v < s2Loop.numVertices(); ++v) {
        SpherePoint spherePoint = s2Loop.vertex(v);
        sb.append(spherePoint.toDegreesString());
        sb.append("\n"); // end of vertex
      }
      sb.append(">\n"); // end of loop
    }
    return sb.toString();
  }

  private static final class UndirectedEdge {
    // Note: An UndirectedEdge and an S2Edge can never be considered equal (in
    // terms of the equals() method) and hence they re not be related types.
    // If you need to convert between the types then separate conversion
    // methods should be introduced.

    private final SpherePoint a;
    private final SpherePoint b;

    public UndirectedEdge(SpherePoint start, SpherePoint end) {
      this.a = start;
      this.b = end;
    }

    public SpherePoint getStart() {
      return a;
    }

    public SpherePoint getEnd() {
      return b;
    }

    @Override
    public String toString() {
      return String.format("Edge: (%s <-> %s)\n   or [%s <-> %s]",
          a.toDegreesString(), b.toDegreesString(), a, b);
    }

    @Override
    public boolean equals(Object o) {
      if (o == null || !(o instanceof UndirectedEdge)) {
        return false;
      }
      UndirectedEdge other = (UndirectedEdge) o;
      return ((getStart().equals(other.getStart()) && getEnd().equals(other.getEnd()))
          || (getStart().equals(other.getEnd()) && getEnd().equals(other.getStart())));
    }

    @Override
    public int hashCode() {
      return getStart().hashCode() + getEnd().hashCode();
    }
  }

  private static final class LoopVertexIndexPair {
    private final int loopIndex;
    private final int vertexIndex;

    public LoopVertexIndexPair(int loopIndex, int vertexIndex) {
      this.loopIndex = loopIndex;
      this.vertexIndex = vertexIndex;
    }

    public int getLoopIndex() {
      return loopIndex;
    }

    public int getVertexIndex() {
      return vertexIndex;
    }
  }

  /**
   * An S2Point that also has a parameter associated with it, which corresponds
   * to a time-like order on the points.
   */
  private static final class ParametrizedS2Point implements Comparable<ParametrizedS2Point> {
    private final double time;
    private final SpherePoint point;

    public ParametrizedS2Point(double time, SpherePoint point) {
      this.time = time;
      this.point = point;
    }

    public double getTime() {
      return time;
    }

    public SpherePoint getPoint() {
      return point;
    }

    @Override
    public int compareTo(ParametrizedS2Point o) {
      int compareTime = Double.compare(time, o.time);
      if (compareTime != 0) {
        return compareTime;
      }
      return point.compareTo(o.point);
    }
  }
}
