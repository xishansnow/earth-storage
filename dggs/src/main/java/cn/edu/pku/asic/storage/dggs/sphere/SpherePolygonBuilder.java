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

import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * 基于球面有向边，组装球面多边形的工具类
 * This is a simple class for assembling polygons out of edges. It requires that
 * no two edges cross. It can handle both directed and undirected edges, and
 * optionally it can also remove duplicate edge pairs (consisting of two
 * identical edges or an edge and its reverse edge). This is useful for
 * computing seamless unions of polygons that have been cut into pieces.
 *
 *  Here are some of the situations this class was designed to handle:
 *
 *  1. Computing the union of disjoint polygons that may share part of their
 * boundaries. For example, reassembling a lake that has been split into two
 * loops by a state boundary.
 *
 *  2. Constructing polygons from input data that does not follow S2
 * conventions, i.e. where loops may have repeated vertices, or distinct loops
 * may share edges, or shells and holes have opposite or unspecified
 * orientations.
 *
 *  3. Computing the symmetric difference of a set of polygons whose edges
 * intersect only at vertices. This can be used to implement a limited form of
 * polygon intersection or subtraction as well as unions.
 *
 *  4. As a tool for implementing other polygon operations by generating a
 * collection of directed edges and then assembling them into loops.
 *
 */
public strictfp class SpherePolygonBuilder {
  private static final Logger log = Logger.getLogger(SpherePolygonBuilder.class.getCanonicalName());

  private Options options;

  /**
   * The current set of edges, grouped by origin. The set of destination
   * vertices is a multiset so that the same edge can be present more than once.
   */
  private Map<SpherePoint, Multiset<SpherePoint>> edges;

  /**
   * Default constructor for well-behaved polygons. Uses the DIRECTED_XOR
   * options.
   */
  public SpherePolygonBuilder() {
    this(Options.DIRECTED_XOR);

  }

  public SpherePolygonBuilder(Options options) {
    this.options = options;
    this.edges = Maps.newHashMap();
  }

  public enum Options {

    /**
     * These are the options that should be used for assembling well-behaved
     * input data into polygons. All edges should be directed such that "shells"
     * and "holes" have opposite orientations (typically CCW shells and
     * clockwise holes), unless it is known that shells and holes do not share
     * any edges.
     */
    DIRECTED_XOR(false, true),

    /**
     * These are the options that should be used for assembling polygons that do
     * not follow the conventions above, e.g. where edge directions may vary
     * within a single loop, or shells and holes are not oppositely oriented.
     */
    UNDIRECTED_XOR(true, true),

    /**
     * These are the options that should be used for assembling edges where the
     * desired output is a collection of loops rather than a polygon, and edges
     * may occur more than once. Edges are treated as undirected and are not
     * XORed together, in particular, adding edge A->B also adds B->A.
     */
    UNDIRECTED_UNION(true, false),

    /**
     * Finally, select this option when the desired output is a collection of
     * loops rather than a polygon, but your input edges are directed and you do
     * not want reverse edges to be added implicitly as above.
     */
    DIRECTED_UNION(false, false);

    private boolean undirectedEdges;
    private boolean xorEdges;
    private boolean validate;
    private S1Angle mergeDistance;

    private Options(boolean undirectedEdges, boolean xorEdges) {
      this.undirectedEdges = undirectedEdges;
      this.xorEdges = xorEdges;
      this.validate = false;
      this.mergeDistance = S1Angle.radians(0);
    }

    /**
     * If "undirected_edges" is false, then the input is assumed to consist of
     * edges that can be assembled into oriented loops without reversing any of
     * the edges. Otherwise, "undirected_edges" should be set to true.
     */
    public boolean getUndirectedEdges() {
      return undirectedEdges;
    }

    /**
     * If "xor_edges" is true, then any duplicate edge pairs are removed. This
     * is useful for computing the union of a collection of polygons whose
     * interiors are disjoint but whose boundaries may share some common edges
     * (e.g. computing the union of South Africa, Lesotho, and Swaziland).
     *
     *  Note that for directed edges, a "duplicate edge pair" consists of an
     * edge and its corresponding reverse edge. This means that either (a)
     * "shells" and "holes" must have opposite orientations, or (b) shells and
     * holes do not share edges. Otherwise undirected_edges() should be
     * specified.
     *
     *  There are only two reasons to turn off xor_edges():
     *
     *  (1) assemblePolygon() will be called, and you want to assert that there
     * are no duplicate edge pairs in the input.
     *
     *  (2) assembleLoops() will be called, and you want to keep abutting loops
     * separate in the output rather than merging their regions together (e.g.
     * assembling loops for Kansas City, KS and Kansas City, MO simultaneously).
     */
    public boolean getXorEdges() {
      return xorEdges;
    }

    /**
     * Default value: false
     */
    public boolean getValidate() {
      return validate;
    }

    /**
     * Default value: 0
     */
    public S1Angle getMergeDistance() {
      return mergeDistance;
    }

    /**
     * If true, isValid() is called on all loops and polygons before
     * constructing them. If any loop is invalid (e.g. self-intersecting), it is
     * rejected and returned as a set of "unused edges". Any remaining valid
     * loops are kept. If the entire polygon is invalid (e.g. two loops
     * intersect), then all loops are rejected and returned as unused edges.
     */
    public void setValidate(boolean validate) {
      this.validate = validate;
    }

    /**
     * If set to a positive value, all vertices that are separated by at most
     * this distance will be merged together. In addition, vertices that are
     * closer than this distance to a non-incident edge will be spliced into it
     * (TODO).
     *
     *  The merging is done in such a way that all vertex-vertex and vertex-edge
     * distances in the output are greater than 'merge_distance'.
     *
     *  This method is useful for assembling polygons out of input data where
     * vertices and/or edges may not be perfectly aligned.
     */
    public void setMergeDistance(S1Angle mergeDistance) {
      this.mergeDistance = mergeDistance;
    }

    // Used for testing only
    void setUndirectedEdges(boolean undirectedEdges) {
      this.undirectedEdges = undirectedEdges;
    }

    // Used for testing only
    void setXorEdges(boolean xorEdges) {
      this.xorEdges = xorEdges;
    }
  }

  public Options options() {
    return options;
  }

  /**
   * Add the given edge to the polygon builder. This method should be used for
   * input data that may not follow S2 polygon conventions. Note that edges are
   * not allowed to cross each other. Also note that as a convenience, edges
   * where v0 == v1 are ignored.
   */
  public void addEdge(SpherePoint v0, SpherePoint v1) {
    // If xor_edges is true, we look for an existing edge in the opposite
    // direction. We either delete that edge or insert a new one.

    if (v0.equals(v1)) {
      return;
    }

    if (options.getXorEdges()) {
      Multiset<SpherePoint> candidates = edges.get(v1);
      if (candidates != null && candidates.count(v0) > 0) {
        eraseEdge(v1, v0);
        return;
      }
    }

    if (edges.get(v0) == null) {
      edges.put(v0, HashMultiset.<SpherePoint>create());
    }

    edges.get(v0).add(v1);
    if (options.getUndirectedEdges()) {
      if (edges.get(v1) == null) {
        edges.put(v1, HashMultiset.<SpherePoint>create());
      }
      edges.get(v1).add(v0);
    }
  }

  /**
   * Add all edges in the given loop. If the sign() of the loop is negative
   * (i.e. this loop represents a hole), the reverse edges are added instead.
   * This implies that "shells" are CCW and "holes" are CW, as required for the
   * directed edges convention described above.
   *
   * This method does not take ownership of the loop.
   */
  public void addLoop(SphereLoop loop) {
    int sign = loop.sign();
    for (int i = loop.numVertices(); i > 0; --i) {
      // Vertex indices need to be in the range [0, 2*num_vertices()-1].
      addEdge(loop.vertex(i), loop.vertex(i + sign));
    }
  }

  /**
   * Add all loops in the given polygon. Shells and holes are added with
   * opposite orientations as described for AddLoop(). This method does not take
   * ownership of the polygon.
   */
  public void addPolygon(SpherePolygon polygon) {
    for (int i = 0; i < polygon.numLoops(); ++i) {
      addLoop(polygon.loop(i));
    }
  }

  /**
   * Assembles the given edges into as many non-crossing loops as possible. When
   * there is a choice about how to assemble the loops, then CCW loops are
   * preferred. Returns true if all edges were assembled. If "unused_edges" is
   * not NULL, it is initialized to the set of edges that could not be assembled
   * into loops.
   *
   *  Note that if xor_edges() is false and duplicate edge pairs may be present,
   * then undirected_edges() should be specified unless all loops can be
   * assembled in a counter-clockwise direction. Otherwise this method may not
   * be able to assemble all loops due to its preference for CCW loops.
   *
   * This method resets the S2PolygonBuilder state so that it can be reused.
   */
  public boolean assembleLoops(List<SphereLoop> loops, List<SphereEdge> unusedEdges) {
    if (options.getMergeDistance().radians() > 0) {
//      mergeVertices();
    }

    List<SphereEdge> dummyUnusedEdges = Lists.newArrayList();
    if (unusedEdges == null) {
      unusedEdges = dummyUnusedEdges;
    }

    // We repeatedly choose an arbitrary edge and attempt to assemble a loop
    // starting from that edge. (This is always possible unless the input
    // includes extra edges that are not part of any loop.)

    unusedEdges.clear();
    while (!edges.isEmpty()) {
      Map.Entry<SpherePoint, Multiset<SpherePoint>> edge = edges.entrySet().iterator().next();

      SpherePoint v0 = edge.getKey();
      SpherePoint v1 = edge.getValue().iterator().next();

      SphereLoop loop = assembleLoop(v0, v1, unusedEdges);
      if (loop == null) {
        continue;
      }

      // In the case of undirected edges, we may have assembled a clockwise
      // loop while trying to assemble a CCW loop. To fix this, we assemble
      // a new loop starting with an arbitrary edge in the reverse direction.
      // This is guaranteed to assemble a loop that is interior to the previous
      // one and will therefore eventually terminate.

      while (options.getUndirectedEdges() && !loop.isNormalized()) {
        loop = assembleLoop(loop.vertex(1), loop.vertex(0), unusedEdges);
      }
      loops.add(loop);
      eraseLoop(loop, loop.numVertices());
    }
    return unusedEdges.isEmpty();
  }

  /**
   * Like AssembleLoops, but normalizes all the loops so that they enclose less
   * than half the sphere, and then assembles the loops into a polygon.
   *
   *  For this method to succeed, there should be no duplicate edges in the
   * input. If this is not known to be true, then the "xor_edges" option should
   * be set (which is true by default).
   *
   *  Note that S2Polygons cannot represent arbitrary regions on the sphere,
   * because of the limitation that no loop encloses more than half of the
   * sphere. For example, an S2Polygon cannot represent a 100km wide band around
   * the equator. In such cases, this method will return the *complement* of the
   * expected region. So for example if all the world's coastlines were
   * assembled, the output S2Polygon would represent the land area (irrespective
   * of the input edge or loop orientations).
   */
  public boolean assemblePolygon(SpherePolygon polygon, List<SphereEdge> unusedEdges) {
    List<SphereLoop> loops = Lists.newArrayList();
    boolean success = assembleLoops(loops, unusedEdges);

    // If edges are undirected, then all loops are already CCW. Otherwise we
    // need to make sure the loops are normalized.
    if (!options.getUndirectedEdges()) {
      for (int i = 0; i < loops.size(); ++i) {
        loops.get(i).normalize();
      }
    }
    if (options.getValidate() && !SpherePolygon.isValid(loops)) {
      if (unusedEdges != null) {
        for (SphereLoop loop : loops) {
          rejectLoop(loop, loop.numVertices(), unusedEdges);
        }
      }
      return false;
    }
    polygon.init(loops);
    return success;
  }

  /**
   * Convenience method for when you don't care about unused edges.
   */
  public SpherePolygon assemblePolygon() {
    SpherePolygon polygon = new SpherePolygon();
    List<SphereEdge> unusedEdges = Lists.newArrayList();

    assemblePolygon(polygon, unusedEdges);

    return polygon;
  }

  // Debugging functions:

  protected void dumpEdges(SpherePoint v0) {
    log.info(v0.toString());
    Multiset<SpherePoint> vset = edges.get(v0);
    if (vset != null) {
      for (SpherePoint v : vset) {
        log.info("    " + v.toString());
      }
    }
  }

  protected void dump() {
    for (SpherePoint v : edges.keySet()) {
      dumpEdges(v);
    }
  }

  private void eraseEdge(SpherePoint v0, SpherePoint v1) {
    // Note that there may be more than one copy of an edge if we are not XORing
    // them, so a VertexSet is a multiset.

    Multiset<SpherePoint> vset = edges.get(v0);
    // assert (vset.count(v1) > 0);
    vset.remove(v1);
    if (vset.isEmpty()) {
      edges.remove(v0);
    }

    if (options.getUndirectedEdges()) {
      vset = edges.get(v1);
      // assert (vset.count(v0) > 0);
      vset.remove(v0);
      if (vset.isEmpty()) {
        edges.remove(v1);
      }
    }
  }

  private void eraseLoop(List<SpherePoint> v, int n) {
    for (int i = n - 1, j = 0; j < n; i = j++) {
      eraseEdge(v.get(i), v.get(j));
    }
  }

  private void eraseLoop(SphereLoop v, int n) {
    for (int i = n - 1, j = 0; j < n; i = j++) {
      eraseEdge(v.vertex(i), v.vertex(j));
    }
  }

  /**
   * We start at the given edge and assemble a loop taking left turns whenever
   * possible. We stop the loop as soon as we encounter any vertex that we have
   * seen before *except* for the first vertex (v0). This ensures that only CCW
   * loops are constructed when possible.
   */
  private SphereLoop assembleLoop(SpherePoint v0, SpherePoint v1, List<SphereEdge> unusedEdges) {

    // The path so far.
    List<SpherePoint> path = Lists.newArrayList();

    // Maps a vertex to its index in "path".
    Map<SpherePoint, Integer> index = Maps.newHashMap();
    path.add(v0);
    path.add(v1);

    index.put(v1, 1);

    while (path.size() >= 2) {
      // Note that "v0" and "v1" become invalid if "path" is modified.
      v0 = path.get(path.size() - 2);
      v1 = path.get(path.size() - 1);

      SpherePoint v2 = null;
      boolean v2Found = false;
      Multiset<SpherePoint> vset = edges.get(v1);
      if (vset != null) {
        for (SpherePoint v : vset) {
          // We prefer the leftmost outgoing edge, ignoring any reverse edges.
          if (v.equals(v0)) {
            continue;
          }
          if (!v2Found || Sphere.orderedCCW(v0, v2, v, v1)) {
            v2 = v;
          }
          v2Found = true;
        }
      }
      if (!v2Found) {
        // We've hit a dead end. Remove this edge and backtrack.
        unusedEdges.add(new SphereEdge(v0, v1));
        eraseEdge(v0, v1);
        index.remove(v1);
        path.remove(path.size() - 1);
      } else if (index.get(v2) == null) {
        // This is the first time we've visited this vertex.
        index.put(v2, path.size());
        path.add(v2);
      } else {
        // We've completed a loop. Throw away any initial vertices that
        // are not part of the loop.
        path = path.subList(index.get(v2), path.size());

        if (options.getValidate() && !SphereLoop.isValid(path)) {
          // We've constructed a loop that crosses itself, which can only happen
          // if there is bad input data. Throw away the whole loop.
          rejectLoop(path, path.size(), unusedEdges);
          eraseLoop(path, path.size());
          return null;
        }
        return new SphereLoop(path);
      }
    }
    return null;
  }

  /** Erases all edges of the given loop and marks them as unused. */
  private void rejectLoop(SphereLoop v, int n, List<SphereEdge> unusedEdges) {
    for (int i = n - 1, j = 0; j < n; i = j++) {
      unusedEdges.add(new SphereEdge(v.vertex(i), v.vertex(j)));
    }
  }

  /** Erases all edges of the given loop and marks them as unused. */
  private void rejectLoop(List<SpherePoint> v, int n, List<SphereEdge> unusedEdges) {
    for (int i = n - 1, j = 0; j < n; i = j++) {
      unusedEdges.add(new SphereEdge(v.get(i), v.get(j)));
    }
  }

  /** Moves a set of vertices from old to new positions. */
  private void moveVertices(Map<SpherePoint, SpherePoint> mergeMap) {
    if (mergeMap.isEmpty()) {
      return;
    }

    // We need to copy the set of edges affected by the move, since
    // this.edges_could be reallocated when we start modifying it.
    List<SphereEdge> edgesCopy = Lists.newArrayList();
    for (Map.Entry<SpherePoint, Multiset<SpherePoint>> edge : this.edges.entrySet()) {
      SpherePoint v0 = edge.getKey();
      Multiset<SpherePoint> vset = edge.getValue();
      for (SpherePoint v1 : vset) {
        if (mergeMap.get(v0) != null || mergeMap.get(v1) != null) {

          // We only need to modify one copy of each undirected edge.
          if (!options.getUndirectedEdges() || v0.lessThan(v1)) {
            edgesCopy.add(new SphereEdge(v0, v1));
          }
        }
      }
    }

    // Now erase all the old edges, and add all the new edges. This will
    // automatically take care of any XORing that needs to be done, because
    // EraseEdge also erases the sibiling of undirected edges.
    for (int i = 0; i < edgesCopy.size(); ++i) {
      SpherePoint v0 = edgesCopy.get(i).getStart();
      SpherePoint v1 = edgesCopy.get(i).getEnd();
      eraseEdge(v0, v1);
      if (mergeMap.get(v0) != null) {
        v0 = mergeMap.get(v0);
      }
      if (mergeMap.get(v1) != null) {
        v1 = mergeMap.get(v1);
      }
      addEdge(v0, v1);
    }
  }

  /**
   * An S2Point that can be marked. Used in PointIndex.
   */
  private class MarkedS2Point {
    private SpherePoint point;
    private boolean mark;

    public MarkedS2Point(SpherePoint point) {
      this.point = point;
      this.mark = false;
    }

    public boolean isMarked() {
      return mark;
    }

    public SpherePoint getPoint() {
      return point;
    }

    public void mark() {
      // assert (!isMarked());
      this.mark = true;
    }
  }
}
