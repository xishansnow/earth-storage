/*
 * Copyright 2011 Google Inc.
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
 * 球面有向边类
 * An abstract directed edge from one S2Point to another S2Point.
 *
 * @author kirilll@google.com (Kirill Levin)
 */
public final class SphereEdge {

  private final SpherePoint start;
  private final SpherePoint end;

  public SphereEdge(SpherePoint start, SpherePoint end) {
    this.start = start;
    this.end = end;
  }

  public SpherePoint getStart() {
    return start;
  }

  public SpherePoint getEnd() {
    return end;
  }

  @Override
  public String toString() {
    return String.format("Edge: (%s -> %s)\n   or [%s -> %s]",
        start.toDegreesString(), end.toDegreesString(), start, end);
  }

  @Override
  public int hashCode() {
    return getStart().hashCode() - getEnd().hashCode();
  }

  @Override
  public boolean equals(Object o) {
    if (o == null || !(o instanceof SphereEdge)) {
      return false;
    }
    SphereEdge other = (SphereEdge) o;
    return getStart().equals(other.getStart()) && getEnd().equals(other.getEnd());
  }
}
