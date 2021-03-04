/*
 * Copyright 2018 University of California, Riverside
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
package cn.edu.pku.asic.earthstorage.common.io;

import org.locationtech.jts.geom.Envelope;

import java.util.function.BiFunction;

/**
 * A function that parses a CSV representation of envelopes represented as (minimum coordinates, maximum coordinates)
 * Each of them is a coordinate of a point as (x, y, z, ...).
 */
public class CSVEnvelopeDecoderJTS implements BiFunction<String, Envelope, Envelope> {

  private CSVEnvelopeDecoderJTS() {}

  public static final CSVEnvelopeDecoderJTS instance = new CSVEnvelopeDecoderJTS();

  @Override
  public Envelope apply(String s, Envelope env) {
    String[] parts = s.split(",");
    assert (parts.length % 2) == 0;
    if (env == null)
      env = new Envelope();
    double minx = Double.parseDouble(parts[0]);
    double miny = Double.parseDouble(parts[1]);
    double maxx = Double.parseDouble(parts[2]);
    double maxy = Double.parseDouble(parts[3]);
    env.init(minx, maxx, miny, maxy);
    return env;
  }
}
