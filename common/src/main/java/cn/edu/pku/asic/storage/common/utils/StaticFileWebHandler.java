/*
 * Copyright 2021 University of California, Riverside
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
package cn.edu.pku.asic.storage.common.utils;

import cn.edu.pku.asic.storage.common.cli.BeastOptions;
import cn.edu.pku.asic.storage.common.cli.WebMethod;
import org.apache.spark.SparkContext;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;

/**
 * A web handler for static files that works like a simple web server.
 * Enable this web handler with caution as it might expose unwanted files to the public.
 * For security reasons, this web handler is disabled by default, to enable, add the
 * parameter "-enableStaticFileHandling" when starting the server.
 */
public class StaticFileWebHandler extends cn.edu.pku.asic.storage.common.utils.AbstractWebHandler {

  private boolean enabled;

  @cn.edu.pku.asic.storage.common.utils.OperationParam(
      description = "Enable static file handling",
      defaultValue = "false"
  )
  public static final String EnableStaticFileHandling = "enableStaticFileHandling";

  @Override
  public void setup(SparkContext sc, BeastOptions opts) {
    super.setup(sc, opts);
    this.enabled = opts.getBoolean(EnableStaticFileHandling, false);
  }

  @WebMethod(url = "/(.*)", order = Integer.MAX_VALUE)
  public boolean handleStaticResource(String target, HttpServletRequest request, HttpServletResponse response) throws IOException {
    if (enabled)
      return super.handleStaticResource(target, request, response);
    else
      return false;
  }
}
