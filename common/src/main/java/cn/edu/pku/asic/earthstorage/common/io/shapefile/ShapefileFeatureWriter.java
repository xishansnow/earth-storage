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
package cn.edu.pku.asic.earthstorage.common.io.shapefile;

import cn.edu.pku.asic.earthstorage.common.geolite.IFeature;
import cn.edu.pku.asic.earthstorage.common.io.FeatureWriter;
import cn.edu.pku.asic.earthstorage.common.utils.FileUtil;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.spark.beast.CRSServer;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;

/**
 * Writes {@link IFeature} values in a shapefile format which consists of mainly three files, .shp, .shx, and .dbf files.
 */
@FeatureWriter.Metadata(extension = ".shp", shortName = "shapefile")
public class ShapefileFeatureWriter extends FeatureWriter {

  /**Writes the geometry part of the output*/
  protected ShapefileGeometryWriter shpWriter;

  /**Writes the database (DBF) file*/
  protected DBFWriter dbfWriter;

  /**The output shape path*/
  protected Path shpPath;

  /**A flag that indicates that the projection file has been written to avoid writing it twice.*/
  protected boolean prjFileWritten;

  /**System configuration. Used with the CRSServer to the CRS associated with the SRID*/
  protected Configuration conf;

  @Override
  public void initialize(Path shpPath, Configuration conf) throws IOException {
    this.shpPath = shpPath;
    this.conf = conf;
    shpWriter = new ShapefileGeometryWriter();
    shpWriter.initialize(shpPath, conf);
    Path dbfPath = new Path(shpPath.getParent(), FileUtil.replaceExtension(shpPath.getName(), ".dbf"));
    dbfWriter = new DBFWriter();
    dbfWriter.initialize(dbfPath, conf);
    this.prjFileWritten = false;
  }

  @Override
  public void initialize(OutputStream out, Configuration conf) {
  }

  @Override
  public void write(Object dummy, IFeature iFeature) throws IOException {
    if (!prjFileWritten) {
      if (shpPath != null) {
        CoordinateReferenceSystem crs = CRSServer.sridToCRS(iFeature.getGeometry().getSRID(),
            CRSServer.sparkConfFromHadoopConf(conf));
        if (crs != null) {
          Path prjPath = new Path(shpPath.getParent(), FileUtil.replaceExtension(shpPath.getName(), ".prj"));
          FileSystem fileSystem = prjPath.getFileSystem(conf);
          PrintWriter writer = new PrintWriter(fileSystem.create(prjPath));
          writer.print(crs.toWKT());
          writer.close();
        }
      }
      prjFileWritten = true;
    }
    shpWriter.write(dummy, iFeature.getGeometry());
    dbfWriter.write(dummy, iFeature);
  }

  @Override
  public int estimateSize(IFeature f) {
    // Binary storage size + MBR + record id and length + record in index file
    return f.getStorageSize() + 4 * 8 + 8 + 8;
  }

  @Override
  public void close(TaskAttemptContext taskAttemptContext) throws IOException {
    if (shpWriter != null)
      shpWriter.close(taskAttemptContext);
    if (dbfWriter != null)
      dbfWriter.close(taskAttemptContext);
  }

  public long getShapefileSize() {
    return shpWriter.getCurrentSize();
  }
}
