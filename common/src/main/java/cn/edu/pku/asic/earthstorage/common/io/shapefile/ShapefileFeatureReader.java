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

import cn.edu.pku.asic.earthstorage.common.geolite.Feature;
import cn.edu.pku.asic.earthstorage.common.geolite.IFeature;
import cn.edu.pku.asic.earthstorage.common.io.FeatureReader;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.LocalFileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.apache.spark.beast.CRSServer;
import org.geotools.referencing.CRS;
import org.locationtech.jts.geom.Geometry;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import java.io.*;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 * Reads shapefiles that contain geometries in .shp files and data in .dbf files.
 */
@FeatureReader.Metadata(
    description = "Reads Esri-formatted shapefiles. The input can be a .shp file, a .zip file that contains" +
        ".shp and .dbf files or a directory that contains a mix of these two",
    extension = ".shp",
    shortName = "shapefile",
    filter = "*.zip\n*.shp",
    noSplit = true
)
public class ShapefileFeatureReader extends FeatureReader {
  /**Logger for this class*/
  private static final Log LOG = LogFactory.getLog(ShapefileFeatureReader.class);

  /**Used if the input is a ZIP file*/
  protected ZipFile zipFile;

  /**
   * Only if the input is a ZIP file, this list contains the names of the .shp files (including the extension) that
   * are contained within the ZIP file and have not been read yet. Once a file is opened, it is removed from this list.
   * This means that if the ZIP file contains one shapefile, this list will be empty as soon as the first record is
   * retrieved.
   */
  protected List<String> containedShapefiles;

  /**An underlying reader for the shape file*/
  protected ShapefileGeometryReader shapefileReader;

  /**An underlying reader for the DBF file*/
  protected DBFReader dbfReader;

  /**Configuration of this reader*/
  protected Configuration conf;

  /**The local temporary ZIP file if the input is a remote ZIP file*/
  protected File tempZipFile;

  /**The coordinate reference system defined in the shapefile*/
  protected CoordinateReferenceSystem crs;

  /**Path to the .prj file if needed*/
  protected String prjFileName;

  /**SRID to set in created geometries*/
  protected int srid;

  /**The current value pointed by the reader*/
  protected IFeature currentValue;

  public ShapefileFeatureReader() {
    this.shapefileReader = new ShapefileGeometryReader();
    this.dbfReader = new DBFReader();
  }

  @Override
  public void initialize(InputSplit inputSplit, Configuration conf) throws IOException {
    FileSplit fsplit = (FileSplit) inputSplit;
    Path shapefilePath = fsplit.getPath();
    this.initialize(shapefilePath, conf);
  }

  public void initialize(Path shapefilePath, Configuration conf) throws IOException {
    this.conf = conf;
    String shapeFileName = shapefilePath.getName();
    int lastDot = shapeFileName.lastIndexOf('.');
    String extension = shapeFileName.substring(lastDot).toLowerCase();
    if (extension.equals(".shp")) {
      // The input split points to the shapefile directly, assume the DBF file is in the same path
      Path dbfFilePath = new Path(shapefilePath.getParent(), shapeFileName.substring(0, lastDot) + ".dbf");
      this.prjFileName = new Path(shapefilePath.getParent(), shapeFileName.substring(0, lastDot) + ".prj").toString();
      this.shapefileReader.initialize(shapefilePath, conf);
      this.dbfReader.initialize(dbfFilePath, conf);
    } else if (extension.equals(".zip")) {
      // A compressed zip file, find the .shp and .dbf entries inside the file
      FileSystem fs = shapefilePath.getFileSystem(conf);
      if (! (fs instanceof LocalFileSystem)) {
        // The ZIP file stored in a remote file system, copy the file locally to use the java.utils.zip library
        tempZipFile = File.createTempFile(shapeFileName, "temp");
        fs.copyToLocalFile(shapefilePath, new Path(tempZipFile.toString()));
        zipFile = new ZipFile(tempZipFile);
      } else {
        // The ZIP file is stored locally, open it directly
        zipFile = new ZipFile(shapefilePath.toUri().getPath());
      }
      // Find the .dbf and .shp entries in the ZIP file
      containedShapefiles = new ArrayList<>();
      Enumeration<? extends ZipEntry> entries = zipFile.entries();
      while (entries.hasMoreElements()) {
        ZipEntry entry = entries.nextElement();
        if (entry.getName().toLowerCase().endsWith(".shp"))
          containedShapefiles.add(entry.getName());
      }
      if (containedShapefiles.isEmpty())
        LOG.warn(String.format("No .shp files found in the zip file '%s'", shapefilePath));
      openNextShapefile();
    } else {
      throw new RuntimeException(String.format("Cannot open the shapefile '%s'. Unrecognized extension.", shapeFileName));
    }
    // Retrieve the spatial reference identifier and assume the default CRS to be EPSG:4326
    this.srid = 4326;
    CoordinateReferenceSystem crs = getCRS();
    if (crs != null) {
      int srid = CRSServer.crsToSRID(crs, CRSServer.sparkConfFromHadoopConf(conf));
      if (srid != 0)
        this.srid = srid;
    }
    this.shapefileReader.setSRID(this.srid);
  }

  @Override
  public CoordinateReferenceSystem getCRS() {
    if (this.crs == null) {
      FileSystem fileSystem = null;
      InputStream in = null;
      try {
        // Read and parse the associated .prj file
        if (zipFile == null) {
          // Read a local file
          Path prjFilePath = new Path(this.prjFileName);
          fileSystem = prjFilePath.getFileSystem(conf);
          if (!fileSystem.exists(prjFilePath))
            return null;
          in = fileSystem.open(prjFilePath);
        } else {
          // Read from ZIP file
          ZipEntry prjZipEntry = zipFile.getEntry(prjFileName);
          if (prjZipEntry == null)
            return null;
          in = zipFile.getInputStream(prjZipEntry);
        }
        String wktCRS = "";
        byte[] buffer = new byte[1024];
        int length = 0;
        while ((length = in.read(buffer, 0, buffer.length)) > 0) {
          wktCRS += new String(buffer, 0, length);
        }
        this.crs = CRS.parseWKT(wktCRS);
      } catch (IOException e) {
        this.crs = super.getCRS();
      } catch (FactoryException e) {
        e.printStackTrace();
      } finally {
        if (fileSystem != null) {
          if (in != null) {
            try {
              in.close();
            } catch (IOException e) {
              e.printStackTrace();
            }
          }
        }
      }
    }
    return this.crs;
  }

  /**
   * Opens the last shapefile in the list {@link #containedShapefiles} and removes the opened file from the list.
   * If the file is opened correctly, a {@code true} is returned. If the list is empty or the file could not be opened,
   * a {@code false is returned}
   * @return {@code true} if a new file was opened.
   * @throws IOException if an error happens while opening the shapefile
   */
  protected boolean openNextShapefile() throws IOException {
    if (containedShapefiles == null || containedShapefiles.isEmpty())
      return false;
    // TODO handle the case if the DBF, SHX, and SHP files are not all lower case
    // Locate the SHP file in the ZIP file
    String shpFileName = containedShapefiles.remove(containedShapefiles.size() - 1);
    ZipEntry shpEntry = zipFile.getEntry(shpFileName);
    if (shpEntry == null) {
      LOG.warn(String.format("Could not retrieve the shapefile '%s' from the ZIP file '%s'", shpFileName, zipFile.getName()));
      return false;
    }
    // Locate the DBF file in the ZIP file
    String dbfFileName = shpFileName.replace(".shp", ".dbf");
    ZipEntry dbfEntry = zipFile.getEntry(dbfFileName);
    if (dbfEntry == null) {
      LOG.warn(String.format("Could not retrieve the file '%s' from the ZIP file '%s'", dbfFileName, zipFile.getName()));
      return false;
    }
    // Locate the SHX file in the ZIP file
    String shxFileName = shpFileName.replace(".shp", ".shx");
    ZipEntry shxEntry = zipFile.getEntry(shxFileName);
    if (shxEntry == null) {
      LOG.warn(String.format("Could not retrieve the file '%s' from the ZIP file '%s'", shxFileName, zipFile.getName()));
      // Continue with a .shx file by assuming consecutive features in the shape file
    }
    // Locate the PRJ file but do not open it yet
    prjFileName = shpFileName.replace(".shp", ".prj");
    // 1- Load all shape indexes
    DataInputStream shxIn = new DataInputStream(new BufferedInputStream(zipFile.getInputStream(shxEntry)));
    shapefileReader.readIndexFile(shxIn);
    shxIn.close();
    // 2- Initialize the shape file
    DataInputStream shpIn = new DataInputStream(new BufferedInputStream(zipFile.getInputStream(shpEntry)));
    shapefileReader.initialize(shpIn, conf);
    // 3- Initialize the DBF file
    DataInputStream dbfIn = new DataInputStream(new BufferedInputStream(zipFile.getInputStream(dbfEntry)));
    dbfReader.initialize(dbfIn, conf);
    return true;
  }

  @Override
  public boolean nextKeyValue() throws IOException {
    if (!shapefileReader.nextKeyValue()) {
      // Current file terminated. Try to open the next file
      shapefileReader.close();
      dbfReader.close();
      if (!openNextShapefile())
        return false;
      shapefileReader.nextKeyValue();
    }
    // The shapefileReader might skip some records based on the filter MBR.
    // Make sure we skip similar records in the DBF file
    while (dbfReader.iRecord < shapefileReader.iShape) {
      if (!dbfReader.nextKeyValue())
        throw new RuntimeException(String.format("Could not reach record #%d in the DBF file", shapefileReader.iShape));
    }
    Geometry geometry = shapefileReader.getCurrentValue();
    this.currentValue = Feature.create(dbfReader.getCurrentValue(), geometry);
    return true;
  }

  @Override
  public IFeature getCurrentValue() throws IOException {
    return currentValue;
  }

  @Override
  public float getProgress() throws IOException {
    return dbfReader.getProgress();
  }

  @Override
  public void close() throws IOException {
    shapefileReader.close();
    dbfReader.close();
    if (zipFile != null) {
      zipFile.close();
      zipFile = null;
    }
    if (tempZipFile != null) {
      tempZipFile.delete();
      tempZipFile = null;
    }
  }
}
