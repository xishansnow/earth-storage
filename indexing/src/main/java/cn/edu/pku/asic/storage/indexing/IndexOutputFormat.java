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
package cn.edu.pku.asic.storage.indexing;

import edu.ucr.cs.bdlab.beast.geolite.IFeature;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.OutputCommitter;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

import java.io.IOException;

/**
 * An output format that writes pairs that represent a partition ID and a feature. All records in the same partition ID
 * are written to a separate file using the configured output format. Additionally, a master file is added that
 * stores metadata about partitions including partition ID, MBR, file name, number of records, and total size in bytes.
 * @author Ahmed Eldawy
 *
 */
public class IndexOutputFormat extends FileOutputFormat<Integer, IFeature> {

  @Override
  public RecordWriter<Integer, IFeature> getRecordWriter(TaskAttemptContext task) throws IOException {
    Path file = getDefaultWorkFile(task, "").getParent();
    return new edu.ucr.cs.bdlab.beast.indexing.IndexRecordWriter(task, file);
  }

  @Override
  public synchronized OutputCommitter getOutputCommitter(TaskAttemptContext task) throws IOException {
    Path jobOutputPath = getOutputPath(task);
    return new edu.ucr.cs.bdlab.beast.indexing.IndexMasterFileCommitter(jobOutputPath, task);
  }
}
