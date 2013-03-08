package org.apache.hadoop.io;

import java.io.IOException;
import java.io.OutputStream;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.io.SequenceFile.Metadata;


public class FakeWriter extends SequenceFile.Writer 
{
	
	@SuppressWarnings("deprecation")
	public FakeWriter(Configuration conf, OutputStream o) throws IOException 
	{
		super(conf, new FSDataOutputStream(o), NullWritable.class, DoubleWritable.class, new Metadata());
		
	}
}
