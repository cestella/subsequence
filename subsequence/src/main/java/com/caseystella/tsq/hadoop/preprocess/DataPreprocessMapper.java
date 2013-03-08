package com.caseystella.tsq.hadoop.preprocess;

import java.io.IOException;

import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.Writable;
import org.apache.hadoop.mapred.MapReduceBase;
import org.apache.hadoop.mapred.Mapper;
import org.apache.hadoop.mapred.OutputCollector;
import org.apache.hadoop.mapred.Reporter;



public class DataPreprocessMapper extends MapReduceBase implements Mapper<Writable, DoubleWritable, PartitionKey, DoubleWritable>
{

	@Override
	public void map(Writable key, DoubleWritable value,
			OutputCollector<PartitionKey, DoubleWritable> output,
			Reporter reporter) throws IOException 
			
	{
		//reporter.g
		// TODO Auto-generated method stub
		
	}

}
