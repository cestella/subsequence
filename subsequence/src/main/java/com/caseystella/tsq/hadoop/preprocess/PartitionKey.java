package com.caseystella.tsq.hadoop.preprocess;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.WritableComparable;

public class PartitionKey implements WritableComparable<PartitionKey>
{
	private LongWritable partitionId;
	private LongWritable offset;
	
	public PartitionKey()
	{
		partitionId = new LongWritable();
		offset = new LongWritable();
	}
	
	public PartitionKey(LongWritable partitionId, LongWritable offset)
	{
		this.partitionId = partitionId;
		this.offset = offset;
	}
	
	@Override
	public void write(DataOutput out) throws IOException {
		out.writeLong(partitionId.get());
		out.writeLong(offset.get());
		
	}

	@Override
	public void readFields(DataInput in) throws IOException {
		partitionId = new LongWritable(in.readLong());
		offset = new LongWritable(in.readLong());
		
	}

	@Override
	public int compareTo(PartitionKey key) 
	{
		if(partitionId.get() == key.partitionId.get())
		{
			if(offset.get() == key.offset.get())
			{
				return 0;
			}
			else if(offset.get() < key.offset.get())
			{
				return -1;
			}
			else
			{
				return 1;
			}
			
		}
		else if(partitionId.get() < key.partitionId.get())
		{
			return -1;
		}
		else 
		{
			return 1;
		}
	}

}
