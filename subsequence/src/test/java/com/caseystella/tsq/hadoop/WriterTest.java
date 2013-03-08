package com.caseystella.tsq.hadoop;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import junit.framework.Assert;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.DoubleWritable;
import org.apache.hadoop.io.FakeWriter;
import org.apache.hadoop.io.NullWritable;
import org.junit.Test;

public class WriterTest 
{
	int getSizeOnDisk(int k) throws IOException
	{
		Configuration conf = new Configuration();
		conf.setInt("fs.local.block.size", 32 * 1024 * 1024);
		{
			ByteArrayOutputStream os = new ByteArrayOutputStream();
			FakeWriter writer = new FakeWriter(conf, os);
			for(int i = 0;i < k;++i)
			{
				writer.append(NullWritable.get(), new DoubleWritable(0));
			}
			writer.syncFs();
			writer.close();
			return os.toByteArray().length;
		}
	}
	
	@Test
	public void test() throws IOException 
	{
		
		int sizeOfOne = getSizeOnDisk(1);
		int sizeOfTwo = getSizeOnDisk(2);
		int sizeOfThree = getSizeOnDisk(3);
		
		System.out.println("Size of one = " + sizeOfOne);

		System.out.println("Size of two = " + sizeOfTwo + ", diff = " + (sizeOfTwo - sizeOfOne));

		System.out.println("Size of three = " + sizeOfThree+ ", diff = " + (sizeOfThree - sizeOfTwo));
		
		Assert.assertEquals("Ensure size of key is internally consistent", sizeOfTwo - sizeOfOne, sizeOfThree - sizeOfTwo);
		Assert.assertEquals("Ensure size is what we presume", sizeOfTwo - sizeOfOne, Constants.KV_SIZE);
		Assert.assertEquals("Ensure size of header is what we presume", sizeOfOne - Constants.KV_SIZE, Constants.BLOCK_HEADER_SIZE);
	}

}
