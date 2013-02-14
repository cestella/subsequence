package com.caseystella.tsq;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Random;

import junit.framework.Assert;

import org.junit.Test;

import com.caseystella.tsq.dtw.DTWQueryEngine;


public class DTWQueryTest 
{
	public static int DATA_SIZE = 1000000;
	public static int QUERY_SIZE = 100;
	
	public static double PI = 3.14159;
	private byte[] createSyntheticData(int size, int offset)
	{
		double delta = 2*PI / DATA_SIZE;
		byte[] ret = new byte[8*DATA_SIZE];
		int i = 0;
		Random r = new Random(0);
		for( double pt = 0;pt < 2*PI;pt += delta)
		{
			double multiplier = (i/8 < offset) ? 0:7*r.nextGaussian();
			double y = multiplier*Math.sin(pt);
			byte[] bytes = toByteArray(y);
			for(int k = 0;k < 8;++k)
			{
				ret[i++] = bytes[k];
			}
		}
		return ret;
	}
	public static byte[] toByteArray(double value) {
	    byte[] bytes = new byte[8];
	    ByteBuffer.wrap(bytes).putDouble(value);
	    return bytes;
	}

	public static double toDouble(byte[] bytes) {
	    return ByteBuffer.wrap(bytes).getDouble();
	}
	
	@Test
	public void test() throws IOException 
	{
		int offset = 4000;
		byte[] data = createSyntheticData(DATA_SIZE, 0);
		double[] dataDouble = new double[DATA_SIZE];
		double[] query = new double[QUERY_SIZE];
		{
			DataInputStream stream = new DataInputStream(new ByteArrayInputStream(data));
			for(int i = 0;i < DATA_SIZE;++i)
			{
				dataDouble[i] = stream.readDouble();
			}
			
			
			for(int i = 0;i < query.length;++i)
			{
				query[i] = dataDouble[i + offset];
			}
		}
		DTWQueryEngine engine = new DTWQueryEngine();
		Location loc = engine.query( query
								   , new DataInputStream(new ByteArrayInputStream(data))
								   , 0.008
								   );
		System.out.println("Loc = " + loc.offset + ", dist = " + loc.distance + ", data scanned = " + loc.dataScanned);
		Assert.assertEquals(offset, loc.offset);
	}

}
