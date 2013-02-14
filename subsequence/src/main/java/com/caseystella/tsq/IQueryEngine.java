package com.caseystella.tsq;

import java.io.DataInputStream;

public interface IQueryEngine 
{
	
	public Location query(double[] q, DataInputStream data, double warpingWindow);
}
