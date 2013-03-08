package com.caseystella.tsq.hadoop;

public class DataPreprocessor 
{
	private long numDataPoints;
	private int numTasks;
	private long numDataPointsPerTaskSansOverlap;
	private long overlap;
	
	public DataPreprocessor( long numDataPoints
						   , int numTasks
						   , float overlapPercentage
						   )
	{
		this(numDataPoints,numDataPoints / numTasks, overlapPercentage);
	}
	
	public DataPreprocessor( long numDataPoints
						   , long numDataPointsPerTask
						   , float overlapPercentage
						   )
	{
		this(numDataPoints, numDataPointsPerTask, (int)(overlapPercentage*numDataPointsPerTask) );
	}
	
	public DataPreprocessor( long numDataPoints
						   , int numTasks
						   , long overlap
						   )
	{
		this(numDataPoints,numDataPoints / numTasks, overlap);
		//don't trust the other constructor
		this.numTasks = numTasks;
	}
	
	public DataPreprocessor( long numDataPoints
						   , long numDataPointsPerTask
						   , long overlap
						   )
	{
		this.numDataPoints = numDataPoints;
		this.numDataPointsPerTaskSansOverlap = numDataPointsPerTask;
		//TODO: check overflow here; dangerous cast
		this.numTasks = (int)(numDataPoints / numDataPointsPerTask);
		this.overlap = overlap;
	}
	
	public long getNumDataPoints() {
		return numDataPoints;
	}
	public long getNumDataPointsPerTask() {
		return numDataPointsPerTaskSansOverlap + overlap;
	}
	public int getNumTasks() {
		return numTasks;
	}
	public long getOverlap() {
		return overlap;
	}
	
	public long getBlockSize()
	{
		return Constants.BLOCK_HEADER_SIZE + getNumDataPointsPerTask()*Constants.KV_SIZE;
	}
}
