package com.caseystella.tsq;

/**
 * Stateful circular array for use in finding LB_Keogh envelopes
 * @author cstella
 *
 */
public class CircularArray 
{
	private int size;
	private int capacity;
	private int frontIdx;
	private int rearIdx;
	private int[] data;
	public CircularArray(int capacity)
	{
		this.capacity = capacity;
		this.size = 0;
		this.data = new int[capacity];
		this.frontIdx = 0;
		this.rearIdx = capacity - 1;
	}
	
	public void pushBack(int value)
	{
		data[rearIdx--] = value;
		if(rearIdx < 0)
		{
			rearIdx = capacity - 1;
		}
		size++;
	}
	public void popFront()
	{
		frontIdx--;
		if(frontIdx < 0)
		{
			frontIdx = capacity - 1;
		}
		size--;
	}
	public void popBack()
	{
		rearIdx = (rearIdx + 1) % capacity;
		size--;
	}
	public int front()
	{
		int aux = frontIdx - 1;
		if(aux < 0)
		{
			aux = capacity - 1;
		}
		return data[aux];
	}
	public int back()
	{
		int aux = (rearIdx + 1) % capacity;
		return data[aux];
	}
	public boolean isEmpty()
	{
		return size == 0;
	}
}
