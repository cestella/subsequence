package com.caseystella.tsq.dtw;

import java.io.DataInputStream;
import java.io.IOException;
import java.nio.DoubleBuffer;
import java.util.Arrays;

import com.caseystella.tsq.CircularArray;
import com.caseystella.tsq.IQueryEngine;
import com.caseystella.tsq.Location;

public class DTWQueryEngine implements IQueryEngine
{
	public static double INF = 1e20;
	/// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
	public static int EPOCH = 100000;
	private static class Index implements Comparable<Index>
	{
		public double value;
		public int index;
		@Override
		public int compareTo(Index y) {
			return(int)(Math.abs(y.value) - Math.abs(value));
		}
	}
	public static class Bounds
	{
		final double[] lower;
		final double[] upper;
		final int querySize;
		public Bounds(double[] lower, double[] upper)
		{
			this.lower = lower;
			this.upper = upper;
			querySize = lower.length;
		}
		public Bounds(int querySize)
		{
			lower = new double[querySize];
			upper = new double[querySize];
			this.querySize = querySize;
		}
		public double[] getLower()
		{
			return lower;
		}
		public double[] getUpper()
		{
			return upper;
		}
		public int getQuerySize()
		{
			return querySize;
		}
	}
	
	/**
	 * Finding the envelop of min and max value for LB_Keogh
	 * Implementation idea is intoruduced by Danial Lemire in his paper
	 * "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
	 * @param query
	 * @param warpingWindow
	 * @return
	 */
	private Bounds computeLemireBounds(double[] buff, int len, int warpingWindow, Bounds bounds)
	{
		CircularArray upperDeque = new CircularArray(2*warpingWindow + 2);
		CircularArray  lowerDeque = new CircularArray(2*warpingWindow + 2);
		double[] lower = bounds.getLower();
		double[] upper = bounds.getUpper();
		lowerDeque.pushBack(0);
		upperDeque.pushBack(0);
		int doubleWindow = 2*warpingWindow + 1;
		for(int i = 1;i < len;++i)
		{
			if(i > warpingWindow)
			{
				int idx = i - warpingWindow - 1;
				upper[idx] = buff[lowerDeque.front()];
				lower[idx] = buff[upperDeque.front()];
				
				double data_i = buff[i];
				if(data_i > buff[i - 1])
				{
					upperDeque.popBack();
					while(!upperDeque.isEmpty() && data_i > buff[upperDeque.back()])
					{
						upperDeque.popBack();
					}
				}
				else
				{
					lowerDeque.popBack();
					while(!lowerDeque.isEmpty() && data_i < buff[lowerDeque.back()])
					{
						lowerDeque.popBack();
					}
				}
				
				upperDeque.pushBack(i);
				lowerDeque.pushBack(i);
				
				
				
				if(i == doubleWindow + upperDeque.front())
				{
					upperDeque.popFront();
				}
				else if(i == doubleWindow + lowerDeque.front())
				{
					lowerDeque.popFront();
				}
		        
		    }
		}
		for(int i = bounds.getQuerySize(); i < bounds.getQuerySize() + warpingWindow + 1;++i)
		{
			int idx = i - warpingWindow - 1;
			int upperFront = upperDeque.front();
			int lowerFront = lowerDeque.front();
			upper[idx] = buff[upperFront];
			lower[idx] = buff[lowerFront];
			if(i - upperFront >= doubleWindow)
			{
				upperDeque.popFront();
			}
			if(i - lowerFront >= doubleWindow)
			{
				lowerDeque.popFront();
			}
		}
		return bounds;
	}
	
	
	
	
	/** Calculate quick lower bound
	 * Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
	 * However, because of z-normalization the top and bottom cannot give siginifant benefits.
	 * And using the first and last points can be computed in constant time.
	 * The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
	*/
	private double lowerBoundKimHierarchy(double[] t
										 , double[] q
										 , int j //start position of the data in the circular buffer
										 , int len
										 , double mean
										 , double std
										 , double bsf
										 )
	{
		   /// 1 point at front and back
	    double d, lb; 
	    double x0 = (t[j] - mean) / std;
	    double y0 = (t[(len-1+j)] - mean) / std;
	    lb = dist(x0,q[0]) + dist(y0,q[len-1]);
	    if (lb >= bsf)   return lb; 

	    /// 2 points at front
	    double x1 = (t[(j+1)] - mean) / std;
	    d = Math.min(dist(x1,q[0]), dist(x0,q[1]));
	    d = Math.min(d, dist(x1,q[1]));
	    lb += d;
	    if (lb >= bsf)   return lb; 

	    /// 2 points at back
	    double y1 = (t[(len-2+j)] - mean) / std;
	    d = Math.min(dist(y1,q[len-1]), dist(y0, q[len-2]) );
	    d = Math.min(d, dist(y1,q[len-2]));
	    lb += d;
	    if (lb >= bsf)   return lb; 

	    /// 3 points at front
	    double x2 = (t[(j+2)] - mean) / std;
	    d = Math.min(dist(x0,q[2]), dist(x1, q[2]));
	    d = Math.min(d, dist(x2,q[2]));
	    d = Math.min(d, dist(x2,q[1]));
	    d = Math.min(d, dist(x2,q[0]));
	    lb += d;
	    if (lb >= bsf)   return lb;

	    /// 3 points at back
	    double y2 = (t[(len-3+j)] - mean) / std;
	    d = Math.min(dist(y0,q[len-3]), dist(y1, q[len-3]));
	    d = Math.min(d, dist(y2,q[len-3]));
	    d = Math.min(d, dist(y2,q[len-2]));
	    d = Math.min(d, dist(y2,q[len-1]));
	    lb += d;

	    return lb;
	}
	
	private static double dist(double x, double y)
	{
		double sub = x - y;
		return sub*sub;
	}

	

	/// LB_Keogh 1: Create Envelop for the query
	/// Note that because the query is known, envelop can be created once at the begenining.
	///
	/// Variable Explanation,
	/// order : sorted indices for the query.
	/// uo, lo: upper and lower envelops for the query, which already sorted.
	/// t     : a circular array keeping the current data.
	/// j     : index of the starting location in t
	/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
	double lb_keogh_cumulative( int[] order
							  , double[] t
							  , double[] uo
							  , double[] lo
							  , double[] cb
							  , int j
							  , int len
							  , double mean
							  , double std
							  , double best_so_far
							  )
	{
	    double lb = 0;
	    double x, d;

	    for (int i = 0; i < len && lb < best_so_far; i++)
	    {
	        x = (t[(order[i]+j)] - mean) / std;
	        d = 0;
	        if (x > uo[i])
	            d = dist(x,uo[i]);
	        else if(x < lo[i])
	            d = dist(x,lo[i]);
	        lb += d;
	        cb[order[i]] = d;
	    }
	    return lb;
	}

	/// LB_Keogh 2: Create Envelop for the data
	/// Note that the envelops have been created (in main function) when each data point has been read.
	///
	/// Variable Explanation,
	/// tz: Z-normalized data
	/// qo: sorted query
	/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
	/// l,u: lower and upper envelop of the current data
	double lb_keogh_data_cumulative( int[] order
								   , double[] tz
								   , double[] qo
								   , double[] cb
								   , DoubleBuffer l
								   , DoubleBuffer u
								   , int len
								   , double mean
								   , double std
								   , double best_so_far
								   )
	{
	    double lb = 0;
	    double uu,ll,d;

	    for (int i = 0; i < len && lb < best_so_far; i++)
	    {
	    	double u_i = u.get(order[i]);
	    	double l_i = l.get(order[i]);
	        uu = (u_i-mean)/std;
	        ll = (l_i-mean)/std;
	        d = 0;
	        if (qo[i] > uu)
	            d = dist(qo[i], uu);
	        else
	        {   if(qo[i] < ll)
	            d = dist(qo[i], ll);
	        }
	        lb += d;
	        cb[order[i]] = d;
	    }
	    return lb;
	}

	/// Calculate Dynamic Time Wrapping distance
	/// A,B: data and query, respectively
	/// cb : cummulative bound used for early abandoning
	/// r  : size of Sakoe-Chiba warpping band
	double dtw( double[] A
			  , double[] B
			  , double[] cb
			  , int m
			  , int r
			  , double bsf
			  )
	{

	    double[] cost_tmp;
	    int i,j,k;
	    double x,y,z,min_cost;

	    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
	    double[] cost = new double[2*r + 1];//(double*)malloc(sizeof(double)*(2*r+1));
	    for(k=0; k<2*r+1; k++)    cost[k]=INF;

	    double[] cost_prev = new double[2*r + 1];//(double*)malloc(sizeof(double)*(2*r+1));
	    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;

	    for (i=0; i<m; i++)
	    {
	        k = Math.max(0,r-i);
	        min_cost = INF;

	        for(j=Math.max(0,i-r); j<=Math.min(m-1,i+r); j++, k++)
	        {
	            /// Initialize all row and column
	            if ((i==0)&&(j==0))
	            {
	                cost[k]=dist(A[0],B[0]);
	                min_cost = cost[k];
	                continue;
	            }

	            if ((j-1<0)||(k-1<0))     y = INF;
	            else                      y = cost[k-1];
	            if ((i-1<0)||(k+1>2*r))   x = INF;
	            else                      x = cost_prev[k+1];
	            if ((i-1<0)||(j-1<0))     z = INF;
	            else                      z = cost_prev[k];

	            /// Classic DTW calculation
	            cost[k] = Math.min( Math.min( x, y) , z) + dist(A[i],B[j]);

	            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
	            if (cost[k] < min_cost)
	            {   min_cost = cost[k];
	            }
	        }

	        /// We can abandon early if the current cumulative distance with lower bound together are larger than bsf
	        if (i+r < m-1 && min_cost + cb[i+r+1] >= bsf)
	        {   
	            return min_cost + cb[i+r+1];
	        }

	        /// Move current array to previous array.
	        cost_tmp = cost;
	        cost = cost_prev;
	        cost_prev = cost_tmp;
	    }
	    k--;

	    // the DTW distance is in the last cell in the matrix of size O(m^2) 
	    // or at the middle of our array.
	    double final_dtw = cost_prev[k];
	   
	    return final_dtw;
	}

	

	@Override
	public Location query(double[] q
						 , DataInputStream data
						 , double warpingWindow
						 ) 
	{
		
	
	   
	    double bsf;          /// best-so-far
	    
	    int[] order;          ///new order of the query
	    double[] u;
	    double[] l;
	    double[] qo;
	    double[] uo;
	    double[] lo;
	    double[] tz;
	    double[] cb;
	    double[] cb1;
	    double[] cb2;
	    double[] t;
	    double[] l_buff;
	    double[] u_buff;

	    double d;
	    int i , j;
	    double ex , ex2 , mean, std;
	    int m=-1, r=-1;
	    long loc = 0;
	    int kim = 0,keogh = 0, keogh2 = 0;
	    double dist=0, lb_kim=0, lb_k=0, lb_k2=0;
	    double[] buffer;
	    
	    Index[] Q_tmp;

	    

	   
	    /// read size of the query
	    m = q.length;

	    /// read warping windows
	    
	    {   double R = warpingWindow;
	        if (R<=1)
	            r = (int) Math.floor(R*m);
	        else
	            r = (int) Math.floor(R);
	    }

	    
	    qo = new double[m];
	    uo = new double[m];
	    lo = new double[m];

	    order = new int[m];

	    Q_tmp = new Index[m];
	    for(int k = 0;k < Q_tmp.length;++k)
	    {
	    	Q_tmp[k] = new Index();
	    }

	    u = new double[m];

	    l = new double[m];

	    cb = new double[m];

	    cb1 = new double[m];

	    cb2 = new double[m];

	   

	    t = new double[m*2];

	    tz = new double[m];

	    buffer = new double[EPOCH];

	    u_buff = new double[EPOCH];

	    l_buff = new double[EPOCH];


	    /// Read query file
	    bsf = INF;
	    i = 0;
	    j = 0;
	    ex = ex2 = 0;
	    for(;i < m;++i)
	    {
	    	d = q[i];
	    	ex += d;
	    	ex2 += d*d;
	    	
	    }
	  

	    /// Do z-normalize the query, keep in same array, q
	    mean = ex/m;
	    std = ex2/m;
	    std = Math.sqrt(std-mean*mean);
	    for( i = 0 ; i < m ; i++ )
	    {
	    	
	         q[i] = (q[i] - mean)/std;
	    }

	    /// Create envelop of the query: lower envelop, l, and upper envelop, u
	    Bounds bounds = computeLemireBounds(q, m, r, new Bounds(m));
	    l = bounds.getLower();
	    u = bounds.getUpper();

	    /// Sort the query one time by abs(z-norm(q[i]))
	    for( i = 0; i<m; i++)
	    {
	        Q_tmp[i].value = q[i];
	        Q_tmp[i].index = i;
	    }
	    
	    Arrays.sort(Q_tmp);
	    
	    /// also create another arrays for keeping sorted envelop
	    for( i=0; i<m; i++)
	    {   int o = Q_tmp[i].index;
	        order[i] = o;
	        qo[i] = q[o];
	        uo[i] = u[o];
	        lo[i] = l[o];
	    }
	    
	    /// Initial the cummulative lower bound
	    for( i=0; i<m; i++)
	    {   cb[i]=0;
	        cb1[i]=0;
	        cb2[i]=0;
	    }

	    i = 0;          /// current index of the data in current chunk of size EPOCH
	    j = 0;          /// the starting index of the data in the circular array, t
	    ex = ex2 = 0;
	    boolean done = false;
	    int it=0, ep=0, k=0;
	    int I;    /// the starting index of the data in current chunk of size EPOCH

	    while(!done)
	    {
	        /// Read first m-1 points
	        ep=0;
	        if (it==0)
	        {   
	        	for(k=0; k<m-1; k++)
	        	{
	        		try
	        		{
	        			d = data.readDouble();
	        			buffer[k] = d;
	        		}
	        		catch(IOException ioe)
	        		{
	        			
	        		}
	                
	        	}
	        }
	        else
	        {   for(k=0; k<m-1; k++)
	                buffer[k] = buffer[EPOCH-m+1+k];
	        }

	        /// Read buffer of size EPOCH or when all data has been read.
	        ep=m-1;
	        while(ep<EPOCH)
	        {   
	        	
	        	try
	        	{
	        		d = data.readDouble();
	        		buffer[ep] = d;
	        		ep++;
	        	}
	        	catch(IOException ioe)
	        	{
	        		break;
	        	}
	        }

	        /// Data are read in chunk of size EPOCH.
	        /// When there is nothing to read, the loop is end.
	        if (ep<=m-1)
	        {   done = true;
	        } else
	        {   
	        	
	        	Bounds bnds = computeLemireBounds(buffer, ep, r, new Bounds(l_buff, u_buff));
	        	l_buff = bnds.getLower();
	        	u_buff = bnds.getUpper();
	        	
	            /// Just for printing a dot for approximate a million point. Not much accurate.
//	            if (it%(1000000/(EPOCH-m+1))==0)
//	                fprintf(stderr,".");

	            /// Do main task here..
	            ex=0; ex2=0;
	            for(i=0; i<ep; i++)
	            {
	                /// A bunch of data has been read and pick one of them at a time to use
	                d = buffer[i];

	                /// Calcualte sum and sum square
	                ex += d;
	                ex2 += d*d;

	                /// t is a circular array for keeping current data
	                t[i%m] = d;

	                /// Double the size for avoiding using modulo "%" operator
	                t[(i%m)+m] = d;

	                /// Start the task when there are more than m-1 points in the current chunk
	                if( i >= m-1 )
	                {
	                    mean = ex/m;
	                    std = ex2/m;
	                    std = Math.sqrt(std-mean*mean);

	                    /// compute the start location of the data in the current circular array, t
	                    j = (i+1)%m;
	                    /// the start location of the data in the current chunk
	                    I = i-(m-1);

	                    /// Use a constant lower bound to prune the obvious subsequence
	                    lb_kim = lowerBoundKimHierarchy(t, q, j, m, mean, std, bsf);

	                    if (lb_kim < bsf)
	                    {
	                        /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
	                        /// uo, lo are envelop of the query.
	                        lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, m, mean, std, bsf);
	                        if (lb_k < bsf)
	                        {
	                            /// Take another linear time to compute z_normalization of t.
	                            /// Note that for better optimization, this can merge to the previous function.
	                            for(k=0;k<m;k++)
	                            {   tz[k] = (t[(k+j)] - mean)/std;
	                            }

	                            /// Use another lb_keogh to prune
	                            /// qo is the sorted query. tz is unsorted z_normalized data.
	                            /// l_buff, u_buff are big envelop for all data in this chunk
	                            DoubleBuffer l_buff_i = DoubleBuffer.wrap(l_buff, I, l_buff.length - I);
	                            DoubleBuffer u_buff_i = DoubleBuffer.wrap(u_buff, I, u_buff.length - I);
	                            lb_k2 = lb_keogh_data_cumulative( order
	                            								, tz
	                            								, qo
	                            								, cb2
	                            								, l_buff_i
	                            								, u_buff_i
	                            								, m
	                            								, mean
	                            								, std
	                            								, bsf
	                            								);
	                            if (lb_k2 < bsf)
	                            {
	                                /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
	                                /// Note that cb and cb2 will be cumulative summed here.
	                                if (lb_k > lb_k2)
	                                {
	                                    cb[m-1]=cb1[m-1];
	                                    for(k=m-2; k>=0; k--)
	                                        cb[k] = cb[k+1]+cb1[k];
	                                }
	                                else
	                                {
	                                    cb[m-1]=cb2[m-1];
	                                    for(k=m-2; k>=0; k--)
	                                        cb[k] = cb[k+1]+cb2[k];
	                                }

	                                /// Compute DTW and early abandoning if possible
	                                dist = dtw(tz, q, cb, m, r, bsf);

	                                if( dist < bsf )
	                                {   /// Update bsf
	                                    /// loc is the real starting location of the nearest neighbor in the file
	                                    bsf = dist;
	                                    loc = (it)*(EPOCH-m+1) + i-m+1;
	                                }
	                            } else
	                                keogh2++;
	                        } else
	                            keogh++;
	                    } else
	                        kim++;

	                    /// Reduce obsolute points from sum and sum square
	                    ex -= t[j];
	                    ex2 -= t[j]*t[j];
	                }
	            }

	            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
	            if (ep<EPOCH)
	                done=true;
	            else
	                it++;
	        }
	    }

	    i = (it)*(EPOCH-m+1) + ep;
	    
	    Location ret = new Location();
	    ret.offset = loc;
	    ret.distance = Math.sqrt(bsf);
	    ret.dataScanned = i;
	    return ret;
/*
	    /// Note that loc and i are long long.
	    cout << "Location : " << loc << endl;
	    cout << "Distance : " << sqrt(bsf) << endl;
	    cout << "Data Scanned : " << i << endl;
	    cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;

	    /// printf is just easier for formating ;)
	    printf("\n");
	    printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i)*100);
	    printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i)*100);
	    printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i)*100);
	    printf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2)/i*100));
	    return 0;
	    */
	}
}
