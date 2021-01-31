#include <iostream> 
#include "functions.h"
#include <math.h>
#include <sys/time.h>
#include "class_args.h"
#include <pthread.h>

using namespace std;

void mult_row0(int s, args * arg)
{
	int q, y, m = arg->m, n = arg->n, k = arg->q, p = arg->p, r, z, j, q1 = arg->k;
	double c00 = 0.0, c01 = 0.0, c02 = 0.0, c10 = 0.0, c11 = 0.0, c12 = 0.0, c20 = 0.0, c21 = 0.0, c22 = 0.0;
	double *p_a = arg->block_row, *v_1 = arg->v1, *v_2 = arg->v2, *a = arg->a;
	int integer = s/p;
	int begin = p * (integer) + q1;
	int last_k = arg->last_k;
	int min = arg->min, t = 0;

	(void)(p_a);

	if(s % p >= q1) 
	{
		begin = p * (integer + 1) + q1;
	}

	if(begin >= last_k) return;

	for( j = begin; j < k; j += p )
	{
		//p_a = arg->block_row + j*m;
		a = arg->a + s*n*m + j*m;
					
		for(r = 0; r < m; r += 3)
		{		
			for(y = 0; y < m; y += 3)
			{
				for(z = 0; z < m; z ++)
				{
					c00 += v_1[r * m + z] * a[z * n + y];
					c01 += v_1[r * m + z] * a[z * n + y + 1];
					c02 += v_1[r * m + z] * a[z * n + y + 2];
				
					c10 += v_1[r * m + m + z] * a[z * n + y];
					c11 += v_1[r * m + m + z] * a[z * n + y + 1];
					c12 += v_1[r * m + m + z] * a[z * n + y + 2];
				
					c20 += v_1[r * m + 2*m + z] * a[z * n + y];
					c21 += v_1[r * m + 2*m + z] * a[z * n + y + 1];
					c22 += v_1[r * m + 2*m + z] * a[z * n + y + 2];
				}
				
				v_2[r * m + y]     = c00;
				v_2[r * m + y + 1] = c01;
				v_2[r * m + y + 2] = c02;
							
				v_2[r * m + m + y]     = c10;
				v_2[r * m + m + y + 1] = c11;
				v_2[r * m + m + y + 2] = c12;
				
				v_2[r * m + 2*m + y]     = c20;
				v_2[r * m + 2*m + y + 1] = c21;
				v_2[r * m + 2*m + y + 2] = c22;
				
				c00 = 0.0;
				c01 = 0.0;
				c02 = 0.0;
				
				c10 = 0.0;
				c11 = 0.0;
				c12 = 0.0;
				
				c20 = 0.0;
				c21 = 0.0;
				c22 = 0.0;
			}
		}

		a = arg->a + s*n*m + j*m;	

		t = 0;

		for(r = 0; r <  m; r ++) 
			for(y = 0; y < m; y ++)	
			{
				a[r * n + y] = v_2[t];	
				t++;
			}	
	}

	for( j = k; j < last_k; j += p )
	{
		//printf("ssssssssss\n");
		//p_a = arg->block_row + j*m;
		a = arg->a + s*n*m + j*m;

		for(r = 0; r < m; r ++)
		{
			for(z = 0; z < min; z ++ )
			{
				c00 = 0.0;
	
				for(q = 0; q < m; q ++ ) c00 += (v_1[r*m + q]) * (a[ q*n + z ]); //(v_1[q*min + t]);
	
				v_2[r*min + z] = c00;
			}
		}	

		a = arg->a + s*n*m + j*m;

		t = 0;

		for(r = 0; r <  m; r ++) 
			for(y = 0; y < min; y ++)	
			{
				a[r * n + y] = v_2[t];
				t++;
			}
	}
}


void do_diff0(int s, args * arg)
{
	int q, y, t, m = arg->m, s1, n = arg->n, k = arg->q, p = arg->p, r, z, j, q1 = arg->k, i;
	double c00 = 0.0, c01 = 0.0, c02 = 0.0, c10 = 0.0, c11 = 0.0, c12 = 0.0, c20 = 0.0, c21 = 0.0, c22 = 0.0;
	double *p_a, *v_1 = arg->v1, *v_2 = arg->v2, *b = arg->b;
	double *part_b = arg->part_b;
	double *block_row = arg->block_row;
	int integer = s/p;
	int begin = p * (integer) + q1;
	double sum = 0.0;
	int last_k = arg->last_k;
	int min = arg->min;
	int last_k_do_dif = arg->last_k_do_dif;

	if(s % p >= q1 ) 
	{
		begin = p * (integer + 1) + q1;
	}

	if(begin >= last_k_do_dif) return;

	for( i = begin; i < k; i += p )
	{	
		t = 0;

		p_a = arg->a + (i*n+s)*m;
		
		for( s1 = 0; s1 < m; s1++)
		{
			for(q = 0; q < m; q++) 
			{
				v_1[t] = p_a[ s1*n + q ];
				p_a[ s1*n + q ] = 0.0;
				t++;
			}
		}
		
		for(s1 = 0; s1 < m; s1 ++)
		{
			sum = 0.0;
		
			for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (part_b[q]);
		
			v_2[s1] = sum;
		}
		
		for(s1 = 0; s1 < m; s1 ++)
			b[i*m + s1] -= v_2[s1];
		
		for( j = s + 1; j < k; j ++ )
		{
			block_row = arg->block_row + j*m;

			p_a = arg->a + (i*n+j)*m;

			for(r = 0; r < m; r += 3)
			{
					
				for(y = 0; y < m; y += 3)
				{
					for(z = 0; z < m; z ++)
					{
					
						c00 += v_1[r * m + z] * block_row[z * n + y];
						c01 += v_1[r * m + z] * block_row[z * n + y + 1];
						c02 += v_1[r * m + z] * block_row[z * n + y + 2];
					
						c10 += v_1[r * m + m + z] * block_row[z * n + y];
						c11 += v_1[r * m + m + z] * block_row[z * n + y + 1];
						c12 += v_1[r * m + m + z] * block_row[z * n + y + 2];
					
						c20 += v_1[r * m + 2*m + z] * block_row[z * n + y];
						c21 += v_1[r * m + 2*m + z] * block_row[z * n + y + 1];
						c22 += v_1[r * m + 2*m + z] * block_row[z * n + y + 2];
					}
					
					p_a[r * n + y]     -= c00;
					p_a[r * n + y + 1] -= c01;
					p_a[r * n + y + 2] -= c02;
								
					p_a[r * n + n + y]     -= c10;
					p_a[r * n + n + y + 1] -= c11;
					p_a[r * n + n + y + 2] -= c12;
					
					p_a[r * n + 2*n + y]     -= c20;
					p_a[r * n + 2*n + y + 1] -= c21;
					p_a[r * n + 2*n + y + 2] -= c22;
					
					c00 = 0.0;
					c01 = 0.0;
					c02 = 0.0;
					
					c10 = 0.0;
					c11 = 0.0;
					c12 = 0.0;
					
					c20 = 0.0;
					c21 = 0.0;
					c22 = 0.0;
				}
			}
		}

		for( j = k; j < last_k_do_dif; j ++ )
		{
			block_row = arg->block_row + j*m;

			p_a = arg->a + (i*n+j)*m;

			for(s1 = 0; s1 < m; s1 ++)
			{
				for(t = 0; t < min; t ++ )
				{
					sum = 0.0;
	
					for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (block_row[q*n + t]);
	
					p_a[s1*n + t] -= sum;
				}
			}
		}
	}

	// row l by m

	for( i = k; i < last_k; i ++ )
	{	
		t = 0;

		p_a = arg->a + (i*n+s)*m;
		
		for( s1 = 0; s1 < min; s1++)
		{
			for(q = 0; q < m; q++) 
			{
				v_1[t] = p_a[ s1*n + q ];
				p_a[ s1*n + q ] = 0.0;
				t++;
			}
		}
		
		for(s1 = 0; s1 < min; s1 ++)
		{
			sum = 0.0;
		
			for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (part_b[q]);
		
			v_2[s1] = sum;
		}
		
		for(s1 = 0; s1 < min; s1 ++)
			b[i*m + s1] -= v_2[s1];
		
		for( j = s + 1; j < k; j ++ )
		{
			block_row = arg->block_row + j*m;
			p_a = arg->a + (i*n+j)*m;

			for(s1 = 0; s1 < min; s1 ++)
			{
				for(t = 0; t < m; t ++ )
				{
					sum = 0.0;
	
					for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (block_row[q*n + t]);
	
					p_a[s1*n + t] -= sum;
				}
			}
		}
		for( j = k; j < last_k_do_dif; j ++ )
		{
			block_row = arg->block_row + j*m;
			p_a = arg->a + (i*n+j)*m;

			for(s1 = 0; s1 < min; s1 ++)
			{
				for(t = 0; t < min; t ++ )
				{
					sum = 0.0;
	
					for(q = 0; q < m; q ++ ) sum += (v_1[s1*m + q]) * (block_row[q*n + t]);
	
					p_a[s1*n + t] -= sum;
				}
			}
		}
	}
}