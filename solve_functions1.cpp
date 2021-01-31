#include <iostream> 
#include "functions.h"
#include <math.h>
#include <sys/time.h>
#include "class_args.h"
#include <pthread.h>

using namespace std;

void mult_row1(int s, args * arg)
{
	int q, y, r_m, z_m, m = arg->m, n = arg->n, k = arg->q, p = arg->p, r, z, j, q1 = arg->k;
	double c00 = 0.0, c01 = 0.0, c02 = 0.0, c10 = 0.0, c11 = 0.0, c12 = 0.0, c20 = 0.0, c21 = 0.0, c22 = 0.0;
	double *v_1 = arg->v1, *a = arg->a, *v_2 = arg->v2;
	int integer = s/p;
	int begin = p * (integer) + q1;
	int last_k = arg->last_k;
	int min = arg->min;
	int m_new = m - 1;
	//int r_n;
	int t, s1;

	if(s % p >= q1) 
	{
		begin = p * (integer + 1) + q1;
	}

	if(begin >= last_k) return;

	for( j = begin; j < k; j += p )
	{
		//p_a = arg->block_row + j*m;
		a = arg->a + s*n*m + j*m;
					
		for(r = 0; r < m_new; r += 3)
		{
			r_m = r * m;
			//r_n = r * n;
				
			for(y = 0; y < m_new; y += 3)
			{
				for(z = 0; z < m; z ++)
				{
					z_m = z * n;
				
					c00 += v_1[r_m + z] * a[z_m + y];
					c01 += v_1[r_m + z] * a[z_m + y + 1];
					c02 += v_1[r_m + z] * a[z_m + y + 2];
				
					c10 += v_1[r_m + m + z] * a[z_m + y];
					c11 += v_1[r_m + m + z] * a[z_m + y + 1];
					c12 += v_1[r_m + m + z] * a[z_m + y + 2];
				
					c20 += v_1[r_m + 2*m + z] * a[z_m + y];
					c21 += v_1[r_m + 2*m + z] * a[z_m + y + 1];
					c22 += v_1[r_m + 2*m + z] * a[z_m + y + 2];
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

			y = m_new;
				
			for(z = 0; z < m; z ++)
			{
				z_m = z * n;
				
				c00 += v_1[r_m + z] * a[z_m + y];
				
				c10 += v_1[r_m + m + z] * a[z_m + y];
				
				c20 += v_1[r_m + 2*m + z] * a[z_m + y];
			}
				
			v_2[r * m + y]     = c00;
						
			v_2[r * m + m + y] = c10;
						
			v_2[r * m + 2*m + y] = c20;
						
			c00 = 0.0;
						
			c10 = 0.0;
						
			c20 = 0.0;
		}
	
		r_m = m_new * m;
		//r_n = m_new * n;
					
		for(y = 0; y < m_new; y += 3)
		{
			for(z = 0; z < m; z ++)
			{
				z_m = z * n;
					
				c00 += v_1[r_m + z] * a[z_m + y];
				c01 += v_1[r_m + z] * a[z_m + y + 1];
				c02 += v_1[r_m + z] * a[z_m + y + 2];
			}
					
			v_2[m_new * m + y]     = c00;
			v_2[m_new * m + y + 1] = c01;
			v_2[m_new * m + y + 2] = c02;
					
			c00 = 0.0;
			c01 = 0.0;
			c02 = 0.0;
		}

		y = m_new;
		
		for(z = 0; z < m; z ++)
		{
			z_m = z * n;
					
			c00 += v_1[r_m + z] * a[z_m + y];
					
			//c10 += v_1[r_m + m + z] * p_a[z_m + y];
		}
					
		v_2[m_new * m + y] = c00;
					
		c00 = 0.0;
					
		c10 = 0.0;

		t = 0;
		
		for( s1 = 0; s1 < m; s1++)
		{
			for(q = 0; q < m; q++) 
			{
				a[ s1*n + q ] = v_2[t];
				t++;
			}
		}

		// p_a = arg->block_row + j*m;
	
		// t = 0;
		
		// for( s1 = 0; s1 < m; s1++)
		// {
		// 	for(q = 0; q < m; q++) 
		// 	{
		// 		p_a[ s1*n + q ] = v_2[t];
		// 		t++;
		// 	}
		// }
					//t = 0;
					
	}

	for( j = k; j < last_k; j += p )
	{
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

		t = 0;
		
		for( s1 = 0; s1 < m; s1++)
		{
			for(q = 0; q < min; q++) 
			{
				a[ s1*n + q ] = v_2[t];
				t++;
			}
		}		

		// p_a = arg->block_row + j*m;
	
		// t = 0;
		
		// for( s1 = 0; s1 < m; s1++)
		// {
		// 	for(q = 0; q < min; q++) 
		// 	{
		// 		p_a[ s1*n + q ] = v_2[t];
		// 		t++;
		// 	}
		// }			
	}
}


void do_diff1(int s, args * arg)
{
	int q, y, r_m, z_m, t, m = arg->m, s1, n = arg->n, k = arg->q, p = arg->p, r, z, j, q1 = arg->k, i;
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
	int m_new = m - 1;
	int r_n;

	if(s % p >= q1 ) 
	{
		begin = p * (integer + 1) + q1;
	}

	//if(q1 == 0) printf("s = %d begin %d\n", s, begin);

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
		
					//t = 0;
		
		for( j = s + 1; j < k; j ++ )
		{
			block_row = arg->block_row + j*m;
			p_a = arg->a + (i*n+j)*m;

			for(r = 0; r < m_new; r += 3)
			{
				r_m = r * m;
				r_n = r * n;
					
				for(y = 0; y < m_new; y += 3)
				{
					for(z = 0; z < m; z ++)
					{
						z_m = z * n;
					
						c00 += v_1[r_m + z] * block_row[z_m + y];
						c01 += v_1[r_m + z] * block_row[z_m + y + 1];
						c02 += v_1[r_m + z] * block_row[z_m + y + 2];
					
						c10 += v_1[r_m + m + z] * block_row[z_m + y];
						c11 += v_1[r_m + m + z] * block_row[z_m + y + 1];
						c12 += v_1[r_m + m + z] * block_row[z_m + y + 2];
					
						c20 += v_1[r_m + 2*m + z] * block_row[z_m + y];
						c21 += v_1[r_m + 2*m + z] * block_row[z_m + y + 1];
						c22 += v_1[r_m + 2*m + z] * block_row[z_m + y + 2];
					}
					
					p_a[r_n + y]     -= c00;
					p_a[r_n + y + 1] -= c01;
					p_a[r_n + y + 2] -= c02;
								
					p_a[r_n + n + y]     -= c10;
					p_a[r_n + n + y + 1] -= c11;
					p_a[r_n + n + y + 2] -= c12;
					
					p_a[r_n + 2*n + y]     -= c20;
					p_a[r_n + 2*n + y + 1] -= c21;
					p_a[r_n + 2*n + y + 2] -= c22;
					
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
					
				for(z = 0; z < m; z ++)
				{
					z_m = z * n;
					
					c00 += v_1[r_m + z] * block_row[z_m + y];
					
					c10 += v_1[r_m + m + z] * block_row[z_m + y];
					
					c20 += v_1[r_m + 2*m + z] * block_row[z_m + y];
				}
					
				p_a[r_n + y]     -= c00;
							
				p_a[r_n + n + y] -= c10;
							
				p_a[r_n + 2*n + y] -= c20;
							
				c00 = 0.0;
							
				c10 = 0.0;
							
				c20 = 0.0;
			}
		
			r_m = m_new * m;
			r_n = m_new * n;
						
			for(y = 0; y < m_new; y += 3)
			{
				for(z = 0; z < m; z ++)
				{
					z_m = z * n;
						
					c00 += v_1[r_m + z] * block_row[z_m + y];
					c01 += v_1[r_m + z] * block_row[z_m + y + 1];
					c02 += v_1[r_m + z] * block_row[z_m + y + 2];
				}
						
				p_a[r_n + y]     -= c00;
				p_a[r_n + y + 1] -= c01;
				p_a[r_n + y + 2] -= c02;
						
				c00 = 0.0;
				c01 = 0.0;
				c02 = 0.0;
			}
			for(z = 0; z < m; z ++)
			{
				z_m = z * n;
						
				c00 += v_1[r_m + z] * block_row[z_m + y];
						
				//c10 += v_1[r_m + m + z] * block_row[z_m + y];
			}
						
			p_a[r_n + y] -= c00;
						
			c00 = 0.0;
						
			//c10 = 0.0;
		
			// p_a = arg->a + (i*n+j)*m;
		
			// for(s1 = 0; s1 < m; s1 ++)
			// 	for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
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
		
			// p_a = arg->a + (i*n+j)*m;
		
			// for(s1 = 0; s1 < m; s1 ++)
			// 	for(t = 0; t < min; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
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
		
			// p_a = arg->a + (i*n+j)*m;
		
			// for(s1 = 0; s1 < min; s1 ++)
			// 	for(t = 0; t < m; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*m + t];
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
		
			// p_a = arg->a + (i*n+j)*m;
		
			// for(s1 = 0; s1 < min; s1 ++)
			// 	for(t = 0; t < min; t ++ ) p_a[s1*n + t] = p_a[s1*n + t] - v_2[s1*min + t];
		}
	}
}