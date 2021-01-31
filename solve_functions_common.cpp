#include <iostream> 
#include "functions.h"
#include <math.h>
#include <sys/time.h>
#include "class_args.h"
#include <pthread.h>

using namespace std;

int inverse_matr(int s, args * arg, int min)
{
	int d, t, m = arg->m, s1, n = arg->n, q;
	double norm = arg->norm, eps = 1.0, eps_new = 0.0, divisor = 0.0;
	double *p_a, *a = arg->a, *v_1 = arg->inverse_matr, *b = arg->b, *v_2 = arg->v2;
	//double *mass_main = arg->main_mass;

	while(eps + 1.0 > 1.0) eps = eps / 2.0;

	eps_new = eps * norm;

	p_a = a + (s*n+s)*m;

	for(t = 0; t < min; t ++ )
	{
		for( d = 0; d < t; d ++ ) v_1[t*min + d] = 0.0;
		v_1[t*min + t] = 1.0;
		for( d = t + 1; d < min; d ++ ) v_1[t*min + d] = 0.0;
	}

	for(d = 0; d < min; d++ )  // step number s
	{	
		divisor = p_a[d*n + d];

		if(fabs(divisor) <= eps_new)   // if a_ss == 0
		{
			printf("Used norm is column norm, according to the introduced matrix it equals %10.3e\n", norm);
			printf("One of block - matrices is degenerated: the main element is %10.3e, which is less then %10.3e\n", fabs(divisor), eps * norm );
			return 1;
		}
		
		// futher a_ss != 0
		
		for(t = 0; t < min; t++) // made a_ss equal to 1
		{
			p_a[d*n + t] /= divisor; 
		
			v_1[d*min + t] /= divisor; 
		}
		
		for(s1 = 0; s1 < d; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
		{
			divisor = p_a[s1*n + d];

			for( t = 0; t < min; t++ )	
			{
				p_a[s1*n + t] -= divisor * p_a[d*n + t];
		
				v_1[s1*min + t] -= divisor * v_1[d*min + t];
			}
		}
		for(s1 = d + 1; s1 < min; s1++) // for i from 0 to m - 1, i != s,  making a_is equal to 0
		{
			
			divisor = p_a[s1*n + d];
		
			for( t = 0; t < min; t++ )	
			{
				p_a[s1*n + t] -= divisor * p_a[d*n + t];
		
				v_1[s1*min + t] -= divisor * v_1[d*min + t];
			}
		}
	}

	for(s1 = 0; s1 < min; s1 ++)
	{
		divisor = 0.0;
		
		for(q = 0; q < min; q ++ ) 
		{
			divisor += (v_1[s1*min + q]) * (b[s*m + q]);
		}
	
		v_2[s1] = divisor;
	}
		
	for(s1 = 0; s1 < min; s1++) b[s*m + s1] = v_2[s1];

	return 0;
}



void do_X(int s, args * arg)
{
	double *v_2 = arg->v2, *a, *x = arg->x;
	int m = arg->m, j, p = arg->p, k = arg->q, proc_num = arg->k, n = arg->n, s1, q;
	double sum = 0.0;
	int min = arg->min;
	int last_k_do_dif = arg->last_k_do_dif;

	for(j = 0; j < m; j ++)	v_2[j] = 0.0;	

	for( j = k - proc_num - 1; j > s; j -= p )
	{
		a = arg->a + s*n*m + j*m;
		x = arg->x +j*m;
				
		for(s1 = 0; s1 < m; s1 ++)
		{			
				for(q = 0; q < m; q ++ ) sum += (a[s1*n + q]) * (x[q]);
		
				v_2[s1] += sum;		

				sum = 0.0;
		}
	}
	for( j = last_k_do_dif - proc_num - 1; j > k - 1; j -= p )
	{
		a = arg->a + s*n*m + j*m;
		x = arg->x +j*m;
				
		for(s1 = 0; s1 < m; s1 ++)
		{			
				for(q = 0; q < min; q ++ ) sum += (a[s1*n + q]) * (x[q]);
		
				v_2[s1] += sum;		

				sum = 0.0;
		}
	}
}