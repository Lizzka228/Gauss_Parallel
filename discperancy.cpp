#include <iostream>
#include "functions.h"
#include "class_args.h"
#include <pthread.h>
#include <cstring>
#include <math.h>


using namespace std;


void get_b_copy(args *a, int p)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    double *block_row = a->block_row, *b = a->b, *x = a->x;
    int i, n = a->n;
    
    pthread_mutex_lock(&m);

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &m);
    
    for(i = 0; i < n; i ++)
    {	
    	block_row[i] = b[i];
    	block_row[i + n] = x[i];
    }
    
    t_out ++;
    
    if(t_out >= p)
    {
        t_in = 0;
        //p_a = 0;
        pthread_cond_broadcast(&c_out);
    }
    else while(t_out < p) pthread_cond_wait(&c_out, &m);
    
    pthread_mutex_unlock(&m);
}


void reduce_diff(args *a, int p, double res, double margin)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static double global_margin = 0.0;
    static double global_discperancy = 0.0;
    double *b = a->b;
    int i, n = a->n;
    static double norm_b = 0.0;
    static args *p_a = 0;
    
    pthread_mutex_lock(&m);

    if(!p_a) 
    {
    	p_a = a;

    	for(i = 0; i < n; i ++)	norm_b += b[i] * b[i];

    	norm_b = sqrt(norm_b);
    }

    global_discperancy += res;
    global_margin += margin;

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &m);
 
    a->global_discperancy = sqrt(global_discperancy) / norm_b;
    a->global_margin = sqrt(global_margin);
    
    t_out ++;
    
    if(t_out >= p)
    {
        t_in = 0;
        p_a = 0;
        pthread_cond_broadcast(&c_out);
    }
    else while(t_out < p) pthread_cond_wait(&c_out, &m);
    
    pthread_mutex_unlock(&m);
}



void discperancy_m_greater_then_1(args * arg) // norma nevyazki
{
	int p = arg->p, n = arg->n, i, j, r, k = arg->q, m = arg->m, proc_num = arg->k, last_k = arg->last_k;
	double *a = arg->a, *b = arg->block_row, *x = arg->block_row + n;
	int min = arg->min;
	double sum = 0.0;
	double res = 0.0; // |Ax|^2
	double margin = 0.0; // |x_my - x_origin| ^ 2

	get_b_copy(arg, p);

	int t = 0;

	for( i = proc_num; i < k; i += p )
	{	
		a = arg->a + i * n * m;
		b = arg->block_row + i*m;
		t = 0;

		for( r = 0; r < m; r ++ )
		{
			sum = 0.0;
			for(j = 0; j < n; j += 2) 
			{
				sum += a[r * n + j] * x[j];
				margin += (x[j] - 1.0) * (x[j] - 1.0);
			}

			for(j = 1; j < n; j += 2) 
			{
				sum += a[r * n + j] * x[j];
				margin += (x[j]) * (x[j]);
			}

			sum -= b[t];

			t++;

			res += sum * sum;
		}
	}

	for( i = k; i < last_k; i += p )
	{	
		a = arg->a + i * n * m;
		b = arg->block_row + i*m;

		t = 0;

		for( r = 0; r < min; r ++ )
		{
			sum = 0.0;

			for(j = 0; j < n; j += 2) 
			{
				sum += a[r * n + j] * x[j];
				margin += (x[j] - 1.0) * (x[j] - 1.0);
			}

			for(j = 1; j < n; j += 2) 
			{
				sum += a[r * n + j] * x[j];
				margin += (x[j]) * (x[j]);
			}

			sum -= b[t];

			t++;

			res += sum * sum;
		}
	}

	reduce_diff(arg , p, res, margin);
}


void get_b_copy_less(args *a, int p, double *x1)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    double *block_row = a->block_row, *b = a->b, *x = a->x;;
    int i, n = a->n;
    
    pthread_mutex_lock(&m);

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &m);
    
    for(i = 0; i < n; i ++)
    {   
        block_row[i] = b[i];
        x1[i] = x[i];
    }
    
    t_out ++;
    
    if(t_out >= p)
    {
        t_in = 0;
        //p_a = 0;
        pthread_cond_broadcast(&c_out);
    }
    else while(t_out < p) pthread_cond_wait(&c_out, &m);
    
    pthread_mutex_unlock(&m);
}


void discperancy_m_less_then_1(args * arg) // norma nevyazki
{
    int p = arg->p, n = arg->n, i, j, r, k = arg->q, m = arg->m, proc_num = arg->k, last_k = arg->last_k;
    double *a = arg->a, *b = arg->block_row, *x = new double[n];
    int min = arg->min;
    double sum = 0.0;
    double res = 0.0; // |Ax|^2
    double margin = 0.0; // |x_my - x_origin| ^ 2

    get_b_copy_less(arg, p, x);

    int t = 0;

    for( i = proc_num; i < k; i += p )
    {   
        a = arg->a + i * n * m;
        b = arg->block_row + i*m;
        t = 0;

        for( r = 0; r < m; r ++ )
        {
            sum = 0.0;
            for(j = 0; j < n; j += 2) 
            {
                sum += a[r * n + j] * x[j];
                margin += (x[j] - 1.0) * (x[j] - 1.0);
            }

            for(j = 1; j < n; j += 2) 
            {
                sum += a[r * n + j] * x[j];
                margin += (x[j]) * (x[j]);
            }

            sum -= b[t];

            t++;

            res += sum * sum;
        }
    }

    for( i = k; i < last_k; i += p )
    {   
        a = arg->a + i * n * m;
        b = arg->block_row + i*m;

        t = 0;

        for( r = 0; r < min; r ++ )
        {
            sum = 0.0;

            for(j = 0; j < n; j += 2) 
            {
                sum += a[r * n + j] * x[j];
                margin += (x[j] - 1.0) * (x[j] - 1.0);
            }

            for(j = 1; j < n; j += 2) 
            {
                sum += a[r * n + j] * x[j];
                margin += (x[j]) * (x[j]);
            }

            sum -= b[t];

            t++;

            res += sum * sum;
        }
    }

    reduce_diff(arg , p, res, margin);

    delete [] x;
}


