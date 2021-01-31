#include <iostream>
#include "functions.h"
#include "class_args.h"
#include <sched.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <pthread.h>
#include <cstring>
#include <math.h>
#include <sys/resource.h>


#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

using namespace std;

double get_time()
{
	struct rusage r;
 	
 	getrusage(RUSAGE_THREAD, &r);

	return r.ru_utime.tv_sec + r.ru_utime.tv_usec/1000000.0;	
}


double get_full_time()
{
	struct timeval r;
 	
 	gettimeofday(&r, 0);

	return r.tv_sec + r.tv_usec/1000000.0;	
}


void define_last_k(args * a)
{
	int l = a->l;

	if(l == 0)
	{
		a->last_k = a->q;
		a->last_k_do_dif = a->q;

		return;
	}
	else
	{
		int p = a->p, k = a->k, q = a->q, count = q / p, remain = q - count * p;
		a->last_k_do_dif = a->q + 1;

		if(k == remain) a->last_k = a->q + 1;
		else a->last_k = a->q;
	}
}

void recover_A_and_b_file(int p, args *a)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static args *p_a = 0;
    double *A = a->a, *b = a->b;
    
    pthread_mutex_lock(&m);

    if(!p_a) 
    {
    	p_a = a;
    }

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &m);

    if(p_a == a) 
    {
		enter_matrix_from_file(A, a->n, a->filename);
		vector(a->n, A, b);
    }
    
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

void recover_A_and_b(int p, args *a)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static args *p_a = 0;
    double *A = a->a, *b = a->b;
    
    pthread_mutex_lock(&m);

    if(!p_a) 
    {
    	p_a = a;
    }

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &m);

    if(p_a == a) 
    {
		enter_matrix(a->n, a->m, a->q, a->s, a->l, A);
		vector(a->n, A, b);
    }
    
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


void reduce_norm (int p, args *a)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static args *p_a = 0;
    static int norm;
    
    pthread_mutex_lock(&m);

    if(!p_a) 
    {
    	p_a = a;
		norm = a->norm;
    }
    else 
    {
		if(a->norm > norm) norm = a->norm;
    }

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &m);

    //printf("a->min_max_sum = p_a->min_max_sum is %lf \n\n", p_a->min_max_sum);
    
    a->norm = norm;
    
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


void find_norm(args * a)
{
	double norm = 0;
	double sum = 0.0;
	int p = a->p, k = a->k, n = a->n, m = a->m, l = a->l, q = a->q, count = q / p, remain = q - count * p, i, j, z;
	double *p_a = a->a;

	//printf("remain = %d, l = %d, q = %d, count = %d\n", remain, l, q, count);

	if(k*m < n - l)
	{
		//printf("lessssssssssss for proc %d\n", k);
		for(j = 0; j < n; j ++) norm += fabs(p_a[k * m * n + j]);
	
		for(i = 1; i < m; i ++)
		{
			sum = 0.0;
	
			for(j = 0; j < n; j ++) sum += fabs(p_a[k * m * n + i * n + j]);
	
			if( sum > norm ) norm = sum;
		}
	}

	for(i = 1; i < count; i ++)
	{
		for(z = 1; z < m; z ++)
		{
			sum = 0.0;

			for(j = 0; j < n; j ++) sum += fabs(p_a[ i * p * m * n + k * m * n + z * n + j]);

			if( sum > norm ) norm = sum;
		}
	}

	if(remain != 0 && k < remain && (remain + count * p)*n*m == n * n - l * n)
	{
		//printf("remain != 0 HEEEEEREEEEEEEEEE for proc %d\n", k);

		for(z = 1; z < m; z ++)
		{
			sum = 0.0;

			for(j = 0; j < n; j ++) sum += fabs(p_a[ count * p * m * n + k * m * n + z * n + j]);

			if( sum > norm ) norm = sum;
		}
	}

	if(l != 0 && k == remain)
	{
		//printf("l != 0 HEEEEEREEEEEEEEEE for proc %d\n", k);

		for(z = 0; z < l; z ++)
		{
			sum = 0.0;

			for(j = 0; j < n; j ++) sum += fabs(p_a[ count * p * m * n + k * m * n + z * n + j]);

			if( sum > norm ) norm = sum;
		}
	}

	a->norm = norm;

	//printf("The Norm for the introduced matrix equals %10.3e for proc %d\n", norm, k);
}



void wait(int p)
{
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    
    pthread_mutex_lock(&m);

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &m);
    
    t_out ++;
    
    if(t_out >= p)
    {
        t_in = 0;
        pthread_cond_broadcast(&c_out);
    }
    else while(t_out < p) pthread_cond_wait(&c_out, &m);
    
    pthread_mutex_unlock(&m);
    
    
}



void null_make(args * a)
{
	double *p_a = a->a, *p_b = a->b, *p_x = a->x;
	int p = a->p, k = a->k, n = a->n, m = a->m, l = a->l, q = a->q, count = q / p, remain = q - count * p, i;

	for(i = 0; i < count; i ++)
	{
		memset(p_a + i * p * m * n + k * m * n, 0.0, m * n * sizeof(double));  // a matrix A init in thread k
		memset(p_b + i * p * m + k * m, 0.0, m * sizeof(double));              // a vector B init in thread k
		memset(p_x + i * p * m + k * m, 0.0, m * sizeof(double));              // a vector X init in thread k
	}

	if(remain != 0 && k < remain && (remain + count * p)*n*m == n * n - l * n)
	{
		memset(p_a + count * p * m * n + k * m * n, 0.0, m * n * sizeof(double));   // a matrix A init in thread k
		memset(p_b + count * p * m + k * m, 0.0, m * sizeof(double));               // a vector B init in thread k
		memset(p_x + count * p * m + k * m, 0.0, m * sizeof(double));               // a vector X init in thread k
	}

	if(l != 0 && k == remain)
	{
		memset(p_a + count * p * m * n + k * m * n, 0.0, l * n * sizeof(double));   // a matrix A init in thread k
		memset(p_b + count * p * m + k * m, 0.0, l * sizeof(double));               // a vector B init in thread k
		memset(p_x + count * p * m + k * m, 0.0, l * sizeof(double));               // a vector X init in thread k
	}
}



void * initialization_formula(void * p1)
{	
	cpu_set_t cpu;
	int nprocs = get_nprocs();
	CPU_ZERO(&cpu);
	args * a = (args * )p1;
	CPU_SET(nprocs - 1 - (a->k), &cpu);
	sched_setaffinity(getpid(), sizeof(cpu), &cpu);

	static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
	static int flag1 = 0;
	double *p_a = a->a, *p_b = a->b, *p_inverse_matr = a->inverse_matr, *p_v1, *p_v2, *block_row, *part_b = a->part_b;
	int p = a->p, n = a->n, m = a->m, l = a->l, q = a->q, s = a->s, i;
	static double eps = 1.0;

	p_v1 = new double[m*m];  // additional block for an each thread
	p_v2 = new double[m*m];
	block_row = new double[n*m];
	part_b = new double[m];

	for (i = 0; i < m; ++i)
	{
		p_v1[i] = 0.0;
		p_v2[i] = 0.0;
		block_row[i] = 0.0;
		part_b[i] = 0.0;
	}
	for (i = m; i < m * m; ++i)
	{
		p_v1[i] = 0.0;
		p_v2[i] = 0.0;
		block_row[i] = 0.0;
	}
	for (i = m * m; i < n * m; ++i)	block_row[i] = 0.0;

	a->v1 = p_v1;
	a->v2 = p_v2;
	a->block_row = block_row;
	a->part_b = part_b;

	null_make(a);	//  each thread initialize it's own part of matrix, vector b and x

	define_last_k(a);

	wait(p);         //  wait for all threads to initialize their data, below I put elements using formula

	pthread_mutex_lock(&mut);

	if(flag1 == 0)
	{
		enter_matrix(n, m, q, s, l, p_a);

		printf(" Matrix\n");

		print_matrix(n, n, a->r, p_a);

		vector(n, p_a, p_b);

		std::cout<<" Vector b:"<<endl;                    // vector b

		for(int i = 0; i < a->r; i++) printf("%10.3e\n", p_b[i]);
		printf("\n");

		flag1 ++;

		for (int i = 0; i < m * m; ++i)
		{
			p_inverse_matr[i] = 0.0;
		}

		while( 1.0 + eps > 1.0 )
		{
			eps = eps / 2.0;
		}
	}

	a->eps = eps;

	pthread_mutex_unlock(&mut);

	find_norm(a);

	reduce_norm (p, a);

	a->elapsed_time = get_time();
    a->wall_clock_time = get_full_time();

	if(m % 3 == 0)	solve_0(a);
	if(m % 3 == 1)	solve_1(a);
	if(m % 3 == 2)	solve_2(a);

	a->elapsed_time = get_time() - a->elapsed_time;
    a->wall_clock_time = get_full_time() - a->wall_clock_time;

	recover_A_and_b(p, a);

	a->elapsed_time_discperancy = get_time();
    a->wall_clock_time_discperancy = get_full_time();

	if(m > 1 && (p > 1 || n < 50000))	discperancy_m_greater_then_1(a);
	if(m == 1 && (p > 1 || n < 50000))	discperancy_m_less_then_1(a);
	if(!(p > 1 || n < 50000))	a->global_discperancy = -1;

	a->elapsed_time_discperancy = get_time() - a->elapsed_time_discperancy;
    a->wall_clock_time_discperancy = get_full_time() - a->wall_clock_time_discperancy;

	wait(p);

	delete [] p_v1;
	delete [] p_v2;
	delete [] block_row;
	delete [] part_b;

    return 0;
}

void * initialization_file(void * p1)
{	
	cpu_set_t cpu;
	int nprocs = get_nprocs();
	CPU_ZERO(&cpu);
	args * a = (args * )p1;
	CPU_SET(nprocs - 1 - (a->k), &cpu);
	sched_setaffinity(getpid(), sizeof(cpu), &cpu);

	static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
	static int flag1 = 0;
	double *p_a = a->a, *p_b = a->b, *p_inverse_matr = a->inverse_matr, *p_v1, *p_v2, *block_row, *part_b = a->part_b;
	int p = a->p, n = a->n, m = a->m, i;
	static double eps = 1.0;

	p_v1 = new double[m*m];  // additional block for an each thread
	p_v2 = new double[m*m];
	block_row = new double[n*m];
	part_b = new double[m];

	for (i = 0; i < m; ++i)
	{
		p_v1[i] = 0.0;
		p_v2[i] = 0.0;
		block_row[i] = 0.0;
		part_b[i] = 0.0;
	}
	for (i = m; i < m * m; ++i)
	{
		p_v1[i] = 0.0;
		p_v2[i] = 0.0;
		block_row[i] = 0.0;
	}
	for (i = m * m; i < n * m; ++i)	block_row[i] = 0.0;

	a->v1 = p_v1;
	a->v2 = p_v2;
	a->block_row = block_row;
	a->part_b = part_b;

	null_make(a);	//  each thread initialize it's own part of matrix, vector b and x

	define_last_k(a);

	wait(p);         //  wait for all threads to initialize their data, below I put elements using formula

	pthread_mutex_lock(&mut);

	if(flag1 == 0)
	{
		flag1 = enter_matrix_from_file(p_a, n, a->filename);

		if(flag1 == 0)	

		{
			printf(" Matrix\n");
		
			print_matrix(n, n, a->r, p_a);
	
			vector(n, p_a, p_b);
	
			std::cout<<" Vector b:"<<endl;                    // vector b
	
			for(int i = 0; i < a->r; i++) printf("%10.3e\n", p_b[i]);
			printf("\n");
	
			flag1 ++;
	
			for (int i = 0; i < m * m; ++i)
			{
				p_inverse_matr[i] = 0.0;
			}
	
			while( 1.0 + eps > 1.0 )
			{
				eps = eps / 2.0;
			}
		}
	}

	a->eps = eps;

	pthread_mutex_unlock(&mut);

	wait(p);

	if(flag1 != 0 && flag1 != 1)
	{
		a->error = flag1;
		
		delete [] p_v1;
		delete [] p_v2;
		delete [] block_row;
		delete [] part_b;

	    return 0;
	}

	find_norm(a);

	reduce_norm (p, a);

	a->elapsed_time = get_time();
    a->wall_clock_time = get_full_time();

	if(m % 3 == 0)	solve_0(a);
	if(m % 3 == 1)	solve_1(a);
	if(m % 3 == 2)	solve_2(a);

	a->elapsed_time = get_time() - a->elapsed_time;
    a->wall_clock_time = get_full_time() - a->wall_clock_time;

	recover_A_and_b_file(p, a);

	a->elapsed_time_discperancy = get_time();
    a->wall_clock_time_discperancy = get_full_time();

	if(m > 1 && (p > 1 || n < 50000))	discperancy_m_greater_then_1(a);
	if(m == 1 && (p > 1 || n < 50000))	discperancy_m_less_then_1(a);
	if(!(p > 1 || n < 50000))	a->global_discperancy = -1;

	a->elapsed_time_discperancy = get_time() - a->elapsed_time_discperancy;
    a->wall_clock_time_discperancy = get_full_time() - a->wall_clock_time_discperancy;

	wait(p);

	delete [] p_v1;
	delete [] p_v2;
	delete [] block_row;
	delete [] part_b;

    return 0;
}




void vector(int n, double * a, double * b)  // make vector b
{
	for(int i = 0; i < n; i++)
	{
		b[i] = 0;

		for(int j = 0; j < (n+1)/2  ; j++)
		{
			b[i] += a[i*n + 2*j];
		}
	}
}


int print_matrix(int l, int n, int r, double * a) // print matrix l * n
{
	if ( r > n )
	{
		r = n;
	}

	if ( r > l )
	{
		r = l;
	}

	for(int i = 0; i < r; i ++)
	{
		for(int j = 0; j < r; j ++) printf("%10.3e ", a[i*n + j]);
		cout<<endl;
	}
	cout<<endl;

	return 0;
}

















