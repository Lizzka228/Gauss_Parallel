#include <iostream>
#include "functions.h"
#include "class_args.h"
#include <pthread.h>
#include <cstring>
#include <math.h>


using namespace std;

void reduce_inverse_min(int p, args *arg, int s, int flag)
{
    static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static args * p_a = 0;
    double * A, *block_row = arg->block_row, *inverse_matr, *v1 = arg->v1, *part_b = arg->part_b, *b;
    int proc_num = arg->k;
    int n, m , i, j, r, min;
    int static flag_global;
    
    pthread_mutex_lock(&mut);

    if(!p_a && proc_num == s % p)
    {
        p_a = arg;
        n = arg->n;
        m = arg->m;
        min = arg->min;
        A = arg->a + s*n*m;
        b = arg->b + s*m;

        inverse_matr = arg->inverse_matr;

        flag_global = flag;

        arg->error = flag_global;
    }

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &mut);

    n = arg->n;
    m = arg->m;
    min = arg->min;
    A = p_a->a + s*n*m;
    b = p_a->b + s*m;

    inverse_matr = arg->inverse_matr;

    arg->error = flag_global;

    for(i = 0; i < min; i ++)
    {
        part_b[i] = b[i];
        for(j = 0; j < n; j ++) block_row[i*n + j] = A[i*n + j];
        for(r = 0; r < min; r ++) v1[i*min + r] = inverse_matr[i*min + r];
    }
    
    t_out ++;
    
    if(t_out >= p)
    {
        p_a = 0;
        t_in = 0;
        pthread_cond_broadcast(&c_out);
    }
    else while(t_out < p) pthread_cond_wait(&c_out, &mut);
    
    pthread_mutex_unlock(&mut); 
}


void reduce_inverse (int p, args *arg, int s, int flag)
{
    static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static args * p_a = 0;
    double /** A, *block_row = arg->block_row, */ *inverse_matrix, *v1 = arg->v1, *part_b = arg->part_b, *b;
    //int n;
    int proc_num = arg->k;
    int n, m , i, j, r;
    int static flag_global;

    (void)(j);
    (void)(n);
    
    pthread_mutex_lock(&mut);

    if(!p_a && proc_num == s % p)
    {
        p_a = arg;

        m = arg->m;

        flag = inverse_matr(s, arg, m);

        //A = arg->a + s*n*m;

        //inverse_matrix = arg->inverse_matr;

        flag_global = flag;

        //arg->error = flag_global;
    }

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &mut);

    n = arg->n;
    m = arg->m;
    //A = p_a->a + s*n*m;
    b = p_a->b + s*m;

    inverse_matrix = arg->inverse_matr;

    arg->error = flag_global;

    for(i = 0; i < m; i ++)
    {
        part_b[i] = b[i];
        //for(j = 0; j < n; j ++) block_row[i*n + j] = A[i*n + j];
        for(r = 0; r < m; r ++) v1[i*m + r] = inverse_matrix[i*m + r];
    }
    
    t_out ++;
    
    if(t_out >= p)
    {
        p_a = 0;
        t_in = 0;
        pthread_cond_broadcast(&c_out);
    }
    else while(t_out < p) pthread_cond_wait(&c_out, &mut);
    
    pthread_mutex_unlock(&mut); 
    
}



void reduce_mult_row (int p, args *arg, int s)
{
    static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    //static args * p_a = 0;
    double * A, *block_row = arg->block_row;
    int n, m , i, r;
    n = arg->n;
    m = arg->m;
    
    pthread_mutex_lock(&mut);

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &mut);

    n = arg->n;
    m = arg->m;
    A = arg->a + s*n*m;
    block_row = arg->block_row;

    for(i = 0; i < m; i ++)
    {
        for(r = 0; r < n; r ++)
        {
            block_row[i*n + r] = A[i*n + r];
        }
    }
 
    t_out ++;

    //do_diff0(s, arg);
    
    if(t_out >= p)
    {
        t_in = 0;
        pthread_cond_broadcast(&c_out);
    }
    else while(t_out < p) pthread_cond_wait(&c_out, &mut);  
    
    pthread_mutex_unlock(&mut); 

    //do_diff0(s, arg);
}

void reduce_b(args * arg, int s, int p)
{
    static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static args *p_a = 0;
    int m = arg->m, i;
    double *b, *v_2 = arg->v2, *x = arg->x + s*m;
    //int proc_num = arg->k;
    
    pthread_mutex_lock(&mut);

    if(!p_a )
    {
        p_a = arg;
        b = arg->b + s*m;
        x = arg->x + s*m;

        for(i = 0; i < m; i++) 
        {
            x[i] = b[i] - v_2[i];
        }
    }

    t_in ++;

    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else while(t_in < p) pthread_cond_wait(&c_in, &mut);

    if(p_a != arg) for(i = 0; i < m; i++) x[i] -= v_2[i];

    t_out ++;
    
    if(t_out >= p)
    {
        t_in = 0;
        p_a = 0;
        pthread_cond_broadcast(&c_out);
    }
    else while(t_out < p) pthread_cond_wait(&c_out, &mut);
    
    pthread_mutex_unlock(&mut);

}