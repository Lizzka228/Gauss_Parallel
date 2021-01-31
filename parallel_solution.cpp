#include <iostream> 
#include "functions.h"
#include <math.h>
#include <sys/time.h>
#include "class_args.h"
#include <pthread.h>

using namespace std;

int solve_0(args * arg) // m % 3 == 0
{
	int p = arg->p, k = arg->q, proc_num = arg->k, s1 = 0, last; // local k is n/m, local q - thread number
	int flag = 0 ; 
	int min = arg->min, i, m = arg->m;

	if(arg->l == 0) last = k - 1;
	else last = k;

	for(s1 = 0; s1 < last; s1 ++)
	{
		reduce_inverse (p, arg, s1, flag);

		if(arg->error != 0) 
		{
			return -1;
		}

		mult_row0(s1, arg);

		reduce_mult_row (p, arg, s1);
		
		do_diff0(s1, arg);

		wait(p);
	} 

	(void)(flag);

	if(proc_num == s1 % p) 
	{

		flag = inverse_matr(s1, arg, min);
	
		double *x = arg->x + last*m, *b = arg->b + last*m;
	
		for(i = min - 1; i >= 0; i --) x[i] = b[i];
	}

	//reduce_inverse_min(p, arg, s1, flag);
	wait(p);

	if(arg->error != 0) 
	{
		return -1;
	}

	//wait(p);

	s1 = last - 1;

	//BEGIN X_i

	while(s1 >= 0)
	{
		do_X(s1, arg);

		reduce_b(arg, s1, p);

		s1--;
	}
	return 0;
}


int solve_1(args * arg) // m % 3 == 1
{
	int p = arg->p, k = arg->q, proc_num = arg->k, s1 = 0, last; // local k is n/m, local q - thread number
	int flag = 0 ; 
	int min = arg->min, i, m = arg->m;

	if(arg->l == 0) last = k - 1;
	else last = k;

	for(s1 = 0; s1 < last; s1 ++)
	{
		// if(proc_num == s1 % p) flag = inverse_matr(s1, arg, m);

		reduce_inverse (p, arg, s1, flag);

		if(arg->error != 0) 
		{
			return -1;
		}

		mult_row1(s1, arg);

		reduce_mult_row (p, arg, s1);
		
		do_diff1(s1, arg);

		//wait(p);
	} 

	if(proc_num == s1 % p) flag = inverse_matr(s1, arg, min);

	//reduce_inverse_min(p, arg, s1, flag);
	//wait(p);

	if(arg->error != 0) 
	{
		return -1;
	}

	if(proc_num == s1 % p) 
	{
		double *x = arg->x + last*m, *b = arg->b + last*m;

		for(i = min - 1; i >= 0; i --) x[i] = b[i];
	}

	wait(p);

	s1 = last - 1;

	//BEGIN X_i

	while(s1 >= 0)
	{
		do_X(s1, arg);

		reduce_b(arg, s1, p);

		s1--;
	}
	return 0;
}


int solve_2(args * arg) // m % 3 == 1
{
	int p = arg->p, k = arg->q, proc_num = arg->k, s1 = 0, last; // local k is n/m, local q - thread number
	int flag = 0 ; 
	int min = arg->min, i, m = arg->m;

	if(arg->l == 0) last = k - 1;
	else last = k;

	for(s1 = 0; s1 < last; s1 ++)
	{
		//if(proc_num == s1 % p) flag = inverse_matr(s1, arg, m);

		reduce_inverse (p, arg, s1, flag);

		if(arg->error != 0) 
		{
			return -1;
		}

		mult_row2(s1, arg);

		reduce_mult_row (p, arg, s1);
		
		do_diff2(s1, arg);

		//wait(p);
	} 

	if(proc_num == s1 % p) flag = inverse_matr(s1, arg, min);

	//reduce_inverse_min(p, arg, s1, flag);

	if(arg->error != 0) 
	{
		return -1;
	}

	if(proc_num == s1 % p) 
	{
		double *x = arg->x + last*m, *b = arg->b + last*m;

		for(i = min - 1; i >= 0; i --) x[i] = b[i];
	}

	wait(p);

	s1 = last - 1;

	//BEGIN X_i

	while(s1 >= 0)
	{
		do_X(s1, arg);

		reduce_b(arg, s1, p);

		s1--;
	}
	return 0;
}