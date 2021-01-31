#include <iostream>
#include <stdio.h>
#include "functions.h"
#include "class_args.h"
//#include <sched.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/types.h>
//#include <unistd.h>
#include <pthread.h>

using namespace std;



int main(int argc, char* argv[])
{	
	if( argc < 6 || argc > 7) 
	{
		std::cout<<"argc is less or greater then needed"<<endl;

		return 1;
	}

	int r = 0, n = 0, m = 0, l = 0, q = 0, s = 0, p = 0, k = 0; 

	if( !( sscanf(argv[1], "%d", &n) == 1 && n > 0 
		&& sscanf(argv[2], "%d", &m) == 1 && m > 0 
		&& sscanf(argv[3], "%d", &p) == 1 && p > 0 
		&& sscanf(argv[4], "%d", &r) == 1 && r >= 0 
		&& sscanf(argv[5], "%d", &s) == 1 && s >= 0 && s < 5) )
	{

		if( sscanf(argv[1], "%d", &n) == 0 )
		{
			printf("n should be an integer\n");
			
			return 1;
		}

		if( sscanf(argv[2], "%d", &m) == 0 )
		{
			printf("m should be an integer\n");
			
			return 1;
		}

		if( sscanf(argv[4], "%d", &r) == 0 )
		{
			printf("r should be an integer\n");
			
			return 1;
		}

		if( sscanf(argv[5], "%d", &s) == 0 )
		{
			printf("s should be an integer\n");
			
			return 1;
		}

		if( n == 0) 
		{
			printf("As n = 0 there are no equations to solve\n");
			
			return 1;
		}

		if( m == 0) printf("m, which equals 0, is not an acceptable parameter for the algorithm\n");
		return 1;
	}

	if( (atoi(argv[2]) > atoi(argv[1])) )
	{
		if( atoi(argv[2]) > atoi(argv[1]) ) std::cout<<"m is greater then n, m = "<<m<<", n = "<<n<<endl;

		return 1;
	}

	if( atoi(argv[4]) > atoi(argv[1]) ) 
	{
		if(n <= 8) r = n;
		else r = 8;
	}

	double * b, *a, *x, *inverse_matr; //, *main_mass;
	args *arg;
	pthread_t *tids;

	arg = new args[p];
	if(!arg)
    {
        printf("Can't allocate memory for args\n");
        return 1;
    }

	q = n/m;
	l = n - m*q;

	int min;

	if(l == 0) min = m;
	else min = l;

	if(p > q )	p = q;

	tids = new pthread_t[p - 1];
	
	if(!tids)
	{
		delete [] arg;
		printf("Can't allocate memory for tids \n");
		return 1;
	}

	a = new double[n*n]; 

	if(!a)
	{
		delete [] arg;
		delete [] tids;
		printf("Can't allocate memory for matrix \n");
		return 1;
	}

	// main_mass = new double[n*m];

	// if(!main_mass)
	// {
	// 	delete [] arg;
	// 	delete [] tids;
	// 	delete [] a;

	// 	printf("Can't allocate memory for shared buffer \n");
	// 	return 1;
	// }

	b = new double[n];

	if(!b)
	{
		delete [] arg;
		delete [] tids;
		delete [] a;
		//delete [] main_mass;

		printf("Can't allocate memory for vector b \n");
		return 1;
	}

	x = new double[n];

	if(!x)
	{
		delete [] arg;
		delete [] tids;
		delete [] a;
		//delete [] main_mass;
		delete [] b;

		printf("Can't allocate memory for vector of solution x \n");
		return 1;
	}

	inverse_matr = new double[m*m]; // shared block for A_ss ^ (-1) 

	if(!inverse_matr)
	{
		delete [] arg;
		delete [] tids;
		delete [] a;
		//delete [] main_mass;
		delete [] b;
		delete [] x;

		printf("Can't allocate memory for shared buffer, which is used for an inverse matrix \n");
		return 1;
	}

	for(k = 0; k < p; k ++)
	{
		arg[k].a = a;
		arg[k].b = b;
		arg[k].x = x;
		arg[k].n = n;
		arg[k].m = m;
		arg[k].k = k;
		arg[k].l = l;
		arg[k].min = min;
		arg[k].p = p;
		arg[k].q = q;
		arg[k].s = s;
		arg[k].r = r;
		//arg[k].main_mass = main_mass;
		arg[k].inverse_matr = inverse_matr;
		arg[k].argv = argv[0];
	}
	
	if(argc == 6) // set matrix using formula
	{
		if( s == 0)
		{
			std::cout<<" Wrong s: if you want to input matrix using the function, you should enter s such that 0 < s < 5 "<<endl;
			printf(" Usage: %s n m r s [name].txt\n", argv[0]);

			delete [] a;
			delete [] b;
			delete [] x;
			delete [] inverse_matr;
			delete [] tids;
			delete [] arg;
			//delete [] main_mass;

			return 1;
		}

		int flag = 0;

		for(k = 0; k < p - 1; k ++)
		{
			flag = pthread_create(&tids[k], 0, &initialization_formula, arg + k);
	        if( flag != 0 )
	        {
	        	printf("Can't create thread %d for task 'find_max_min'\n", k);
	        }
		}

		initialization_formula(arg + p - 1);

		for(k = 0; k < p - 1; k ++)
		{
			if(pthread_join(tids[k], 0))
	            fprintf(stderr, "Can't wait thread %d for task 'find_max_min'\n", k);
	
	        //printf("For process number %d the time of calculating is %10.8e sec\n", k, arg[k].elapsed_time);

            if(flag == 0) flag = arg[k].error;
		}

		if(flag == 0) flag = arg[p - 1].error;

		//printf("For process number %d the time of calculating is %10.8e sec\n", p - 1, arg[p - 1].elapsed_time);

		if(flag != 0)
		{
			delete [] a;
			delete [] b;
			delete [] x;
			delete [] inverse_matr;
			delete [] tids;
			delete [] arg;
			//delete [] main_mass;

			return flag;
		}


	}

	if(argc == 7) 
	{
		const char* filename;
		if( s != 0)
		{
			std::cout<<" Wrong s: if you want to input matrix from the txt file, you should enter s as 0 "<<endl;
			printf(" Usage: %s n m r s [name].txt\n", argv[0]);

			delete [] a;
			delete [] b;
			delete [] x;
			delete [] inverse_matr;
			delete [] tids;
			delete [] arg;
			//delete [] main_mass;

			return 1;
		}

		filename = argv[6];

		for(k = 0; k < p; k ++)
		{
			arg[k].filename = filename;
		}

		int flag = 0;

		for(k = 0; k < p - 1; k ++)
		{
			flag = pthread_create(&tids[k], 0, &initialization_file, arg + k);
	        if( flag != 0 )
	        {
	        	printf("Can't create thread %d for task 'find_max_min'\n", k);
	        }
		}

		initialization_file(arg + p - 1);

		for(k = 0; k < p - 1; k ++)
		{
			if(pthread_join(tids[k], 0))
	            fprintf(stderr, "Can't wait thread %d for task 'find_max_min'\n", k);
	
	        //printf("For process number %d the time of calculating is %10.8e sec\n", k, arg[k].elapsed_time);

            if(flag == 0) flag = arg[k].error;
		}

		//printf("For process number %d the time of calculating is %10.8e sec\n", p - 1, arg[p - 1].elapsed_time);

		if(flag == 0) flag = arg[p - 1].error;

		if(flag != 0)
		{
			delete [] a;
			delete [] b;
			delete [] x;
			delete [] inverse_matr;
			delete [] tids;
			delete [] arg;
			//delete [] main_mass;

			return flag;
		}
	}

	printf("\n\n");

	std::cout<<" Vector x:"<<endl;                    // vector x
	for(int i = 0; i < r; i++) printf("%10.3e\n", x[i]);
	printf("\n");   
	printf("\n\n"); 

	for(k = 0; k < p; k ++)
	{
        printf(" For process number %d the time of calculating is %10.8e sec\n", k, arg[k].elapsed_time);
	}
	printf("\n\n");

	printf(" Total time of calculating of a solution is %10.8e sec\n", arg[0].wall_clock_time);
	printf(" Total time of calculating of a discrepancy is %10.8e sec\n", arg[0].wall_clock_time_discperancy);

	printf("\n\n");
	
	printf("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d p = %d\n", argv[0], arg[0].global_discperancy, arg[0].wall_clock_time, s, n, m, p);

	delete [] a;
	delete [] b;
	delete [] x;
	delete [] inverse_matr;
	delete [] tids;
	delete [] arg;
	//delete [] main_mass;

	return 0;
}


