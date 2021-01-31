#include <iostream> 
#include "functions.h"
#include <math.h>
using namespace std;




void enter_matrix(int n, int m, int k, int s, int l, double * a) // set matrix, using function f
{
	int  a_v, a_h, last;
	double * p_a;
	int max;

	if( l!= 0 )
	{
		last = k + 1;
	}
	else last = k;
	for(int i = 0; i < last; i++)
	{
		if( i < k) a_v = m;
		else a_v = l;
		for(int j = 0; j < last; j++)
		{
			if( j < k) a_h = m;
			else a_h = l;

			p_a = a + (i*n+j)*m;

			for(int t = 0; t < a_v; t++)
			{
				for(int q = 0; q < a_h; q++) //p_a[ t*n + q ] = f(s, n, (i )*m + t + 1 , (j)*m + q + 1 );
				{
					if( s == 1) 
					{
						max = 0;
						if( ((i )*m + t + 1) <= ((j)*m + q + 1)) max = ((j)*m + q + 1);
						else max = ((i )*m + t + 1);
						//std::cout<<n - max + 1<<" = n - max + 1"<<endl;
						p_a[ t*n + q ] = n - max + 1;
					}
					if( s == 2) 
					{
						max = 0;
						if( ((i )*m + t + 1) <= ((j)*m + q + 1)) max = ((j)*m + q + 1);
						else max = ((i )*m + t + 1);
						p_a[ t*n + q ] = max;
					}
					if( s == 3 ) 
					{
						if( ((i )*m + t + 1)-((j)*m + q + 1) == 0)
						{
							p_a[ t*n + q ] = 0.0;
						}
						else
						{
							p_a[ t*n + q ] = abs( ( (i )*m + t + 1)-( (j)*m + q + 1) );
						}
					}
					if( s == 4 )
					{ 
						p_a[ t*n + q ] = 1.0/(((i )*m + t + 1)+((j)*m + q + 1)-1);
					}
				}
			}
		}
	}
}


int enter_matrix_from_file(double * a, int n, const char * filename)
{
	double c;
	
	FILE * f = fopen(filename, "r");

	if (f == NULL) 
    {
       printf ("Can't open file\n");
       return -1;
    }

    for (int i = 0; i < n*n; i++)
    {
        if( fscanf(f,"%lf", &c) != 1)
        {
        	std::cout <<"Wrong format of file contant or size of matrix is too big"<<endl;
        	
        	fclose(f);
        	return -1;
		}
        a[i] = c;
    }     
	
	fclose(f); 
	
	return 0;         
	
}