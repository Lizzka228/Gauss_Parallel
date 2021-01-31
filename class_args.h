#ifndef __CLASS_ARGS_H__
#define __CLASS_ARGS_H__
using namespace std;

class args{
	public:
        int p = 0;
        int k = 0;
	int n = 0;
	int m = 0;
	int l = 0;
	int q = 0;
	int s = 0;
	int r = 0;
	int last_k = 0;
	int last_k_do_dif = 0;
	int min = 0;
	int curr_step = 0;
	int additioal_block_str1 = 0;
	int additioal_block_str2 = 0;
	int error = 0;
	double *a = 0;
	double *b = 0;
	double *x = 0;
	double *v1 = 0;
	double *v2 = 0;
	double *v = 0;
	double *main_mass = 0;
	double *part_b = 0;
	double *global_part_b = 0;
	double *block_row = 0;
	double *inverse_matr = 0;
	double elapsed_time = 0.0;
 	double wall_clock_time = 0.0;
 	double elapsed_time_discperancy = 0.0;
 	double wall_clock_time_discperancy = 0.0;
 	double norm = 0.0;
 	double eps = 0.0;
 	double global_discperancy = 0.0;
 	double global_margin = 0.0;
 	const char *filename = 0;
 	const char *argv = 0;
	args() = default;
};
 #endif 
