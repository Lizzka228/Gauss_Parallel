#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "class_args.h"
#include "functions.h"
#include <stdio.h>
#include <pthread.h>

double get_time();

double get_full_time();

void recover_A_and_b_file(int p, args *a);

void recover_A_and_b(int p, args *a);

void define_last_k(args * a);

int inverse_matr(int s, args * arg, int min);

void mult_row0(int s, args * arg);

void do_diff0(int s, args * arg);

void mult_row1(int s, args * arg);

void do_diff1(int s, args * arg);

void mult_row2(int s, args * arg);

void do_diff2(int s, args * arg);

void plus_step(int p, args * a);

void reduce_norm (int p, args *a);

int solve_0(args * arg);

int solve_1(args * arg);

int solve_2(args * arg);

void null_make(args * a);

void * initialization_formula(void * p1);

void * initialization_file(void * p1);

void wait(int p);

void vector(int n, double * a, double * b);

int print_matrix(int l, int n, int r, double * a);

void enter_matrix(int n, int m, int k, int s, int l, double * a);

int enter_matrix_from_file(double * a, int n, const char * filename);

int solution( int n, int m, double *a, double *b, double *x, double *v_1, double *v_2);

void print_discrepancy(int s5, int n,int m, double * a, double * b, double * x, double *v_1, double *v_2, float tv_seconds, const char * name);



void reduce_inverse_min(int p, args *arg, int s, int flag);

void reduce_inverse (int p, args *a, int s, int flag);

void reduce_mult_row (int p, args *arg, int s);

void do_X(int s, args * arg);

void reduce_b(args * arg, int s, int p);


void get_b_copy(args *a, int p);

void reduce_diff(args *a, int p, double res, double margin);

void discperancy_m_greater_then_1(args * arg);

void get_b_copy_less(args *a, int p, double *x);

void discperancy_m_less_then_1(args * arg);







#endif
