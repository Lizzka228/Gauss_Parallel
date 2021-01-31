FLAGS = -Wunused -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -O3 -ffast-math


all: a.out

a.out: functions.o put_matrix.o parallel_solution.o solve_functions0.o solve_functions1.o solve_functions2.o solve_functions_common.o reduce_func_solve.o discperancy.o test.o
	g++ functions.o put_matrix.o parallel_solution.o solve_functions0.o solve_functions1.o solve_functions2.o solve_functions_common.o reduce_func_solve.o discperancy.o test.o -lpthread -o a.out 

discperancy.o: discperancy.cpp functions.h class_args.h
	g++ -c $(FLAGS) discperancy.cpp
	
reduce_func_solve.o: reduce_func_solve.cpp functions.h solve_functions_common.cpp class_args.h
	g++ -c $(FLAGS) reduce_func_solve.cpp
	
solve_functions0.o: solve_functions0.cpp functions.h class_args.h
	g++ -c $(FLAGS) solve_functions0.cpp
	
solve_functions1.o: solve_functions1.cpp functions.h class_args.h
	g++ -c $(FLAGS) solve_functions1.cpp

solve_functions2.o: solve_functions2.cpp functions.h class_args.h
	g++ -c $(FLAGS) solve_functions2.cpp
	
solve_functions_common.o: solve_functions_common.cpp functions.h class_args.h
	g++ -c $(FLAGS) solve_functions_common.cpp

functions.o: functions.cpp discperancy.cpp parallel_solution.cpp functions.h class_args.h 
	g++ -c $(FLAGS) functions.cpp

put_matrix.o: put_matrix.cpp functions.h class_args.h
	g++ -c $(FLAGS) put_matrix.cpp

parallel_solution.o: parallel_solution.cpp solve_functions_common.cpp solve_functions0.cpp solve_functions1.cpp solve_functions2.cpp class_args.h functions.h
	g++ -c $(FLAGS) parallel_solution.cpp

test.o: test.cpp functions.h class_args.h
	g++ -c $(FLAGS) test.cpp

clean:
	rm -rf *.0 a.out
