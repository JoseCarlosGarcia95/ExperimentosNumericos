default: main

main:
	gcc -lm main_noparalell.c -o noparallel -Wall -O1
	gcc -fopenmp -lm main_paralell.c -o parallel -Wall -O1
	nvcc -lm main_paralell_gpu.cu -o parallel_gpu -ccbin clang-3.8


