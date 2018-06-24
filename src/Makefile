main: main.cpp main.h
	g++ main.cpp -o main -O3 -Wall -std=c++11 -fopenmp

main_gpu: main_gpu.cu main.h
	nvcc main_gpu.cu -o main_gpu -v --compiler-options -Wall,-O3 -Xptxas -O3 -arch=compute_61
