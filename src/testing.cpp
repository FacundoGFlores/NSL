#include "cudaParallel.hpp"
#include "nonParallel.hpp"
#include "linearMatrix.hpp"
#include "testing.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RUTA_JACOBI "jacobi.dat"
#define RUTA_BICG "bicg.dat"
#define RUTA_CG "cg.dat"
#define RUTA_BICGSTAB "bicgstab.dat"
#define RUTA_PRECOND "precond.dat"


void show_times(int N, double c_time, double cuda_time)
{
	puts("**********************");
	printf("Tamanio: %d\n", N);
	printf("CUDA: %f ms\n", cuda_time);
	printf("C: %f ms\n", c_time);
	printf("SPEEDUP: %fx\n", c_time/cuda_time);
	puts("**********************");
}

void init_file(const char *PATH)
{
	FILE *fp;
	fp = fopen(PATH, "w");
	if(fp != NULL)
	{
		fprintf(fp, "Orden\t\t\tC\t\t\tCUDA\t\t\tSPEEDUP\n" );
		fclose(fp);
	}
	else
	{
		puts("No se puede abrir el archivo");
	}

}

void save_times(int N, double c_time, double cuda_time, const char *PATH)
{
	FILE *fp;

	fp = fopen(PATH, "a+");
	if(fp != NULL)
	{
		fprintf(fp, "%d\t\t\t%f\t\t\t%f\t\t\t%.2f\n", N, c_time, cuda_time, c_time / cuda_time);
		fclose(fp);
	}
	else
	{
		puts("No se puede abrir el archivo");
	}

}

int check_solutions(linearMatrix &X_C, linearMatrix &X_CUDA)
{
	if(X_C.Get_Width() != X_CUDA.Get_Width())
		return -1;
	for(int i = 0; i < X_C.Get_Width(); i++)
	{
		float a = X_C.Get_Element(i, 0);
		float b = X_CUDA.Get_Element(i, 0);

		if( fabs(a - b) > 0.01f )
		{
			printf("Difference in Pos:%d->C:%f\tCUDA:%f\t dif:%f", i, a, b, fabs(a-b));
			return 0;
		}
	}
	return 1;
}

void prepare_matrices(linearMatrix &A, linearMatrix &B, linearMatrix &X)
{
	A.Generate_Diagonal_Dominant();
	B.Generate_Random_Matrix();
	X.Generate_Null_Matrix();
}

void prepare_matrices(linearMatrix &A, linearMatrix &B, linearMatrix &X, linearMatrix &P)
{
	A.Generate_Diagonal_Dominant();
	B.Generate_Random_Matrix();
	X.Generate_Null_Matrix();
	P.Generate_Idendity_Matrix();
}

void run_jacobi(const int N, const float error, const int iter,
				long &c_time, long &cuda_time)
{
	linearMatrix A(N, N);
	linearMatrix B(N,1);
	linearMatrix X(N,1);
	prepare_matrices(A, B, X);

	nonParallel nP(A, B, X, error, iter);
	nP.With_Jacobi();


	Parallel P(A, B, X, error, iter);
	P.With_Jacobi();

	linearMatrix C_SOL(N, 1);
	linearMatrix CUDA_SOL(N, 1);

	C_SOL = nP.Get_Solution();
	CUDA_SOL = P.Get_Solution();

	int err = check_solutions(C_SOL, CUDA_SOL);
	if(err >= 0)
	{
		printf("Checked solutions!\n");
		c_time = nP.Get_Time_Used();
		cuda_time = P.Get_Time_Used();
	}
	else
	{
		printf("Dimension errors");
	}
}

void test_Jacobi()
{
	long tiempo_cuda;
	long tiempo_c;

	init_file(RUTA_JACOBI);

	const float error = 0.001f;
	const int iter = 500;

	int N = 1000;
	run_jacobi(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_JACOBI);

	N = 2500;
	run_jacobi(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_JACOBI);

	N = 4000;
	run_jacobi(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_JACOBI);

	N = 5500;
	run_jacobi(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_JACOBI);

	N = 7000;
	run_jacobi(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_JACOBI);

	N = 8500;
	run_jacobi(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_JACOBI);

	N = 10000;
	run_jacobi(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_JACOBI);

	printf("Jacobi Test Finished!\n");

}


void run_bicg(const int N, const float error, const int iter,
				long &c_time, long &cuda_time)
{
	linearMatrix A(N, N);
	linearMatrix B(N,1);
	linearMatrix X(N,1);
	prepare_matrices(A, B, X);

	nonParallel nP(A, B, X, error, iter);
	nP.With_BICG();


	Parallel P(A, B, X, error, iter);
	P.With_BICG();

	linearMatrix C_SOL(N, 1);
	linearMatrix CUDA_SOL(N, 1);

	C_SOL = nP.Get_Solution();
	CUDA_SOL = P.Get_Solution();

	int err = check_solutions(C_SOL, CUDA_SOL);
	if(err >= 0)
	{
		printf("Checked solutions!\n");
		c_time = nP.Get_Time_Used();
		cuda_time = P.Get_Time_Used();
	}
	else
	{
		printf("Dimension errors");
	}
}

void test_bicg()
{
	long tiempo_cuda;
	long tiempo_c;

	init_file(RUTA_BICG);

	const float error = 0.001f;
	const int iter = 500;

	int N = 1000;
	run_bicg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICG);

	N = 2500;
	run_bicg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICG);

	N = 4000;
	run_bicg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICG);

	N = 5500;
	run_bicg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICG);

	N = 7000;
	run_bicg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICG);

	N = 8500;
	run_bicg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICG);

	N = 10000;
	run_bicg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICG);

	printf("BICG Test Finished!\n");

}

void run_cg(const int N, const float error, const int iter,
				long &c_time, long &cuda_time)
{
	linearMatrix A(N, N);
	linearMatrix B(N,1);
	linearMatrix X(N,1);
	prepare_matrices(A, B, X);

	nonParallel nP(A, B, X, error, iter);
	nP.With_Conjugate_Gradient();


	Parallel P(A, B, X, error, iter);
	P.With_Conjugate_Gradient();

	linearMatrix C_SOL(N, 1);
	linearMatrix CUDA_SOL(N, 1);

	C_SOL = nP.Get_Solution();
	CUDA_SOL = P.Get_Solution();

	int err = check_solutions(C_SOL, CUDA_SOL);
	if(err >= 0)
	{
		printf("Checked solutions!\n");
		c_time = nP.Get_Time_Used();
		cuda_time = P.Get_Time_Used();
	}
	else
	{
		printf("Dimension errors");
	}
}

void test_cg()
{
	long tiempo_cuda;
	long tiempo_c;

	init_file(RUTA_CG);

	const float error = 0.001f;
	const int iter = 500;

	int N = 1000;
	run_cg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_CG);

	N = 2500;
	run_cg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_CG);

	N = 4000;
	run_cg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_CG);

	N = 5500;
	run_cg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_CG);

	N = 7000;
	run_cg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_CG);

	N = 8500;
	run_cg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_CG);

	N = 10000;
	run_cg(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_CG);

	printf("CG Test Finished!\n");
}

void run_bicgstab(const int N, const float error, const int iter,
				long &c_time, long &cuda_time)
{
	linearMatrix A(N, N);
	linearMatrix B(N,1);
	linearMatrix X(N,1);
	prepare_matrices(A, B, X);

	nonParallel nP(A, B, X, error, iter);
	nP.With_Bicgstab();


	Parallel P(A, B, X, error, iter);
	P.With_BICGSTAB();

	linearMatrix C_SOL(N, 1);
	linearMatrix CUDA_SOL(N, 1);

	C_SOL = nP.Get_Solution();
	CUDA_SOL = P.Get_Solution();

	int err = check_solutions(C_SOL, CUDA_SOL);
	if(err >= 0)
	{
		printf("Checked solutions!\n");
		c_time = nP.Get_Time_Used();
		cuda_time = P.Get_Time_Used();
	}
	else
	{
		printf("Dimension errors");
	}
}

void test_bicgstab()
{
	long tiempo_cuda;
	long tiempo_c;

	init_file(RUTA_BICGSTAB);

	const float error = 0.001f;
	const int iter = 500;

	int N = 1000;
	run_bicgstab(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICGSTAB);

	N = 2500;
	run_bicgstab(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICGSTAB);

	N = 4000;
	run_bicgstab(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICGSTAB);

	N = 5500;
	run_bicgstab(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICGSTAB);

	N = 7000;
	run_bicgstab(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICGSTAB);

	N = 8500;
	run_bicgstab(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICGSTAB);

	N = 10000;
	run_bicgstab(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_BICGSTAB);

	printf("BICGSTAB Test Finished!\n");
}

void run_precond(const int N, const float error, const int iter,
				long &c_time, long &cuda_time)
{
	linearMatrix A(N, N);
	linearMatrix Q(N, N);
	linearMatrix B(N,1);
	linearMatrix X(N,1);

	prepare_matrices(A, B, X,Q);

	nonParallel nP(A, B, X, error, iter);
	nP.With_PCG(Q);

	Parallel P(A, B, X, error, iter);
	P.With_PCG(Q);

	linearMatrix C_SOL(N, 1);
	linearMatrix CUDA_SOL(N, 1);

	C_SOL = nP.Get_Solution();
	CUDA_SOL = P.Get_Solution();

	int err = check_solutions(C_SOL, CUDA_SOL);
	if(err >= 0)
	{
		printf("Checked solutions!\n");
		c_time = nP.Get_Time_Used();
		cuda_time = P.Get_Time_Used();
	}
	else
	{
		printf("Dimension errors");
	}
}

void test_precond()
{
	long tiempo_cuda;
	long tiempo_c;

	init_file(RUTA_PRECOND);

	const float error = 0.001f;
	const int iter = 500;

	int N = 1000;
	run_precond(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_PRECOND);

	N = 2500;
	run_precond(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_PRECOND);

	N = 4000;
	run_precond(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_PRECOND);

	N = 5500;
	run_precond(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_PRECOND);

	N = 7000;
	run_precond(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_PRECOND);

	N = 8500;
	run_precond(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_PRECOND);

	N = 10000;
	run_precond(N, 0.001f, 500, tiempo_c, tiempo_cuda);
	show_times(N, tiempo_c, tiempo_cuda);
	save_times(N, tiempo_c, tiempo_cuda, RUTA_PRECOND);

	printf("PRECONDITIONATED Test Finished!\n");
}

void test_option(const int opt)
{
	switch(opt)
	{
	case 1:
			test_Jacobi();
			break;
	case 2:
			test_bicg();
			break;
	case 3:
			test_cg();
			break;
	case 4:
			test_bicgstab();
			break;
	case 5:
			test_precond();
			break;
	default:
		printf("Incorrect option!\n");
		break;
	}
}

