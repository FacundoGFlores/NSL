Description
-----------

This program was presented as a final project of the subject
called: "Seminar" in April 2013. 

The project implements many iterative methods for solving sys-
tem of linear equations.

Usage
-----

I implemented a generic code for showing the use of the different
created libraries. I use a generic matrix for illustrative purpose.
You can take any of the libraries for solving a specific kind of
matrix. 

Example:

	linearMatrix A(N,N);
	linearMatrix B(N,1);
	linearMatrix X(N,1);
	
	//Call a rutine to prepare matrices
	//Prepare(A,B,X);
	//Remember you can use the following method from "linearMatrix.hpp"
	//void Read_File_Matrix(const char *Path);
	//Read documentation for more info about matrix-file structure

	Parallel P(A,B,X, error, iter); //Object with parallel methods
	P.With_Jacobi(); //P.With_"METHOD_NAME"
	
	linearMatrix SOLUTION(N,1);
	SOLUTION = P.Get_Solution();
	
	//example: SOLUTION.Write_File_Matrix("Solution.dat");

Repo Installation
-----------------

I used the Nsight Eclipse IDE on Linux. I tested the project on the
5 and 5.5 editions.

- Copy the NSL repo on your own cuda-workspace
- Go to File -> Import -> Existing Projects into Workspace
- Select the NSL repo
- Click Finish

Compiling
---------

To configure the compiler please go:

- Project -> Properties -> Build -> Settings 

Then edit into the NVCC Compiler the "All options" section:
	
	-l/home/YOUR_HOME/cuda-workspace/NSL/include -G -g -O0

Errata
------

In the perfomance timing analysis I put the starting clocks
after the cublas structures initializations, so for better
performance timing analysis put the clocks before the cublas
initializations.





