/*****************************************************************************
Name		: 	main.cpp
Author      : 	Flores, Facundo Gabriel
e-mail		:   flores.facundogabriel@gmail.com
Version     : 	0.1
Description : 	This is the NSL main program. It will show to the user different
				options.
License     :
    			Copyright (C) 2012 Flores, Facundo Gabriel

    			This program is free software: you can redistribute it and/or modify
    			it under the terms of the GNU General Public License as published by
    			the Free Software Foundation, either version 3 of the License, or
    			(at your option) any later version.

    			This program is distributed in the hope that it will be useful,
    			but WITHOUT ANY WARRANTY; without even the implied warranty of
    			MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    			GNU General Public License for more details.

    			You should have received a copy of the GNU General Public License
    			along with this program.  If not, see <http://www.gnu.org/licenses/>.

**************************************************************************** */

/*#include <iostream>
#include <vector>
#include <cudaParallel.hpp>
#include <nonParallel.hpp>
#include <linearMatrix.hpp>

#include <pcg.hpp>

#include <sys/time.h>*/

#include <stdio.h>
#include "testing.hpp"

void mypaint(void)
{
	printf("******************************\n");

}

void print_menu()
{
	puts("Elija una opcion. Para salir inserte 0");
	mypaint();
	printf("1- Test Jacobi\n");
	printf("2- Test BiCG\n");
	printf("3- Test CG\n");
	printf("4- Test BiCGStab\n");
	printf("5- Test Preconditionated\n");
	mypaint();
}

#include "testing.hpp"

int main(){

	int opt;
	print_menu();
	scanf("%d", &opt);
	while(opt > 0 && opt < 6)
	{
		test_option(opt);
		print_menu();
		scanf("%d", &opt);
	}
	if(opt == 0)
		puts("Terminado!");

	//test_option(1);
	return 0;
}

/*int main()
{
	float myarray[] = {1,2,3,2,8,4,3,4,11};
	linearMatrix A(3,3);
	linearMatrix L(3,3);

	vector<float> myvector(myarray, myarray + 9);

	A.Set_Matrix(myvector, 3, 3);

	cudaCholesky(A, L);

	L.Show_Little_Matrix();
}*/

/*int main()
{
    timeval start, end;

    long seconds, useconds, final;

	linearMatrix A(2000,2000);
	A.Generate_Diagonal_Dominant();
	linearMatrix L(2000,2000);

	gettimeofday(&start, 0);
	cudaCholesky(A,L);
	gettimeofday(&end, 0);

	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;

	cout << ((seconds) * 1000 + useconds/1000.0) + 0.5 << endl;
}*/

