/*SPH метод, наковальня*/

#include<stdio.h>
#include<conio.h>
#include <iostream>
#include <fstream>
#include "particle.h"

#include <time.h> 

void main() {
	clock_t elapsedTime;
	elapsedTime = clock();

	int N = 1320;
	int N_x;
	int N_y;
	double h_const = 0.076;

	double density_0 = 7.85;
	double V_0 = 200;

	double high = 2.54;
	double thick = 0.76;
	double area = high*thick;

	N_y = sqrt(N*high / thick);
	N_x =  N / N_y;



	std::vector<std::vector<double> > x(N + N_x + N_y + 1);   //координаты каждой частицы газа
	std::cout << N_x << " " << N_y << " " << N_x*N_y << "\n";

	for (int i = 0; i <N; i++) {
		x[i].resize(2);
		x[i][0] = (i - (i / N_x)*N_x + 1)* thick / double(N_x);
		x[i][1] = (i / N_x + 1)*high / double(N_y);
	}

	for (int i = N; i <= N + N_x; i++) {
		x[i].resize(2);
		x[i][0] = (i - N)* thick / double(N_x);
		x[i][1] = 0;
	}

	for (int i = N + N_x + 1; i <= N + N_x + N_y; i++) {
		x[i].resize(2);
		x[i][0] = 0;
		x[i][1] = (i - N - N_x)*high / double(N_y);
	}

	N = N + N_x + N_y + 1;
	double m = density_0*area / double(N);

	iron_particle tmp_particle_initialization;
	iron_particle *X=NULL;

	/*инициализация частиц*/
	// X = tmp_particle_initialization.particle_initialization(x, N, h_const, m); // Old code.
	X = X->particle_initialization(x, N, h_const, m);

	for (int i = 0; i < N; i++) {
		X[i].ro = X[i].ro_particle(X, N);
		std::cout << double(i) / double(N) * 100 << "% \n";
	}

	/*метод SPH*/
	double Viscosity = 0;
	double gradient = 0;

	std::ofstream output_t_0;
	output_t_0.open("output_t_0.txt");


	for (int i = 0; i < N; i++) {
		output_t_0 << X[i].x_coordinat[0] << " " << X[i].x_coordinat[1] << " " << X[i].ro << "\n";
	}
	elapsedTime = clock() - elapsedTime;
	printf("It took me %d clicks (%f seconds).\n", elapsedTime, ((float)elapsedTime) / CLOCKS_PER_SEC);

	_getch();
}