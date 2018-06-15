#include<vector>

class particle {
public:

	double ro;
	std::vector<double>     x_coordinat;
	double              h_radius_kernel;
	double                            m;
	double                            V;

	double                           W(particle x);
	double	                W_gradient(particle x);

};

class iron_particle : public particle {
public:

	double  ro_particle(iron_particle *X, int N);

	double  Viscosity(iron_particle X, double gamma);

	iron_particle* particle_initialization(std::vector<std::vector<double> > x,int N,double h_const,double m){

		iron_particle *X;
		X = new iron_particle[N ];
		for (int i = 0; i < N; i++) {
			X[i].x_coordinat = x[i];
			X[i].m = m;
			X[i].h_radius_kernel = h_const;
		}
		
		return(X);
	}
};






