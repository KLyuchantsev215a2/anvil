#include"particle.h"

const double pi = 3.14159265358979;

double _pow(double base, int exponent) {
	double result = 1.;
	for (int i = 0; i < exponent; ++i) {
		result *= base;
	}

	return result;
}

double norm_vector(std::vector<double> x) {
	double norm_vector_x = 0;

	//for (int i = 0; i < x.size(); i++) 
	//	norm_vector_x += x[i] * x[i];

	return sqrt(x[0]*x[0]+x[1]*x[1]);
	// return sqrt(x[0] * x[0] + x[1] * x[1]);
}

std::vector<double> sub_coordinat(std::vector<double> x, std::vector<double> y) {
	int size = x.size();

	/*if (x.size() > y.size())
		size = y.size();

	for (int i = 0; i < size; i++)
		x[i] = x[i] - y[i];*/

	x[0] -= y[0];
	x[1] -= y[1];

	return x;
}


double particle::W(particle P) {
	std::vector<double> x_i = sub_coordinat(this->x_coordinat, P.x_coordinat);

	double h_a = (this->h_radius_kernel);
	double h_b = P.h_radius_kernel;

	double q_a = norm_vector(x_i) / h_a;
	double q_b = norm_vector(x_i) / h_b;

	double C_a = 1. / (pi*h_a*h_a);
	double C_b = 1. / (pi*h_b*h_b);

	double W_a = 0;
	double W_b = 0;

	if ((q_a <= 1) && (q_a >= 0))
		W_a = C_a     *    (15. / 7.)     *      (2. / 3. - q_a*q_a + 1. / 2.*q_a*q_a*q_a);

	if ((q_a > 1) && (q_a <= 2))
		W_a = C_a      *   (5. / 14.)      *      _pow(2 - q_a, 3);

	if ((q_b <= 1) && (q_b >= 0))
		W_b = C_b   *    (15. / 7.)     *      (2. / 3. - q_b*q_b + 1. / 2.*q_b*q_b*q_b);

	if ((q_b >1) && (q_b <= 2))
		W_b = C_b    *   (5. / 14.)      *      _pow(2 - q_b, 3);

	return 0.5*(W_a + W_b);
}



double particle::W_gradient(particle P) {
	std::vector<double> x_i = sub_coordinat(this->x_coordinat, P.x_coordinat);

	double h_a = (this->h_radius_kernel);
	double h_b = P.h_radius_kernel;

	double pi = 3.14159265358979;

	double q_a = norm_vector(x_i) / h_a;
	double q_b = norm_vector(x_i) / h_b;
	double sign = (x_i[0] > 0) - (x_i[0] < 0);

	double const_C_a = (2. / (3.*h_a))*(sign / h_a);
	double const_C_b = (2. / (3.*h_b))*(sign / h_b);

	double W_grad_a = 0;
	double W_grad_b = 0;

	if ((q_a <= 1) && (q_a >= 0))
		W_grad_a = const_C_a*(-3. * q_a + (9. / 4.)*_pow(q_a, 2));

	if ((q_a > 1) && (q_a <= 2))
		W_grad_a = const_C_a*((-3. / 4.)*    _pow(2 - q_a, 2));

	if ((q_b <= 1) && (q_b >= 0))
		W_grad_b = const_C_b*(-3. * q_b + (9. / 4.)*_pow(q_b, 2));

	if ((q_b > 1) && (q_b <= 2))
		W_grad_b = const_C_b*((-3. / 4.)*    _pow(2 - q_b, 2));

	return 0.5*(W_grad_a + W_grad_b);
}

double iron_particle::ro_particle(iron_particle *X, int N) {
	double ro = 0;
#pragma omp parallel for reduction(+:ro)
	for (int j = 0; j < N; j++)
		ro += this->m* (this->W(X[j]));
	return  ro;
}



double iron_particle::Viscosity(iron_particle X, double gamma) {//одномерная вязкость
	double r_ab = ((this->x_coordinat[0]) - (X.x_coordinat[0]));
	/*узнать точные условия на вязкость*/  /*близкие частицы*/               /*только проблемные частицы*/
	if ((this->v[0] - X.v[0])* r_ab>0)
		return 0;
	/*
	double a = 1.;   //никак в статье одномерный случай
	double b = 2.;

	double h_ab =         (this->h_radius_kernel + X.h_radius_kernel) / 2.;
	double V_ab =                                          (this->V)-(X.V);
	double ro_ab =                                      (X.ro+this->ro)/2.;
	double c_ab =                                      (X.cs+ this->cs)/2.;

	double nu =                                                   h_ab*0.1;
	double mu = (X.h_radius_kernel*V_ab*r_ab) / (_pow(r_ab, 2) + _pow(nu, 2));

	double viscocity=                  (-a*c_ab*mu + b*_pow(mu, 2)) / ro_ab ;//знак


	double a = 1.;
	double ro_ab = (X.ro + this->ro) / 2.;
	double V_ab = (this->V) - (X.V);

	double w_ab = V_ab* ((r_ab > 0) - (r_ab < 0));
	double v_sig = this->cs + X.cs - 3 * w_ab;

	double viscocity = (-a / 2.)*((v_sig*w_ab) / ro_ab);
	*/
	return 0;

}


