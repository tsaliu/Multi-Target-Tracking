#pragma once
#include <armadillo>




class model{
	double init_x, init_y;
	double c_x, c_y, p_x, p_y;
	int ck, pk, dk;

public:
	model();
	~model();

	void init_pos(double, double);
	void getdata(int, int, int, double, arma::mat &, arma::mat &,
		double, double, double, double,
		double &, double &, double &, double &);
};

