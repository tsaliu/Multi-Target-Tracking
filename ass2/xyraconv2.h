#pragma once
#include <armadillo>
#include <math.h>

void ra2xy4(arma::mat ra, arma::mat &xy, int radi, double sigv) {
	//arma::mat out_tmp;
	double a = ra(0, 0);
	double r = ra(0, 1);
	double b = exp(-pow(sigv, 2) / 2);
	double x = r * cos(a) + radi;
	double y = -r * sin(a) + radi;

	xy(0, 0) = x;
	xy(0, 1) = y;

	//xy = out_tmp;
}

void xy2ra4(arma::mat xy, arma::mat &ra, int radi, double sigv) {
	arma::mat out_tmp;
	double b = exp(-pow(sigv, 2) / 2);
	double x = xy(0, 0) - radi;
	double y = -(xy(0, 1) - radi);
	double in_r = sqrt(pow(x, 2) + pow(y, 2));
	double in_a = 0;

	in_a = atan2(y, x);

	ra(0, 0) = in_a;
	ra(0, 1) = in_r;

}


