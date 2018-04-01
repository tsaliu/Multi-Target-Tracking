#pragma once
#include <armadillo>
#include <math.h>


class IPDA{
	int maxdetect = 50;
	double pg = 0.95;
	double ps = 0.9;
	int deg = 2;
	double gamma, lambda, pd;
	double init_quality = 0.2;
	double com_th = 0.95;
	double term_th = 0.05;

	//since likelihood is scalar
	arma::mat dk;
	double dk_mid, dk_tmp;


public:
	void getpara(double, double);
	void ipda(arma::mat,
			arma::mat::col_iterator,
			arma::mat, arma::mat, arma::mat, int, double,
			arma::mat &, arma::mat &,
			int, double, arma::mat &);
	void nn(arma::mat,
		arma::mat::col_iterator,
		arma::mat, arma::mat, arma::mat, int, double,
		arma::mat &, arma::mat &,
		int, double, arma::mat &);

};

