#pragma once
#include <math.h>
#include <armadillo>

double likeicalc(arma::mat sk1, arma::mat zk1, arma::mat hzk1k, int i, int radi, double sigv) {
	/*arma::mat sk1_tmp = sk1(arma::span(i * 2, i * 2 + 1), arma::span(0, 1));

	arma::mat convxy(1, 2);
	arma::mat convar(1, 2);
	convar << zk1(0, 0) << zk1(1, 2) << arma::endr;
	ra2xy5(convar, convxy, radi, sigv);
	zk1(0, 0) = convxy(0, 0);
	zk1(1, 2) = convxy(0, 1);

	arma::mat hz_tmp = arma::zeros(2, 4);
	hz_tmp = hzk1k(arma::span(i * 2, i * 2 + 1), arma::span(0, 3));
	convar << hz_tmp(0, 0) << hz_tmp(1, 2) << arma::endr;
	ra2xy5(convar, convxy, radi, sigv);
	hz_tmp(0, 0) = convxy(0, 0);
	hz_tmp(1, 2) = convxy(0, 1);




	arma::mat hzk1k_tmp = hzk1k(arma::span(i * 2, i * 2 + 1), arma::span(0, 3));
	double siga = sqrt(sk1_tmp(0, 0));
	double sigr = sqrt(sk1_tmp(1, 1));
	double rho = sk1_tmp(0, 1) / (siga*sigr);

	double x_term = pow(zk1(0, 0) - hz_tmp(0, 0), 2) / sk1_tmp(0, 0);
	double y_term = pow(zk1(1, 2) - hz_tmp(1, 2), 2) / sk1_tmp(1, 1);
	double xy_term = 2 * rho*(zk1(0, 0) - hz_tmp(0, 0))*(zk1(1, 2) - hz_tmp(1, 2)) / (siga*sigr);
	double before_exp = -(x_term + y_term - xy_term) / (2 * (1 - pow(rho, 2)));
	double after_exp = exp(before_exp);
	double div = 2 * std::_Pi *siga*sigr*sqrt(1 - pow(rho, 2));
	double like = after_exp / div;
	*/

	arma::mat sk1_tmp = sk1(arma::span(i * 2, i * 2 + 1), arma::span(0, 1));
	arma::mat hzk1k_tmp = hzk1k(arma::span(i * 2, i * 2 + 1), arma::span(0, 3));
	double siga = sqrt(sk1_tmp(0, 0));
	double sigr = sqrt(sk1_tmp(1, 1));
	double rho = sk1_tmp(0, 1) / (siga*sigr);

	double x_term = pow(zk1(0, 0) - hzk1k_tmp(0, 0), 2) / sk1_tmp(0, 0);
	double y_term = pow(zk1(1, 2) - hzk1k_tmp(1, 2), 2) / sk1_tmp(1, 1);
	double xy_term = 2 * rho*(zk1(0, 0) - hzk1k_tmp(0, 0))*(zk1(1, 2) - hzk1k_tmp(1, 2)) / (siga*sigr);
	double before_exp = -(x_term + y_term - xy_term) / (2 * (1 - pow(rho, 2)));
	double after_exp = exp(before_exp);
	double div = 2 * std::_Pi *siga*sigr*sqrt(1 - pow(rho, 2));
	double like = after_exp / div;

	return like;
}