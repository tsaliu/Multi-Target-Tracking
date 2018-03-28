#pragma once
#include <armadillo>
#include "xyraconv3.h"
double calcdis(arma::mat zk1, arma::mat hzk1k, arma::mat sk1, int i, int radi, double sigv) {
	/*arma::mat convxy(1, 2);
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

	arma::mat v = zk1 - hz_tmp;*/

	arma::mat v = zk1 - hzk1k(arma::span(i * 2, i * 2 + 1), arma::span(0, 3));


	arma::mat D = v.t()*arma::inv(sk1(arma::span(i * 2, i * 2 + 1), arma::span(0, 1)))*v;
	std::cout << i << std::endl;
	std::cout << zk1 << std::endl;
	std::cout << hzk1k(arma::span(i * 2, i * 2 + 1), arma::span(0, 3)) << std::endl;
	//std::cout << hz_tmp << std::endl;
	std::cout <<"sk1 "<< sk1(arma::span(i * 2, i * 2 + 1), arma::span(0, 1)) << std::endl;
	return arma::sum(arma::sum(D));
}