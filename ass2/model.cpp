#include "model.h"
#include <math.h>
#include <armadillo>
#include "xyraconvert.h"


model::model(){}
model::~model(){}

void model::init_pos(double initx, double inity) {
	init_x = initx;
	init_y = inity;
}

void model::getdata(int skip, int radi, int k, int st, double sigv, arma::mat &F, arma::mat &Q, double x, double y, double &dx, double &dy, double &da, double &dr) {
	int order = 2;
	//int T = st;
	
	if (k == 0) {
		pk = 0;
		ck = 0;
	}
	else {
		ck = k;
	}

	dk = ck - pk;
	//int T = dk * st;
	 
	int T = st + (st * skip);
	//std::cout << T << std::endl;
	pk = ck;

	arma::mat FF=arma::zeros(4, 4);
	arma::mat QQ=arma::zeros(4, 4);
	arma::mat G(4, 1);
	if (k % st == 0) {
		int n = k / st;

		if (n == 0) {
			p_x = init_x;
			p_y = init_y;
			//T = 0;
		}
		else {
			/*p_x = p_x;
			p_y = p_y;*/
		}

		for (int i = 0; i<(2 * order);) {
			FF(i, i) = 1;
			FF(i, i + 1) = T;
			FF(i + 1, i + 1) = 1;
			G(i, 0) = pow(T, 2) / 2;
			G(i + 1, 0) = T;

			i = i + 2;
		}
		//std::cout << "G   " << G << std::endl;
		QQ = G * pow(sigv, 2)*G.t();
	}
	dx = (x - p_x) / T;
	dy = (y - p_y) / T;

	if (skip != 0) {
		dx = 3;
		dy = 3;
		p_x = x - 3 * T;
		p_y = y - 3 * T;
	}

	arma::mat xy(1, 2);
	xy << x << y << arma::endr;
	arma::mat pxy(1, 2);
	pxy << p_x << p_y << arma::endr;
	/////////////////////////////
	//std::cout << skip << std::endl;
	//std::cout << "getdata xy" << x << "  " << y << std::endl;
	//std::cout << "getdata pxpy " << p_x << "  " << p_y << std::endl;
	////////////////////
	
	arma::mat ar(1, 2);
	xy2ra(xy, ar, radi, sigv);
	arma::mat par(1, 2);
	xy2ra(pxy, par, radi, sigv);

	da = (ar(0, 0) - par(0, 0)) / T;
	dr = (ar(0, 1) - par(0, 1)) / T;
	//////////////////
	//std::cout  <<"rrrrrrrrrrr" << ar << std::endl;
	//////////////////
	F = FF;
	Q = QQ;

	if (skip == 0) {
		p_x = x;
		p_y = y;
	}
	else {
		p_x = p_x;
		p_y = p_y;
	}
	


}
