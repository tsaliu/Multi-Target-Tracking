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

void model::getdata(int radi, int k, int st, double sigv, arma::mat &F, arma::mat &Q,
					double x, double y, double px, double py,
					double &dx, double &dy, double &da, double &dr) {
	int order = 2;
	int T = st;
	if (x == px) {
		T = 0;
	}
	

	arma::mat FF=arma::zeros(4, 4);
	arma::mat QQ=arma::zeros(4, 4);
	arma::mat G(4, 1);
	if (k % st == 0) {
		int n = k / st;


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
	dx = (x - px) / T;
	dy = (y - py) / T;



	arma::mat xy(1, 2);
	xy << x << y << arma::endr;
	arma::mat pxy(1, 2);
	pxy << px << py << arma::endr;
	/////////////////////////////
	//std::cout << "getdata xy" << x << "  " << y << std::endl;
	//std::cout << "getdata pxpy " << px << "  " << py << std::endl;
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


	


}
