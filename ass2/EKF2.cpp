#include "EKF2.h"
#include "xyraconv.h"



void EKF2::kf(int k, double st, int radi, arma::mat meas_data, arma::mat &P, arma::mat R,
	arma::mat phxk1k1, arma::mat &chxk1k1, arma::mat &hzk1k, arma::mat &sk1,
	arma::mat &t_id) {

	int n = k / st * 2;
	std::default_random_engine gene;
	gene.seed(std::chrono::system_clock::now().time_since_epoch().count());

	double b = exp(-pow(sigv, 2) / 2);

	arma::mat::col_iterator cita = meas_data.begin_col(n); 
	arma::mat::col_iterator cita_end = meas_data.end_col(n);  

	arma::mat::col_iterator citr = meas_data.begin_col(n + 1);
	arma::mat::col_iterator citr_end = meas_data.end_col(n + 1);

}

void EKF2::ipda() {

}
