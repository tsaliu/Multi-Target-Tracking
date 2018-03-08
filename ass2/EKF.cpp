#include "EKF.h"




void EKF::kf(int k, double st, int radi, arma::mat meas_data, arma::mat &P,
	arma::mat phxk1k1, arma::mat &chxk1k1, arma::mat &hzk1k, arma::mat &sk1,
	arma::mat &t_id) {
	arma::mat p00 = arma::zeros(4, 4);
	p00(0, 0) = 100;
	p00(2, 2) = 100;
	p00(1, 1) = 20;
	p00(3, 3) = 20;
	arma::mat R = arma::zeros(2, 2);
	R(0, 0) = 0.05;
	R(1, 1) = 1;

	int n = k / st * 2;
	std::default_random_engine gene;
	gene.seed(std::chrono::system_clock::now().time_since_epoch().count());

	double b = exp(-pow(sigv, 2) / 2);

	arma::mat::col_iterator cita = meas_data.begin_col(n); 
	arma::mat::col_iterator cita_end = meas_data.end_col(n);  

	arma::mat::col_iterator citr = meas_data.begin_col(n + 1);
	arma::mat::col_iterator citr_end = meas_data.end_col(n + 1);

	arma::mat cmeas_tmp(1, 2);
	arma::mat cstate_tmp(1, 2);

	arma::uvec id_check_m = arma::find(t_id != 0); //means just started, nothing ran yet
	bool id_check_start = id_check_m.is_empty();
	if (id_check_start) {
		int i = 1;
		for (; cita != cita_end; cita++) {
			t_id(i - 1, 1) = i;
			
			cmeas_tmp << (*cita) << (*citr) << arma::endr;
			ra2xy(cmeas_tmp, cstate_tmp, radi, sigv);
			chxk1k1(k / st * 4 + 0, 0) = cstate_tmp(0, 0);
			chxk1k1(k / st * 4 + 2, 2) = cstate_tmp(0, 1);

			i++;
			citr++;
		}
	}
	

}

void EKF::ipda() {

}
