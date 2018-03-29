#include "Sen.h"
#include "randseed3.h"
#include <random>

void Sen::getdata(arma::mat datain, double pdin) {
	Txydata = datain;
	double x = datain(0, 0)-500;
	double y =-(datain(0, 1)-500);
	double in_r = sqrt(pow(x, 2) + pow(y, 2));
	double in_a = 0;

	in_a = atan2(y, x);

	arma::mat data_tmp(1, 2);
	data_tmp(0, 0) = in_a;
	data_tmp(0, 1) = in_r;
	Tdata = data_tmp;

	pd = pdin;

}

void Sen::getdata2(arma::mat datain) {
	Txydata2 = datain;
	double x = datain(0, 0) - 500;
	double y = -(datain(0, 1) - 500);
	double in_r = sqrt(pow(x, 2) + pow(y, 2));
	double in_a = 0;

	in_a = atan2(y, x);

	arma::mat data_tmp(1, 2);
	data_tmp(0, 0) = in_a;
	data_tmp(0, 1) = in_r;
	Tdata2 = data_tmp;

}

void Sen::pos(int rx, int ry) {
	r_pos_x = rx;
	r_pos_y = ry;
}

void Sen::error(double er, double ea) {
	error_r = er;
	error_a = ea;
}

void Sen::samptime(double st) {
	sample_time = st;
}

void Sen::pdetect(double detection_prob) {
	pd = detection_prob;
}

void Sen::fa_density(double false_dens) {
	fa_d = false_dens;
}

void Sen::fa_gen(double nfa, int radi) {
	num_fa = nfa;

	std::default_random_engine genu;
	genu.seed(std::chrono::system_clock::now().time_since_epoch().count());

	arma::mat fa_tmp(2, num_fa);
	for (int i = 0; i < num_fa; i++) {
		std::uniform_real_distribution<double> distpi(-std::_Pi, std::_Pi);
		std::uniform_real_distribution<double> distr(0, radi);
		
		double a = distpi(genu);
		double r = distr(genu);

		fa_tmp(0, i) = a;
		fa_tmp(1, i) = r;
	}
	fa = fa_tmp.t();
	
}

void Sen::gen(double ea, double er) {

	std::default_random_engine genu;
	genu.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_real_distribution<double> distu(0.0, 1.0);

	pd_dis = distu(genu);
	pd_dis2 = distu(genu);
	int sencase = 0;
	bool test = Tdata.is_empty();
	bool test2 = Tdata2.is_empty();

	if (Tdata.is_empty() && Tdata2.is_empty()) {
		sencase = 1;
	}
	else if(Tdata.is_empty() && !Tdata2.is_empty()){
		sencase = 2;
	}
	else if (!Tdata.is_empty() && Tdata2.is_empty()) {
		sencase = 3;
	}
	else if (!Tdata.is_empty() && !Tdata2.is_empty()) {
		sencase = 4; /////// case 4 is not working !!!!! 3/28
	}

	switch (sencase) {
		case 1:
			sen_data = fa;
			tot_pts = num_fa;
			break;
		case 2:
			if (pd_dis > pd) {
				sen_data = fa;
				tot_pts = num_fa;
			}
			else if (pd_dis <= pd) {
				tot_pts = num_fa + 1;
				arma::mat sen_tmp = arma::resize(fa, tot_pts, 2);
				arma::mat Tsen_data_tmp(1, 2);

				std::normal_distribution<double> distn(0, sqrt(ea));
				std::normal_distribution<double> distn2(0, sqrt(er));
				arma::vec noise = arma::randn(2, 1);
				
				noise(0, 0) = distn(genu);
				noise(1, 0) = distn2(genu);
				

				for (int i = 0; i < tot_pts; i++) {
					if (i == 0) {
						sen_tmp(0, 0) = Tdata2(0, 0) + noise(0, 0);
						sen_tmp(0, 1) = Tdata2(0, 1) + noise(1, 0);
						Tsen_data_tmp(0, 0) = sen_tmp(0, 0);
						Tsen_data_tmp(0, 1) = sen_tmp(0, 1);
					}
					else {
						sen_tmp(i, 0) = fa(i - 1, 0);
						sen_tmp(i, 1) = fa(i - 1, 1);
					}
				}
				sen_data = sen_tmp;
				Tsen_data = Tsen_data_tmp;
			}
			break;
		case 3:
			if (pd_dis2 > pd) {
				sen_data = fa;
				tot_pts = num_fa;
			}
			else if (pd_dis2 <= pd) {
				tot_pts = num_fa + 1;
				arma::mat sen_tmp = arma::resize(fa, tot_pts, 2);
				arma::mat Tsen_data_tmp(1, 2);
				std::normal_distribution<double> distn(0, sqrt(ea));
				std::normal_distribution<double> distn2(0, sqrt(er));
				arma::vec noise = arma::randn(2, 1);
				
				noise(0, 0) = distn(genu);
				noise(1, 0) = distn2(genu);
				

				for (int i = 0; i < tot_pts; i++) {
					if (i == 0) {
						sen_tmp(0, 0) = Tdata(0, 0) + noise(0, 0);
						sen_tmp(0, 1) = Tdata(0, 1) + noise(1, 0);
						Tsen_data_tmp(0, 0) = sen_tmp(0, 0);
						Tsen_data_tmp(0, 1) = sen_tmp(0, 1);
					}
					else {
						sen_tmp(i, 0) = fa(i - 1, 0);
						sen_tmp(i, 1) = fa(i - 1, 1);
					}
				}
				sen_data = sen_tmp;
				Tsen_data = Tsen_data_tmp;
			}
			break;
		case 4:
			
			if (pd_dis > pd && pd_dis2 > pd) {
				sen_data = fa;
				tot_pts = num_fa;
			}
			else if (pd_dis > pd && pd_dis2 <= pd) {
				tot_pts = num_fa + 1;
				arma::mat sen_tmp = arma::resize(fa, tot_pts, 2);
				arma::mat Tsen_data_tmp(1, 2);
				std::normal_distribution<double> distn(0, sqrt(ea));
				std::normal_distribution<double> distn2(0, sqrt(er));
				arma::vec noise = arma::randn(2, 1);
				
				noise(0, 0) = distn(genu);
				noise(1, 0) = distn2(genu);
				

				for (int i = 0; i < tot_pts; i++) {
					if (i == 0) {
						sen_tmp(0, 0) = Tdata(0, 0) + noise(0, 0);
						sen_tmp(0, 1) = Tdata(0, 1) + noise(1, 0);
						Tsen_data_tmp(0, 0) = sen_tmp(0, 0);
						Tsen_data_tmp(0, 1) = sen_tmp(0, 1);
					}
					else {
						sen_tmp(i, 0) = fa(i - 1, 0);
						sen_tmp(i, 1) = fa(i - 1, 1);
					}
				}
				sen_data = sen_tmp;
				Tsen_data = Tsen_data_tmp;
			}
			else if (pd_dis <= pd && pd_dis2 > pd) {
				tot_pts = num_fa + 1;
				arma::mat sen_tmp = arma::resize(fa, tot_pts, 2);
				arma::mat Tsen_data_tmp(1, 2);
				std::normal_distribution<double> distn(0, sqrt(ea));
				std::normal_distribution<double> distn2(0, sqrt(er));
				arma::vec noise = arma::randn(2, 1);
				
				noise(0, 0) = distn(genu);
				noise(1, 0) = distn2(genu);
				

				for (int i = 0; i < tot_pts; i++) {
					if (i == 0) {
						sen_tmp(0, 0) = Tdata2(0, 0) + noise(0, 0);
						sen_tmp(0, 1) = Tdata2(0, 1) + noise(1, 0);
						Tsen_data_tmp(0, 0) = sen_tmp(0, 0);
						Tsen_data_tmp(0, 1) = sen_tmp(0, 1);
					}
					else {
						sen_tmp(i, 0) = fa(i - 1, 0);
						sen_tmp(i, 1) = fa(i - 1, 1);
					}
				}
				sen_data = sen_tmp;
				Tsen_data = Tsen_data_tmp;
			}
			else if (pd_dis <= pd && pd_dis2 <= pd) {
				tot_pts = num_fa + 2;
				arma::mat sen_tmp = arma::resize(fa, tot_pts, 2);
				arma::mat Tsen_data_tmp(1, 2);
				std::normal_distribution<double> distn(0, sqrt(ea));
				std::normal_distribution<double> distn2(0, sqrt(er));
				arma::vec noise = arma::randn(2, 1);
				arma::vec noise2 = arma::randn(2, 1);
				
				noise(0, 0) = distn(genu);
				noise(1, 0) = distn2(genu);
				noise2(0, 0) = distn(genu);
				noise2(1, 0) = distn2(genu);

				for (int i = 0; i < tot_pts; i++) {
					if (i == 0) {
						sen_tmp(i, 0) = Tdata(0, 0) + noise(0, 0);
						sen_tmp(i, 1) = Tdata(0, 1) + noise(1, 0);
						Tsen_data_tmp(0, 0) = sen_tmp(0, 0);
						Tsen_data_tmp(0, 1) = sen_tmp(0, 1);
					}
					if (i == 1) {
						sen_tmp(i, 0) = Tdata2(0, 0) + noise2(0, 0);
						sen_tmp(i, 1) = Tdata2(0, 1) + noise2(1, 0);
						Tsen_data_tmp(0, 0) = sen_tmp(0, 0);
						Tsen_data_tmp(0, 1) = sen_tmp(0, 1);
					}
					if (i > 1) {
						sen_tmp(i, 0) = fa(i - 2, 0);
						sen_tmp(i, 1) = fa(i - 2, 1);
					}
					//std::cout << Tdata2 << std::endl;
				}
				//std::cout << sen_data << std::endl;
				sen_data = sen_tmp;
				Tsen_data = Tsen_data_tmp;
			}
			break;
	}
}

void Sen::graph(cv::Mat &graphon, int radi, int k, int st) {
	if (k%st == 0) {
		for (int i = 0; i < tot_pts; i++) {
			double a = sen_data(i, 0);
			double r = sen_data(i, 1);
			double x = r * cos(a) + radi;
			double y = -r * sin(a) + radi;

			cv::Point senpoint(x, y);
			cv::line(graphon, senpoint, senpoint, cv::Scalar(0, 255, 0), 2);
		}
	}
}

void Sen::meas_data(arma::mat &outd, int k, int st) {
	
	if (k%st == 0) {
		for (int i = 0; i < tot_pts; i++) {
			outd(i, k / st * 2 + 0) = sen_data(i, 0);
			outd(i, k / st * 2 + 1) = sen_data(i, 1);
		}

		//n++;
	}
	
}

void Sen::truth(arma::mat &outt, int k) {
	outt(k, 0) = Tdata(0, 0);
	outt(k, 1) = Tdata(0, 1);
}

void Sen::truthxy(arma::mat &outxy, int k) {
	outxy(k, 0) = Txydata(0, 0);
	outxy(k, 1) = Txydata(0, 1);
}

void Sen::truth2(arma::mat &outt2, int k) {
	outt2(k, 0) = Tdata2(0, 0);
	outt2(k, 1) = Tdata2(0, 1);
}

void Sen::truthxy2(arma::mat &outxy2, int k) {
	outxy2(k, 0) = Txydata2(0, 0);
	outxy2(k, 1) = Txydata2(0, 1);
}

void Sen::graphmeas(cv::Mat &graphon, arma::mat datagraph, int k, int st, int radi) {
	
	if (k % st == 0) {
		int n = k / st;

		double a = datagraph(0, n * 2 + 0);
		double r = datagraph(0, n * 2 + 1);
		double x = r * cos(a) + radi;
		double y = -r * sin(a) + radi;

		cv::Point target_at(x, y);
		if (n == 0) {
			cv::line(graphon, target_at, target_at, cv::Scalar(255, 0, 0), 2);
		}
		if (n > 0) {
			double a2 = datagraph(0, (n - 1) * 2 + 0);
			double r2 = datagraph(0, (n - 1) * 2 + 1);
			double x2 = r2 * cos(a2) + radi;
			double y2 = -r2 * sin(a2) + radi;
			if (sqrt(pow(Txydata(0, 0) - x2, 2) + pow(Txydata(0, 1) - y2, 2)) < 100) {
				cv::Point target_pat(x2, y2);
				cv::line(graphon, target_pat, target_pat, cv::Scalar(255, 0, 0), 2);
			}
		}
		n++;
	}
}