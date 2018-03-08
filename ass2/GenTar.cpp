#include "GenTar.h"
#include "randseed.h"


void GenTar::init_pos(double x, double y) {
	pos_x = x;
	pos_y = y;
}

void GenTar::init_vel(double vx, double vy){
	vel_x = vx;
	vel_y = vy;
}

void GenTar::noise(double nx, double ny) {
	noise_x = nx;
	noise_y = ny;
}

void GenTar::tot_length(int N) {
	length = N;
}

void GenTar::gen() {
	std::default_random_engine geng;
	geng.seed(std::chrono::system_clock::now().time_since_epoch().count());
	arma::mat data(2, length);
	arma::vec noise = arma::randn(2, 1);
	std::normal_distribution<double> distn1(0, sqrt(noise_x));
	std::normal_distribution<double> distn2(0, sqrt(noise_y));
	data(0, 0) = pos_x + distn1(geng);
	data(1, 0) = pos_y + distn2(geng);
	for (int i = 1; i < length; i++) {
		arma::vec noise = arma::randn(2, 1);
		data(0, i) = pos_x + i * vel_x + distn1(geng);
		data(1, i) = pos_y + i * vel_y + distn2(geng);
	}
	sim = data;
}

void GenTar::data(arma::mat &out) {
	out = sim.t();
}

void GenTar::print_sim() {
	std::cout << sim.t() << std::endl;
}

void GenTar::graph(cv::Mat &graphon, arma::mat datagraph, int k) {
	cv::Point target_at(datagraph(k, 0), datagraph(k, 1));
	if (k == 0) {
		cv::line(graphon, target_at, target_at, cv::Scalar(0, 0, 255), 2);
	}
	if (k > 0) {
		cv::Point target_pat(datagraph(k - 1, 0), datagraph(k - 1, 1));
		cv::line(graphon, target_at, target_pat, cv::Scalar(0, 0, 255), 2);
	}
}