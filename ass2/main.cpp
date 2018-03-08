#include <iostream>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <armadillo>
#include <string>
#include <Windows.h>
#include <random>
#include <cmath>

#include "initpara.h"
#include "GenTar.h"
#include "senpara.h"
#include "Sen.h"

#include "EKF.h"
#include "runs.h"

cv::Mat frame;
int maxdetect = 50;
arma::mat target_data;
arma::mat target_data2;
arma::mat truth = arma::zeros(len, 2);
arma::mat truthxy = arma::zeros(len, 2);
arma::mat truth2 = arma::zeros(len2, 2);
arma::mat truthxy2 = arma::zeros(len2, 2);

arma::mat P = arma::zeros(4, 4 * maxdetect);
arma::mat xk1k1 = arma::zeros(4, 4);
arma::mat hzk1k(2, 4);
arma::mat sk1(2, 2);
arma::mat t_id = arma::zeros(maxdetect, 1);
arma::mat chxk1k1 = arma::zeros(4 * maxdetect, 4);
arma::mat phxk1k1 = arma::zeros(4 * maxdetect, 4);


int main(int argc, char *argv[]) {
	srand(time(NULL));
	cv::Mat initbg(bgwidth, bgheight, CV_8UC3, cv::Scalar(0, 0, 0));
	cv::circle(initbg, center, radius, cv::Scalar(255, 255, 255), 1);
	cv::putText(initbg, "R", Rat, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 255, 255), 2);
	cv::namedWindow("Output", cv::WINDOW_AUTOSIZE);

	GenTar target;
	target.init_pos(initx, inity);
	target.init_vel(initvx, initvy);
	target.noise(nx, ny);
	target.tot_length(len);
	target.gen();
	target.data(target_data);

	GenTar target2;
	target2.init_pos(initx2, inity2);
	target2.init_vel(initvx2, initvy2);
	target2.noise(nx2, ny2);
	target2.tot_length(len2);
	target2.gen();
	target2.data(target_data2);

	initbg.copyTo(frame);

	Sen meas;
	meas.pos(rx, ry);
	meas.error(er, ea);
	meas.samptime(st);
	meas.pdetect(pd);
	meas.fa_density(fad);

	cv::Mat graph;

	model kfmodel;
	int nruns = 1;

	runs result;
	std::default_random_engine genp;
	std::poisson_distribution<> poiss(numfam);
	for (int i = 0; i < nruns; i++) {

		for (int k = 0; k < lent; k++) {
			numfa = poiss(genp);
			if (numfa == 0) {
				//numfa = 1;
			}
			int k1 = k - initt;
			int k2 = k - initt2;

			arma::mat meas_data = arma::zeros(numfa * 3, lent * 2);

			if (k >= initt && k < endt) {
				target.graph(frame, target_data, k1);
				meas.getdata(target_data(k1, arma::span(0, 1)));
				meas.truth(truth, k1);				//target with no meas noise ar
				meas.truthxy(truthxy, k1);
			}
			if (k >= initt2 && k < endt2) {
				target2.graph(frame, target_data2, k2);
				meas.getdata2(target_data2(k2, arma::span(0, 1)));
				meas.truth2(truth2, k2);				//target with no meas noise ar
				meas.truthxy2(truthxy2, k2);
			}
			frame.copyTo(graph);

			meas.fa_gen(numfa, radius);
			meas.gen(ea, er);
			//meas.graph(frame, radius, k, st);	//all meas, truth with noise + fa
			meas.meas_data(meas_data, k, st);
			
			arma::uvec meas_at = arma::find(meas_data(arma::span(0, numfa * 3 - 1), (k / st) * 2 + 0) != 0);
			int meas_size = arma::max(meas_at) + 1;

			std::cout << meas_size << std::endl;
			std::cout << meas_at << std::endl;
			std::cout << meas_data(meas_size - 1, (k / st) * 2) << std::endl;
			if (k == 0) {
				
			}





			/*
			EKF ekf1, ekf2;
			if (k >= initt && k < endt) {
				ekf1.init(k1, truth, radius);
				ekf1.start(k1, st, radius, p00, R, P, xk1k1, hzk1k, sk1);
				p00 = P;
				ekf1.asso(k1+initt, st, meas_data, numfa, hzk1k, sk1, radius);
				ekf1.graph_est(frame, k1, st);
			}
			if (k >= initt2 && k < endt2) {
				ekf2.init(k2, truth2, radius);
				ekf2.start(k2, st, radius, p002, R2, P2, xk1k12, hzk1k2, sk12);
				p002 = P2;
				ekf2.asso(k2+initt2, st, meas_data, numfa, hzk1k2, sk12, radius);
				ekf2.graph_est(frame, k2, st);
			}
			
			*/
			
			//std::cout << xk1k12 << std::endl;
		}
		
	}
	//result.graphavg(graph, avg_xy, len, st);
	cv::imshow("Output", frame);
	cv::waitKey(1);

	std::cout << "Done" << std::endl;
	while (!GetAsyncKeyState(VK_ESCAPE)) {}
}