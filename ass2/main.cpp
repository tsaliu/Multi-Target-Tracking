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

arma::mat meas_data = arma::zeros(maxdetect, lent * 2);

arma::mat P = arma::zeros(4 * maxdetect, 4);
//arma::mat xk1k1 = arma::zeros(4, 4);
arma::mat hzk1k = arma::zeros(2 * maxdetect, 4);
arma::mat sk1 = arma::zeros(2 * maxdetect, 2);
arma::mat t_id = arma::zeros(maxdetect, lent);
arma::mat chxk1k1 = arma::zeros(4 * maxdetect, 4);
arma::mat phxk1k1 = arma::zeros(4 * maxdetect, 4);
arma::mat q = arma::zeros(maxdetect, lent); //quality, need change size for all length

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

	EKF ekf;

	model kfmodel;
	int nruns = 1;

	runs result;
	std::default_random_engine genp;
	std::poisson_distribution<> poiss(numfam);
	arma::mat run_tar1 = arma::zeros(endt - initt - 1, 7);
	arma::mat run_tar2 = arma::zeros(endt2 - initt2, 7);
	arma::mat run_tar1_count = arma::zeros(endt - initt - 1, 1);
	arma::mat run_tar2_count = arma::zeros(endt2 - initt2, 1);
	for (int i = 0; i < nruns; i++) {
		arma::mat save_com_all;
		arma::mat tar1;
		arma::mat tar2;
		std::cout << "run " << i + 1 << std::endl;
		for (int k = 0; k < lent; k++) {
			numfa = poiss(genp);
			if (numfa == 0) {
				//numfa = 1;
			}
			int k1 = k - initt;
			int k2 = k - initt2;

			

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

			//std::cout << meas_size << std::endl;
			//std::cout << meas_at << std::endl;
			//std::cout << meas_data(meas_size - 1, (k / st) * 2) << std::endl;
			
			arma::mat save_com;
			//may need to clear hzk1k, sk1
			ekf.kf(frame, k, st, radius, meas_data, P, phxk1k1, chxk1k1, hzk1k, sk1, t_id, fad, pd, q, save_com);
			

			
			if (k >= initt && k < initt2) {
				arma::mat sum = arma::sum(save_com) / save_com.n_rows;
				std::cout << "sum " << sum << std::endl;
				if (!sum.is_empty()) {
					tar1.insert_rows(tar1.n_rows, sum);
				}
				else {
					std::cout << "sum empt" << std::endl;
				}
				
			}
			else if (k>=initt2 && k< endt) {
				int num1 = 0;
				int num2 = 0;
				arma::mat sum1 = arma::zeros(1, 7);
				arma::mat sum2 = arma::zeros(1, 7);
				for (int ii = 0; ii < save_com.n_rows; ii++) {
					double dis1 = sqrt(pow((truth(k1, 0) - save_com(ii, 5)), 2) + pow((truth(k1, 1) - save_com(ii, 6)), 2));
					double dis2 = sqrt(pow((truth2(k2, 0) - save_com(ii, 5)), 2) + pow((truth2(k2, 1) - save_com(ii, 6)), 2));
					if (dis1 < dis2) {
						num1++;
						sum1 = sum1 + save_com.row(ii);
					}
					else if (dis1 > dis2) {
						num2++;
						sum2 = sum2 + save_com.row(ii);
					}
				}
				sum1 = sum1 / num1;
				sum2 = sum2 / num2;
				tar1.insert_rows(tar1.n_rows, sum1);
				tar2.insert_rows(tar2.n_rows, sum2);
			}
			else if (k >= endt && k <= endt2) {
				arma::mat sum = arma::sum(save_com) / save_com.n_rows;
				//std::cout << "sum " << sum << std::endl;
				tar2.insert_rows(tar2.n_rows, sum);
			}
			
			//save_com_all.insert_rows(save_com_all.n_rows, save_com);
			//std::cout << "save " << save_com_all << std::endl;

			
			



		}
		P.resize(4 * maxdetect, 4);
		P.fill(0);
		phxk1k1.fill(0);
		chxk1k1.fill(0);
		hzk1k.fill(0);
		sk1.fill(0);
		t_id.fill(0);
		q.fill(0);
		tar1.elem(arma::find_nonfinite(tar1)).zeros();
		tar2.elem(arma::find_nonfinite(tar2)).zeros();
		run_tar1.elem(arma::find_nonfinite(run_tar1)).zeros();
		run_tar2.elem(arma::find_nonfinite(run_tar2)).zeros();
		for (int ii = 0; ii < tar1.n_rows; ii++) {
			if (arma::sum(tar1.row(ii)) != 0) {
				run_tar1_count(ii, 0) += 1;
			}
		}
		for (int ii = 0; ii < tar2.n_rows; ii++) {
			if (arma::sum(tar2.row(ii)) != 0) {
				run_tar2_count(ii, 0) += 1;
			}
		}
		run_tar1 = (run_tar1 + tar1);
		run_tar2 = (run_tar2 + tar2);
		
		run_tar1.elem(arma::find_nonfinite(run_tar1)).zeros();
		run_tar2.elem(arma::find_nonfinite(run_tar2)).zeros();
		std::cout << "tar1 " << tar1 << std::endl;
		std::cout << "tar2 " << tar2 << std::endl;
	}
	for (int ii = 0; ii < run_tar1.n_rows; ii++) {
		run_tar1.row(ii) = run_tar1.row(ii) / (run_tar1_count(ii, 0));
	}
	for (int ii = 0; ii < run_tar2.n_rows; ii++) {
		run_tar2.row(ii) = run_tar2.row(ii) / (run_tar2_count(ii, 0));
	}
	run_tar1.elem(arma::find_nonfinite(run_tar1)).zeros();
	run_tar2.elem(arma::find_nonfinite(run_tar2)).zeros();
	//result.graphavg(graph, avg_xy, len, st);
	cv::imshow("Output", frame);
	cv::waitKey(0);


	std::cout << run_tar1 << std::endl;
	std::cout << run_tar2 << std::endl;
	//meas_data.save("meas.txt", arma::arma_ascii);
	std::cout << "Done" << std::endl;
	while (!GetAsyncKeyState(VK_ESCAPE)) {}
}