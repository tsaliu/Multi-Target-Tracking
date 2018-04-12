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

#include "xyraconv4.h"

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
	//cv::namedWindow("Output", cv::WINDOW_AUTOSIZE);

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
	frame.copyTo(graph);

	EKF ekf;

	model kfmodel;

	runs result;
	std::default_random_engine genp;
	std::poisson_distribution<> poiss(numfam);
	arma::mat run_tar1 = arma::zeros(endt - initt, 3);
	arma::mat run_tar2 = arma::zeros(endt2 - initt2, 3);
	arma::mat run_tarxy1 = arma::zeros(endt - initt, 3);
	arma::mat run_tarxy2 = arma::zeros(endt2 - initt2, 3);
	arma::mat run_tarxyv1 = arma::zeros(endt - initt - 1, 3);
	arma::mat run_tarxyv2 = arma::zeros(endt2 - initt2 - 1, 3);
	arma::mat run_tar1_count = arma::zeros(endt - initt, 3);
	arma::mat run_tar2_count = arma::zeros(endt2 - initt2, 3);
	arma::mat run_tar1v_count = arma::zeros(endt - initt - 1, 3);
	arma::mat run_tar2v_count = arma::zeros(endt2 - initt2 - 1, 3);
	arma::mat run_tar1_avg = arma::zeros(endt - initt, 3);
	arma::mat run_tar2_avg = arma::zeros(endt2 - initt2, 3);
	arma::mat save_com_all;
	double avg_late1 = 0;
	double avg_late2 = 0;
	double avg_ft_rate = 0;
	for (int i = 0; i < nruns; i++) {
		
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
				if (i == 0) {
					target.graph(graph, target_data, k1);
				}
				meas.getdata(target_data(k1, arma::span(0, 1)), pd);
				meas.truth(truth, k1);				//target with no meas noise ar
				meas.truthxy(truthxy, k1);
			}
			if (k >= initt2 && k < endt2) {
				target2.graph(frame, target_data2, k2);
				if (i == 0) {
					target2.graph(graph, target_data2, k2);
				}
				meas.getdata2(target_data2(k2, arma::span(0, 1)));
				meas.truth2(truth2, k2);				//target with no meas noise ar
				meas.truthxy2(truthxy2, k2);
			}
			
			//std::cout << "HERE" << std::endl;

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

			ekf.kf(frame, k, st, radius, meas_data, P, phxk1k1, chxk1k1, hzk1k, sk1, t_id, fad, pd, q, save_com, ipda_mode);
			//std::cout << chxk1k1 << std::endl;
			if (k >= initt2) {
				//while (!GetAsyncKeyState(VK_SPACE)) {}
			}
			
			if (k >= initt && k < initt2) {
				arma::mat sum = arma::sum(save_com) / save_com.n_rows;
				//std::cout << "sum " << sum << std::endl;
				if (!sum.is_empty()) {
					tar1.insert_rows(tar1.n_rows, sum);
				}
				else {
					//std::cout << "sum empt" << std::endl;
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
					if (dis1 <= 50 || dis2 <= 50) {
						if (dis1 < dis2) {
							num1++;
							sum1 = sum1 + save_com.row(ii);
						}
						else if (dis1 > dis2) {
							num2++;
							sum2 = sum2 + save_com.row(ii);
						}
					}
				}
				sum1 = sum1 / num1;
				sum2 = sum2 / num2;
				tar1.insert_rows(tar1.n_rows, sum1);
				tar2.insert_rows(tar2.n_rows, sum2);
			}
			else if (k >= endt && k < endt2) {
				arma::mat sum = arma::sum(save_com) / save_com.n_rows;
				//std::cout << "sum " << sum << std::endl;
				tar2.insert_rows(tar2.n_rows, sum);
			}
			
			save_com_all.insert_rows(save_com_all.n_rows, save_com);
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

		arma::mat truthk = truth;
		arma::mat truthk2 = truth2;
		arma::mat truthkxy = truthxy;
		arma::mat truthkxy2 = truthxy2;
		arma::mat truthkxyv = arma::zeros(endt - initt - 1, 3);
		arma::mat truthkxyv2 = arma::zeros(endt2 - initt2 - 1, 3);
		truthkxyv.col(1).fill(initvx);
		truthkxyv.col(2).fill(initvy);
		truthkxyv2.col(1).fill(initvx2);
		truthkxyv2.col(2).fill(initvy2);
		arma::mat ks = arma::zeros(lent, 1);
		for (int ii = 0; ii < ks.n_rows; ii++) {
			ks(ii, 0) = ii + 1;
		}
		truthk.insert_cols(0, ks(arma::span(initt - 1, initt - 1 + truthk.n_rows - 1), 0));
		truthkxy.insert_cols(0, ks(arma::span(initt - 1, initt - 1 + truthkxy.n_rows - 1), 0));
		truthk2.insert_cols(0, ks(arma::span(initt2 - 1, initt2 - 1 + truthk2.n_rows - 1), 0));
		truthkxy2.insert_cols(0, ks(arma::span(initt2 - 1, initt2 - 1 + truthkxy2.n_rows - 1), 0));
		truthkxyv.col(0) = truthk(arma::span(1, truthk.n_rows - 1), 0);
		truthkxyv2.col(0) = truthk2(arma::span(1, truthk2.n_rows - 1), 0);
		double min_r = arma::min(truthk.col(2));
		double max_r = arma::max(truthk.col(2));
		double min_a = arma::min(truthk.col(1));
		double max_a = arma::max(truthk.col(1));
		double min_r2 = arma::min(truthk2.col(2));
		double max_r2 = arma::max(truthk2.col(2));
		double min_a2 = arma::min(truthk2.col(1));
		double max_a2 = arma::max(truthk2.col(1));
		

		int num_com = 0;
		int num_ft = 0;
		arma::mat tar1_avg = arma::zeros(endt - initt, 3);
		arma::mat tar2_avg = arma::zeros(endt2 - initt2, 3);
		arma::mat tar1_avgxy = arma::zeros(endt - initt, 3);
		arma::mat tar2_avgxy = arma::zeros(endt2 - initt2, 3);
		arma::mat tar1_avgxyv = arma::zeros(endt - initt - 1, 3);
		arma::mat tar2_avgxyv = arma::zeros(endt2 - initt2 - 1, 3);
		arma::mat tar1_ct = arma::zeros(endt - initt, 1);
		arma::mat tar2_ct = arma::zeros(endt2 - initt2, 1);
		for (int ii = 0; ii < save_com_all.n_rows; ii++) {
			arma::uvec find_same_id = arma::find(save_com_all.col(1) == save_com_all(ii, 1));
			arma::mat extract_track;
			//for every unique id
			if (save_com_all(find_same_id(0), 1) != 0) {
				//std::cout << find_same_id << std::endl;
				if (!find_same_id.is_empty()) {
					for (int iii = 0; iii < find_same_id.n_rows; iii++) {

						extract_track.insert_rows(extract_track.n_rows, save_com_all.row(find_same_id(iii)));

						save_com_all.row(find_same_id(iii)).fill(0);
					}
				}
				//std::cout << extract_track << std::endl;

				if ((arma::min(extract_track.col(5)) < (min_a) * 0.8 && arma::min(extract_track.col(5)) > (max_a2) * 0.8) ||
					(arma::max(extract_track.col(5)) > (max_a2)*0.8 && arma::max(extract_track.col(5)) < (min_a)*0.8) ||
					(arma::min(extract_track.col(6)) < (min_r)*0.8 || arma::max(extract_track.col(6)) > (max_r*1.10)) ||
					(arma::min(extract_track.col(6)) < (min_r2)*0.8 || arma::max(extract_track.col(6)) > (max_r2*1.10))) {
					std::cout << "ft id  " << extract_track(0,1)<< std::endl;
					num_ft++;
				}
				else {
					bool choose_tar = 0;
					if (arma::mean(extract_track.col(5)) > 0) {
						choose_tar = 1;
						std::cout << "tar 1 track " << extract_track(0, 1) << std::endl;
						for (int iii = 0; iii < extract_track.n_rows; iii++) {
							arma::uvec match = arma::find(truthk.col(0) == extract_track(iii, 0));
							if (!match.is_empty()) {
								tar1_ct(match(0), 0) = tar1_ct(match(0), 0) + 1;
								tar1_avg(match(0), 1) += extract_track(iii, 5); 
								tar1_avg(match(0), 2) += extract_track(iii, 6); 

								arma::mat convar1(1, 2);
								arma::mat convxy1(1, 2);
								convar1 << extract_track(iii, 5) << extract_track(iii, 6) << arma::endr;
								ra2xy6(convar1, convxy1, radius);
								tar1_avgxy(match(0), 1) += convxy1(0, 0);
								tar1_avgxy(match(0), 2) += convxy1(0, 1);
								
							}
						}
					}
					else if (arma::mean(extract_track.col(5)) < 0) {
						choose_tar = 2;
						std::cout << "tar 2 track " << extract_track(0, 1) << std::endl;
						for (int iii = 0; iii < extract_track.n_rows; iii++) {
							arma::uvec match = arma::find(truthk2.col(0) == extract_track(iii, 0));
							if (!match.is_empty()) {
								tar2_ct(match(0), 0) = tar2_ct(match(0), 0) + 1;
								tar2_avg(match(0), 1) += extract_track(iii, 5);
								tar2_avg(match(0), 2) += extract_track(iii, 6);

								arma::mat convar1(1, 2);
								arma::mat convxy1(1, 2);
								convar1 << extract_track(iii, 5) << extract_track(iii, 6) << arma::endr;
								ra2xy6(convar1, convxy1, radius);
								tar2_avgxy(match(0), 1) += convxy1(0, 0);
								tar2_avgxy(match(0), 2) += convxy1(0, 1);

							}
						}
					}
				}
				
				
				num_com++;
			}
			
			//if()

		}
		int late_count1 = 0;
		bool late_first1 = true;
		for (int ii = 0; ii < tar1_avg.n_rows; ii++) {
			if (tar1_ct(ii, 0) != 0) {
				tar1_avg(ii, 1) = tar1_avg(ii, 1) / tar1_ct(ii, 0);
				tar1_avg(ii, 2) = tar1_avg(ii, 2) / tar1_ct(ii, 0);
				tar1_avgxy(ii, 1) = tar1_avgxy(ii, 1) / tar1_ct(ii, 0);
				tar1_avgxy(ii, 2) = tar1_avgxy(ii, 2) / tar1_ct(ii, 0);
			}
			if (tar1_ct(ii, 0) == 0 && late_first1) {
				late_count1++;
			}
			else {
				late_first1 = false;
			}
		}
		int late_count2 = 0;
		bool late_first2 = true;
		for (int ii = 0; ii < tar1_avg.n_rows; ii++) {
			if (tar2_ct(ii, 0) != 0) {
				tar2_avg(ii, 1) = tar2_avg(ii, 1) / tar2_ct(ii, 0);
				tar2_avg(ii, 2) = tar2_avg(ii, 2) / tar2_ct(ii, 0);
				tar2_avgxy(ii, 1) = tar2_avgxy(ii, 1) / tar2_ct(ii, 0);
				tar2_avgxy(ii, 2) = tar2_avgxy(ii, 2) / tar2_ct(ii, 0);
			}
			if (tar2_ct(ii, 0) == 0 && late_first2) {
				late_count2++;
			}
			else {
				late_first2 = false;
			}
		}
		tar1_avgxyv = arma::diff(tar1_avgxy);
		tar2_avgxyv = arma::diff(tar2_avgxy);

		tar1_avgxyv.col(0) = truthk(arma::span(1, truthk.n_rows - 1), 0);
		tar2_avgxyv.col(0) = truthk2(arma::span(1, truthk2.n_rows - 1), 0);
		tar1_avg.col(0) = truthk.col(0);
		tar2_avg.col(0) = truthk2.col(0);
		tar1_avgxy.col(0) = truthk.col(0);
		tar2_avgxy.col(0) = truthk2.col(0);
		run_tar1.col(0) = truthk.col(0);
		run_tar2.col(0) = truthk2.col(0);
		run_tarxy1.col(0) = truthk.col(0);
		run_tarxy2.col(0) = truthk2.col(0);
		run_tarxyv1.col(0) = truthk(arma::span(1, truthk.n_rows - 1), 0);
		run_tarxyv2.col(0) = truthk2(arma::span(1, truthk2.n_rows - 1), 0);
		run_tar1_avg.col(0) = truthk.col(0);
		run_tar2_avg.col(0) = truthk2.col(0);
		std::cout << tar1_avg << std::endl;
		std::cout << tar2_avg << std::endl;
		std::cout << tar1_avgxy << std::endl;
		std::cout << tar2_avgxy << std::endl;
		std::cout << tar1_avgxyv << std::endl;
		std::cout << tar2_avgxyv << std::endl;
		

		double ft_rate = (double)num_ft / (double)num_com;
		arma::uvec find_late1 = arma::find(tar1_avg.col(1) == 0);
		double late1 = late_count1;
		arma::uvec find_late2 = arma::find(tar2_avg.col(1) == 0);
		double late2 = late_count2;
		std::cout << "false track rate " << ft_rate << std::endl;
		std::cout << "late " << late1 << "  " << late2 << std::endl;

		for (int ii = 0; ii < tar1_avg.n_rows; ii++) {
			if (tar1_avg(ii, 1) != 0) {
				run_tar1_count(ii, 0) += 1;
				run_tar1(ii, arma::span(1, 2)) += (tar1_avg(ii, arma::span(1, 2)) - truthk(ii, arma::span(1, 2))) % (tar1_avg(ii, arma::span(1, 2)) - truthk(ii, arma::span(1, 2)));
				run_tarxy1(ii, arma::span(1, 2)) += (tar1_avgxy(ii, arma::span(1, 2)) - truthkxy(ii, arma::span(1, 2))) % (tar1_avgxy(ii, arma::span(1, 2)) - truthkxy(ii, arma::span(1, 2)));
			}
		}
		for (int ii = 0; ii < tar2_avg.n_rows; ii++) {
			if (tar2_avg(ii, 1) != 0) {
				run_tar2_count(ii, 0) += 1;
				run_tar2(ii, arma::span(1, 2)) += (tar2_avg(ii, arma::span(1, 2)) - truthk2(ii, arma::span(1, 2))) % (tar2_avg(ii, arma::span(1, 2)) - truthk2(ii, arma::span(1, 2)));
				run_tarxy2(ii, arma::span(1, 2)) += (tar2_avgxy(ii, arma::span(1, 2)) - truthkxy2(ii, arma::span(1, 2))) % (tar2_avgxy(ii, arma::span(1, 2)) - truthkxy2(ii, arma::span(1, 2)));
			}
		}
		run_tar1_avg += tar1_avgxy;
		run_tar2_avg += tar2_avgxy;
		//run_tar1_avg.elem(arma::find_nonfinite(run_tar1_avg)).fill(0);
		//run_tar2_avg.elem(arma::find_nonfinite(run_tar2_avg)).fill(0);
		std::cout << run_tar1 << std::endl;

		for (int ii = 0; ii < tar1_avgxyv.n_rows; ii++) {
			if (tar1_avgxyv(ii, 1) != 0) {
				run_tar1v_count(ii, 0) += 1;
				run_tarxyv1(ii, arma::span(1, 2)) += (tar1_avgxyv(ii, arma::span(1, 2)) - truthkxyv(ii, arma::span(1, 2))) % (tar1_avgxyv(ii, arma::span(1, 2)) - truthkxyv(ii, arma::span(1, 2)));
			}
		}
		for (int ii = 0; ii < tar2_avgxyv.n_rows; ii++) {
			if (tar2_avgxyv(ii, 1) != 0) {
				run_tar2v_count(ii, 0) += 1;
				run_tarxyv2(ii, arma::span(1, 2)) += (tar2_avgxyv(ii, arma::span(1, 2)) - truthkxyv2(ii, arma::span(1, 2))) % (tar2_avgxyv(ii, arma::span(1, 2)) - truthkxyv2(ii, arma::span(1, 2)));
			}
		}
		
		avg_late1 += late1;
		avg_late2 += late2;
		avg_ft_rate += ft_rate;
		/*for (int ii = 0; ii < tar1.n_rows; ii++) {
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
		std::cout << "tar2 " << tar2 << std::endl;*/

	}
	for (int ii = 0; ii < run_tar1.n_rows; ii++) {
		if (run_tar1(ii, 1) != 0) {		
			run_tar1(ii, arma::span(1, 2)) = run_tar1(ii, arma::span(1, 2)) / run_tar1_count(ii, 0);
			run_tarxy1(ii, arma::span(1, 2)) = run_tarxy1(ii, arma::span(1, 2)) / run_tar1_count(ii, 0);
			run_tar1_avg(ii, arma::span(1, 2)) = run_tar1_avg(ii, arma::span(1, 2)) / run_tar1_count(ii, 0);
		}
	}
	for (int ii = 0; ii < run_tar2.n_rows; ii++) {
		if (run_tar2(ii, 1) != 0) {
			run_tar2(ii, arma::span(1, 2)) = run_tar2(ii, arma::span(1, 2)) / run_tar2_count(ii, 0);
			run_tarxy2(ii, arma::span(1, 2)) = run_tarxy2(ii, arma::span(1, 2)) / run_tar2_count(ii, 0);
			run_tar2_avg(ii, arma::span(1, 2)) = run_tar2_avg(ii, arma::span(1, 2)) / run_tar2_count(ii, 0);
			//run_tar2_avg.row(ii)= run_tar2_avg.row(ii) / run_tar2_count(ii, 0);
		}
	}
	for (int ii = 0; ii < run_tarxyv1.n_rows; ii++) {
		if (run_tarxyv1(ii, 1) != 0) {
			run_tarxyv1(ii, arma::span(1, 2)) = run_tarxyv1(ii, arma::span(1, 2)) / run_tarxyv1(ii, 0);
		}
	}
	for (int ii = 0; ii < run_tarxyv2.n_rows; ii++) {
		if (run_tarxyv2(ii, 1) != 0) {
			run_tarxyv2(ii, arma::span(1, 2)) = run_tarxyv2(ii, arma::span(1, 2)) / run_tarxyv2(ii, 0);
		}
	}
	run_tar1(arma::span(0, run_tar1.n_rows - 1), arma::span(1, 2)) = sqrt(run_tar1(arma::span(0, run_tar1.n_rows - 1), arma::span(1, 2)));
	run_tar2(arma::span(0, run_tar2.n_rows - 1), arma::span(1, 2)) = sqrt(run_tar2(arma::span(0, run_tar2.n_rows - 1), arma::span(1, 2)));
	run_tarxy1(arma::span(0, run_tarxy1.n_rows - 1), arma::span(1, 2)) = sqrt(run_tarxy1(arma::span(0, run_tarxy1.n_rows - 1), arma::span(1, 2)));
	run_tarxy2(arma::span(0, run_tarxy2.n_rows - 1), arma::span(1, 2)) = sqrt(run_tarxy2(arma::span(0, run_tarxy2.n_rows - 1), arma::span(1, 2)));
	run_tarxyv1(arma::span(0, run_tarxyv1.n_rows - 1), arma::span(1, 2)) = sqrt(run_tarxyv1(arma::span(0, run_tarxyv1.n_rows - 1), arma::span(1, 2)));
	run_tarxyv2(arma::span(0, run_tarxyv2.n_rows - 1), arma::span(1, 2)) = sqrt(run_tarxyv2(arma::span(0, run_tarxyv2.n_rows - 1), arma::span(1, 2)));
	
	std::cout <<"after "<< run_tar1 << std::endl;
	std::cout << "after xy" << run_tarxy1 << std::endl;
	std::cout << "after xy v" << run_tarxyv1 << std::endl;
	//run_tar1_avg = run_tar1_avg / nruns;
	//run_tar2_avg = run_tar2_avg / nruns;
	run_tar1_avg.col(0) = run_tar1_avg.col(0) / 2;
	run_tar2_avg.col(0) = run_tar2_avg.col(0) / 2;
	avg_late1 = avg_late1 / nruns;
	avg_late2 = avg_late2 / nruns;
	avg_ft_rate = avg_ft_rate / nruns;
	std::cout << "avg late " << avg_late1 << "  " << avg_late2 << std::endl;
	std::cout << "avg ft rate " << avg_ft_rate << std::endl;
	
	std::cout << run_tar1_avg << std::endl;
	std::cout << run_tar2_avg << std::endl;
	for (int ii = 0; ii < run_tar1_avg.n_rows - 1; ii++) {
		if (ii >= avg_late1) {
			if (run_tar1_avg(ii, 1) != 0 && run_tar1_avg(ii + 1, 1) != 0) {
				cv::Point ppoint(run_tar1_avg(ii, 1), run_tar1_avg(ii, 2));
				cv::Point cpoint(run_tar1_avg(ii + 1, 1), run_tar1_avg(ii + 1, 2));
				cv::line(graph, ppoint, cpoint, cv::Scalar(0, 255, 255), 2);
			}
		}
		
	}
	for (int ii = 0; ii < run_tar2_avg.n_rows - 1; ii++) {
		if (ii >= avg_late2) {
			if (run_tar2_avg(ii, 1) != 0 && run_tar2_avg(ii + 1, 1) != 0) {
				cv::Point ppoint(run_tar2_avg(ii, 1), run_tar2_avg(ii, 2));
				cv::Point cpoint(run_tar2_avg(ii + 1, 1), run_tar2_avg(ii + 1, 2));
				cv::line(graph, ppoint, cpoint, cv::Scalar(0, 255, 0), 2);
			}
		}
		
	}
	cv::imshow("Output", frame);
	cv::imshow("Output2", graph);
	cv::imshow("Output3", frame);
	cv::waitKey(0);

	


	//meas_data.save("meas.txt", arma::arma_ascii);
	std::cout << "Done" << std::endl;
	while (!GetAsyncKeyState(VK_ESCAPE)) {}
}