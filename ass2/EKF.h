#pragma once
#include <armadillo>
#include <opencv2\opencv.hpp>
#include <chrono>
#include <random>
#include <math.h>
#include <cmath>

#include "model.h"
#include "IPDA.h"


class EKF{
	model kfmodel;
	IPDA ipda;

	int maxdetect = 50;
	double com_th = 0.9;
	double term_th = 0.05;
	double init_quality = 0.2;

	int start_all_frame_num;
	double sigv = 2;

	arma::mat hz_store;
	arma::mat id_store;
	arma::mat sk_store;

	int com_size;
	arma::mat com_id;
	arma::mat com_color;

	arma::mat asso_data;

public:
	void kf(cv::Mat &,int, double, int, arma::mat, arma::mat &, 
				arma::mat &, arma::mat &, arma::mat &, arma::mat &,
				arma::mat &, double, double,
				arma::mat &, arma::mat &, bool);
	//void ipda();
};

