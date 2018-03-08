#pragma once
#include <armadillo>
#include <opencv2\opencv.hpp>
#include "model.h"
#include <chrono>


class EKF
{
	model kfmodel;
	double init_x, init_y;

	double x, y, dx, dy, da, dr, px, py, plotpx, plotpy, plotx, ploty;
	double hca, hcr, hpa, hpr, hda, hdr;
	double gate_x = 50, gate_y = 50;
	double sigv = 0.2;
	//double R = 0.0001;
	int n;
	arma::mat ra;
	arma::mat xy;

	arma::mat F;
	arma::mat Q;

	//arma::mat Pkk;
	//arma::mat pk1k1;
	//arma::mat hxk1k1;
	arma::mat hX;

	//arma::mat hzk1k;

	

	double plot_est_x, plot_est_y, plot_est_px, plot_est_py;

	double min_dis;
	int min_dis_num;

	double pxx, pyy;

	double gamma = 0.05;
	int skip;
	
	

public:
	EKF();
	~EKF();

	void init(int, arma::mat, int);
	//void asso(int, double, arma::mat, int, int);
	void asso(int, double, arma::mat, int, arma::mat, arma::mat, int);
	void graphasso(cv::Mat &, int, double);
	void start(int, double, int, arma::mat, arma::mat, arma::mat &, arma::mat &, arma::mat &, arma::mat &);
	void graph_est(cv::Mat &, int, double);
	void asso_pda();
};

