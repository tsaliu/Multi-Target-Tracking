#pragma once
#include <armadillo>
#include <math.h>
#include <opencv2\opencv.hpp>
#include <random>

class Sen{
	arma::mat Tdata;
	arma::mat Txydata;
	arma::mat Tdata2;
	arma::mat Txydata2;
	int r_pos_x, r_pos_y;
	double error_r, error_a;
	double sample_time;
	double pd, fa_d;
	int num_fa;
	arma::mat fa;
	double pd_dis, pd_dis2;


	arma::mat sen_data;		//sen data with fa
	arma::mat Tsen_data;	//sen data without fa
	int tot_pts;

	

public:
	void getdata(arma::mat, double);
	void getdata2(arma::mat);
	void pos(int, int);
	void error(double, double);
	void samptime(double);
	void pdetect(double);
	void fa_density(double);
	void fa_gen(double , int);
	void gen(double, double);
	void graph(cv::Mat &, int, int, int);
	void meas_data(arma::mat &, int, int);
	void truth(arma::mat &, int);
	void truthxy(arma::mat &, int);
	void truth2(arma::mat &, int);
	void truthxy2(arma::mat &, int);
	void graphmeas(cv::Mat &, arma::mat, int, int, int);
};

