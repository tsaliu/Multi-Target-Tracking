#pragma once
#include <opencv2\opencv.hpp>
#include <opencv2\highgui.hpp>
#include <opencv2\imgcodecs.hpp>
#include <opencv2\imgproc.hpp>
#include <opencv2\videoio.hpp>
#include <opencv2\video.hpp>
#include <armadillo>


class GenTar {
	double pos_x, pos_y;
	double noise_x, noise_y;
	double vel_x, vel_y;
	int length;
	arma::mat sim;

	

public:
	void init_pos(double, double);
	void noise(double, double);
	void init_vel(double, double);
	void tot_length(int);
	void gen(void);
	void data(arma::mat &);
	void print_sim(void);
	void graph(cv::Mat &, arma::mat, int);
};

