#pragma once
#include <armadillo>
#include <opencv2\opencv.hpp>

class runs
{
public:
	
	~runs();
	void graphavg(cv::Mat &, arma::mat, int, double);
	//void rmse(arma::mat, arma::mat, double);
};

