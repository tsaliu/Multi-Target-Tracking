#include "runs.h"



runs::~runs(){}

void runs::graphavg(cv::Mat &graphon, arma::mat avgxy, int tot_len, double st) {
	int px = 0;
	int py = 0;
	int x = 0;
	int y = 0;
	int len = tot_len / (int)st;
	//if (!avgxy.is_empty()) {
		for (int i = 0; i < len - 1; i++) {
			px = avgxy(i, 0);
			py = avgxy(i, 1);

			if ((i + 1) <= len - 1) {
				x = avgxy(i + 1, 0);
				y = avgxy(i + 1, 1);
			}
			else {
				x = px;
				y = py;
			}


			cv::Point target_at(x, y);
			cv::Point target_pat(px, py);
			cv::line(graphon, target_at, target_pat, cv::Scalar(0, 255, 255), 2);


		}
	//}

}

//void rmse(arma::mat avgxy, arma::mat truth, double st) {
//
//}