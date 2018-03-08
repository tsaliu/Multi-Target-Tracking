#include "EKF.h"
#include "xyraconv.h"
#include "randseed2.h"


EKF::EKF(){}
EKF::~EKF(){}

void EKF::init(int k, arma::mat meas_data, int radi) {
	if (k == 0) {
		arma::mat init_pos_ra(1, 2);
		arma::mat init_pos_xy(1, 2);
		init_pos_ra(0, 0) = meas_data(0, 0);
		init_pos_ra(0, 1) = meas_data(0, 1);
		ra2xy2(init_pos_ra, init_pos_xy, radi, sigv);
		kfmodel.init_pos(init_pos_xy(0, 0), init_pos_xy(0, 1));
		init_x = init_pos_xy(0, 0);
		init_y = init_pos_xy(0, 1);
		x = init_x;
		y = init_y;
		px = x;
		py = y;
		//std::cout << x << "   " << y << std::endl;
	}
}

void EKF::graphasso(cv::Mat &graphon, int k, double st) {
	if (k % (int)st == 0) {
		
		int n = k / st;
		cv::Point target_at(plotx, ploty);
		if (n == 0) {
			cv::line(graphon, target_at, target_at, cv::Scalar(255, 0, 0), 2);
			/*px = x;
			py = y;*/
		}
		if (n > 0) {

			cv::Point target_pat(plotpx, plotpy);
			cv::line(graphon, target_at, target_pat, cv::Scalar(255, 0, 0), 2);
			/*px = x;
			py = y;*/
			//std::cout << plotpx << "  " << plotpy << " " << x << " " << " " << y << std::endl;
		}
		n++;
		
		
	}
}

void EKF::start(int k, double st, int radi, arma::mat p00, arma::mat R, arma::mat &Pkk, arma::mat &hxk1k1, arma::mat &hzk1k, arma::mat &sk1) {

	std::default_random_engine gene;
	gene.seed(std::chrono::system_clock::now().time_since_epoch().count());
	double b = exp(-pow(sigv, 2) / 2);
	if ((k % (int)st == 0)) {
		if (k == 0) {
			skip = 0;
		}
		kfmodel.getdata(skip, radi, k, st, sigv, F, Q, x, y, dx, dy, da, dr);

		arma::mat H = arma::zeros(2, 4);
		H(0, 0) = 1;
		H(1, 2) = 1;
		arma::mat Hx = arma::zeros(2, 4);

			arma::mat w = arma::zeros(4, 4);
			for (int i = 0; i<w.n_rows; i++) {
				for (int j = 0; j<w.n_cols; j++) {
					arma::mat rn = arma::randn(4, 4);
					std::normal_distribution<double> distn(0, sqrt(Q(i, j)));
					//std::cout << rn << std::endl;
					//w(i, j) = rn(i, j)*sqrt(Q(i, j));
					w(i, j) = distn(gene);
				}
			}
			std::normal_distribution<double> distn2(0, sqrt(R(0, 0)));
			std::normal_distribution<double> distn3(0, sqrt(R(1, 1)));
			arma::mat rnn = arma::randn(2, 2);
		
			arma::mat v = arma::zeros(2, 4);
			v(0, 0) = distn2(gene);
			v(1, 2) = distn3(gene);
			
			////////////////////
			//std::cout << "dx  " <<da<< std::endl;
			//std::cout << "dy  " <<dr<< std::endl;
			///////////////////
			arma::mat xx = arma::zeros(4, 4);
			arma::mat x00 = arma::zeros(4, 4);

			arma::mat zk1 = arma::zeros(2, 4);
			arma::mat ra_tmp(1, 2);
			arma::mat xy_tmp(1, 2);
			arma::mat xx_tmp = arma::zeros(4, 4);

			arma::mat hxkk = arma::zeros(4, 4);
			arma::mat pk1k1 = arma::zeros(4, 4);
			if (k == 0) {
				x00(0, 0) = init_x;
				x00(2, 2) = init_y;
				x00(1, 1) = 0;
				x00(3, 3) = 0;

				xx = x00 + w;

				xy_tmp(0, 0) = xx(0, 0);
				xy_tmp(0, 1) = xx(2, 2);
				xy2ra2(xy_tmp, ra_tmp, radi, sigv);
				xx_tmp(0, 0) = ra_tmp(0, 0);
				xx_tmp(2, 2) = ra_tmp(0, 1);

				xy_tmp(0, 0) = xx(1, 1);
				xy_tmp(0, 1) = xx(3, 3);
				xy2ra2(xy_tmp, ra_tmp, radi, sigv);
				xx_tmp(1, 1) = da;
				xx_tmp(3, 3) = dr;
				Hx = H * xx_tmp;
				zk1 = H * xx_tmp + v;

				hxkk = xx;
				///////////////////////
				//std::cout << "00000000" << xx_tmp << std::endl;
				//std::cout << H * xx_tmp << std::endl;
				////////////////////////
				Pkk = p00;
			}
			else {
				xx(0, 0) = x;
				xx(2, 2) = y;
				xx(1, 1) = dx;
				xx(3, 3) = dy;

				xx = xx + w;

				xy_tmp(0, 0) = xx(0, 0);
				xy_tmp(0, 1) = xx(2, 2);
				xy2ra2(xy_tmp, ra_tmp, radi, sigv);
				xx_tmp(0, 0) = ra_tmp(0, 0);
				xx_tmp(2, 2) = ra_tmp(0, 1);

				xy_tmp(0, 0) = xx(1, 1);
				xy_tmp(0, 1) = xx(3, 3);
				xy2ra2(xy_tmp, ra_tmp, radi, sigv);
				xx_tmp(1, 1) = 0;
				xx_tmp(3, 3) = dr;
				Hx = H * xx_tmp;
				zk1 = H * xx_tmp + v;

				hxkk = xx;
				Pkk = p00;
			}
			////////////////////////
			//std::cout<<"hxkk     "<<hxkk<<std::endl;
			///////////////////
			plot_est_x = hxkk(0, 0);
			plot_est_y = hxkk(2, 2);

			arma::mat hxk1k = arma::zeros(4, 4);
			hxk1k = F * hxkk;
			//arma::mat hzk1k = arma::zeros(2, 4);

			arma::mat hxk1k_tmp = arma::zeros(4, 4);
			xy_tmp(0, 0) = hxk1k(0, 0);
			xy_tmp(0, 1) = hxk1k(2, 2);
			
			xy2ra3(xy_tmp, ra_tmp, radi, sigv);
			hxk1k_tmp(0, 0) = ra_tmp(0, 0);
			hxk1k_tmp(2, 2) = ra_tmp(0, 1);
			/////////////////////////////
			//std::cout << hxk1k_tmp << std::endl;
			////////////////////////////
			/*hca = ra_tmp(0, 0);
			hcr = ra_tmp(2, 2);
			if (k == 0) {
				hpa = 0;
				hpr = 0;
			}
*/
			xy_tmp(0, 0) = hxk1k(1, 1);
			xy_tmp(0, 1) = hxk1k(3, 3);
			xy2ra2(xy_tmp, ra_tmp, radi, sigv);
			////////////
			hxk1k_tmp(1, 1) = da;
			hxk1k_tmp(3, 3) = dr;

			Hx = H * hxk1k_tmp;
			//////////////////////////
			//std::cout << "HHH" << Hx<< std::endl;
			/////////////////////////
			hzk1k = H * hxk1k_tmp;
			arma::mat vk1 = arma::zeros(2, 4);
			vk1 = zk1 - hzk1k;
			
			arma::mat pk1k = arma::zeros(4, 4);
			pk1k = F * Pkk*F.t() + Q;

			//float sk1;
			//arma::mat sk1 = arma::zeros(2, 2);
			sk1 = R + Hx * pk1k*Hx.t();
			////sk1 = sk11(0, 0);
			arma::mat wk1 = arma::zeros(4, 2);
			///////////////////////////
			//std::cout <<"S   "<< sk1 << std::endl;
			/////////////////////////
			wk1 = pk1k * Hx.t() * arma::inv(sk1);

			
			pk1k1 = pk1k - wk1 * sk1*wk1.t();
			//P = pk1k1;


			hxk1k1 = hxk1k + wk1 * vk1;
			/*hxk1k1(0,0)=hxk1k1(0,0) *b;
			hxk1k1(2, 2) = hxk1k1(2, 2) * pow(b,1);*/
			//D = vk1.t()*arma::inv(sk1)*vk1;
			//std::cout << pk1k1 << std::endl;
			/////////////////////////////
			//std::cout << hxk1k1 << std::endl;
			//std::cout << "END" << std::endl;
			/////////////////////////////
			
		//}		
	}
	
}

void EKF::asso(int k, double st, arma::mat meas_data, int numfa, arma::mat hzk1k, arma::mat sk1, int radi) {
	if (k == 0) {
		skip = 0;
		
	}
	if ((k % (int)st == 0) && (k != 0)) {
		arma::mat z_tmp(1, 2);
		arma::mat x_tmp(1, 2);
		z_tmp(0, 0) = hzk1k(0, 0);
		z_tmp(0, 1) = hzk1k(1, 2);
		ra2xy2(z_tmp, x_tmp, radi, sigv);
		plotpx = x_tmp(0, 0);
		plotpy = x_tmp(0, 1);


			int max_size = 3 * numfa;
			n = k / st;
			
			/////////////
			//std::cout << hzk1k << std::endl;
			//////////////////
			//extracting usefull data
			//std::cout << meas_data << std::endl;
			arma::uvec meas_a_at = arma::find(meas_data(arma::span(0, max_size - 1), (k / st) * 2 + 0) != 0);
			//std::cout << meas_a_at << std::endl;
			int meas_size = arma::max(meas_a_at);
			arma::vec meas_a;
			arma::vec meas_r;
			
			meas_a = meas_data(arma::span(0, meas_size), (k / st) * 2 + 0);
			meas_r = meas_data(arma::span(0, meas_size), (k / st) * 2 + 1);
			arma::mat zk1 = arma::zeros(2, 4);
			
			int min_d_at = 0;
			double min_d = 0;
			for (int i = 0; i < meas_size; i++) {
				zk1(0, 0) = meas_a(i, 0);
				zk1(1, 2) = meas_r(i, 0);

				arma::mat v = zk1 - hzk1k;
				arma::mat D = v.t() * arma::inv(sk1) * v;
				//////////////////////
				//std::cout <<"D   " << arma::sum(arma::sum(D)) << std::endl;
				////////////////////
				
				double dis = arma::sum(arma::sum(D));
				
				if (i == 0) {
					min_d = dis;
					min_d_at = i;
				}
				else {
					if (dis < min_d) {
						min_d = dis;
						min_d_at = i;
					}
				}
			}

			if (min_d <= gamma) {
				/////////////////
				//std::cout << "DIS    " << min_d_at << "  " << min_d << std::endl;
				//////////////
				arma::mat a_zk1(1, 2);
				a_zk1 << meas_a(min_d_at, 0) << meas_r(min_d_at, 0) << arma::endr;
				/////////////////////
				//std::cout << "asso    " << a_zk1 << std::endl;
				//////////////////
				arma::mat xy_tmp(1, 2);
				ra2xy2(a_zk1, xy_tmp, radi, sigv);

				x = xy_tmp(0, 0);
				y = xy_tmp(0, 1);

				plotx = x;
				ploty = y;
				////////////////
				//std::cout << "xy  " << x << "   " << y << std::endl;
				/////////////////
				skip = 0;
			}
			else {
				arma::mat a_zk1(1, 2);
				a_zk1 << hzk1k(0, 0) << hzk1k(1, 2) << arma::endr;
				arma::mat xy_tmp(1, 2);
				ra2xy2(a_zk1, xy_tmp, radi, sigv);

				x = xy_tmp(0, 0);
				y = xy_tmp(0, 1);
				plotx = x;
				ploty = y;
				//////////////////
				//std::cout << hzk1k << std::endl;
				//std::cout << "skip xy  " << x << "   " << y << std::endl;
				///////////////////
				skip++;
			}
			
			//std::cout << sk1 << std::endl;
			
			
			
			//std::cout <<"meas_ar   "<< meas_a<<"    "<< meas_r << std::endl;
			//std::cout << hzk1k << std::endl;
	}
}

void EKF::graph_est(cv::Mat &graphon, int k, double st) {
	if (k % (int)st == 0) {
		
		int n = k / st;
			
		cv::Point target_at(plot_est_x, plot_est_y);
		if (n == 0) {
			cv::line(graphon, target_at, target_at, cv::Scalar(0, 255, 0), 2);
			plot_est_px = plot_est_x;
			plot_est_py = plot_est_y;
		}
		if (n > 0) {

			cv::Point target_pat(plot_est_px, plot_est_py);
			cv::line(graphon, target_at, target_pat, cv::Scalar(0, 255, 0), 2);
			plot_est_px = plot_est_x;
			plot_est_py = plot_est_y;

			//std::cout << plotpx << "  " << plotpy << " " << x << " " << " " << y << std::endl;
		}
		n++;
	}
}

void EKF::asso_pda() {

}