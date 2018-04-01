#include "EKF.h"
#include "xyraconv.h"
#include <Windows.h>
#include <iostream>
#include <stdio.h>


void EKF::kf(cv::Mat &frame, int k, double st, int radi, arma::mat meas_data, arma::mat &P,
	arma::mat &phxk1k1, arma::mat &chxk1k1, arma::mat &hzk1k, arma::mat &sk1,
	arma::mat &t_id, double lambda, double pd,
	arma::mat &q, arma::mat &out_save, bool ipda_mode) {
	
	ipda.getpara(lambda, pd);
	//arma::mat asso_data;

	arma::mat p00 = arma::zeros(4, 4);
	int init_pos_var = 200;
	int init_vel_var = 100;
	p00 << init_pos_var << 0 << 0 << 0 << arma::endr
		<< 0 << init_vel_var << 0 << 0 << arma::endr
		<< 0 << 0 << init_pos_var << 0 << arma::endr
		<< 0 << 0 << 0 << init_vel_var << arma::endr;
	arma::mat R = arma::zeros(2, 2);
	double init_a_var = 0;
	double init_r_var = 0;
	if (ipda_mode) {
		init_a_var = 0.05;
		init_r_var = 0.5;
	}
	else {
		init_a_var = 0.1;
		init_r_var = 5;
	}
	R << init_a_var << 0 << arma::endr
		<< 0 << init_r_var << arma::endr;

	int n = k / st * 2;

	std::default_random_engine gene;
	gene.seed(std::chrono::system_clock::now().time_since_epoch().count());

	double b = exp(-pow(sigv, 2) / 2);

	arma::mat::col_iterator cita = meas_data.begin_col(n); 
	arma::mat::col_iterator cita_end = meas_data.end_col(n);  

	arma::mat::col_iterator citr = meas_data.begin_col(n + 1);
	arma::mat::col_iterator citr_end = meas_data.end_col(n + 1);

	arma::mat cmeas_tmp(1, 2);
	arma::mat cstate_tmp(1, 2);
	arma::mat sk1_tmp(2, 2);

	arma::uvec id_check_m = arma::find(t_id(arma::span(0, maxdetect - 1), 0) != 0); //means just started, nothing ran yet
	bool id_check_start = id_check_m.is_empty();
	if (id_check_start) {
		start_all_frame_num = k / st;
		int i = 1;
		for (; cita != cita_end; cita++) {		//increment for each row
			if ((*cita != 0)) {				//eliminate the reserved memory (zeros)
				t_id(i - 1, k / st) = i;			//track id record

				//com_id = arma::zeros(1, 1);			//initiate comfirm track parameters
				com_size = 0;


				//std::cout << (*cita) << " what " << (*citr) << std::endl;
				cmeas_tmp << (*cita) << (*citr) << arma::endr;
				ra2xy2(cmeas_tmp, cstate_tmp, radi, sigv);
				chxk1k1((i - 1) * 4 + 0, 0) = cstate_tmp(0, 0);
				chxk1k1((i - 1) * 4 + 2, 2) = cstate_tmp(0, 1);
				
				P(arma::span((i - 1) * 4, (i - 1) * 4 + 3), arma::span(0, 3)) = p00;
				
				hzk1k((i - 1) * 2 + 0, 0) = (*cita);
				hzk1k((i - 1) * 2 + 1, 2) = (*citr);

				

				sk1_tmp << init_a_var   << 0 << arma::endr
					<< 0 << init_r_var  << arma::endr;

				sk1(arma::span((i - 1) * 2, (i - 1) * 2 + 1), arma::span(0, 1)) = sk1_tmp;
			

				i++;
				citr++;
			}
		}
		//asso_data.resize(i - 1, 6);
		//asso_data.col(0) = t_id(arma::span(0, i - 2), k / st);


		phxk1k1 = chxk1k1;
		hz_store = hzk1k;
		id_store = t_id;
		sk_store = sk1;
		
		
	}
	else {
		phxk1k1.fill(0);
		chxk1k1.fill(0);
		arma::mat save_com = arma::zeros(0, 7);
		//if (k / st - 1 == start_all_frame_num) {
			//std::cout << meas_data << std::endl;
		if (ipda_mode) {
			ipda.ipda(meas_data, cita, hz_store, id_store, sk_store, k, st, q, asso_data, radi, sigv, P);
		}
		else {
			ipda.nn(meas_data, cita, hz_store, id_store, sk_store, k, st, q, asso_data, radi, sigv, P);
		}
			
			
			/*if (k == initt2) {
				asso_data.resize(asso_data.n_rows + 1, 6);
				arma::mat::col_iterator cita = meas_data.begin_col(n);
				arma::mat::col_iterator citr = meas_data.begin_col(n + 1);
				cita++;
				citr++;

				asso_data(asso_data.n_rows + 1, 0) = t_id(t_id.index_max()) + 1;
				asso_data(asso_data.n_rows + 1, 2) = (*cita);
				asso_data(asso_data.n_rows + 1, 3) = (*citr);
				asso_data(asso_data.n_rows + 1, 4) = (*cita);
				asso_data(asso_data.n_rows + 1, 5) = (*citr);
				
			}*/
			/* attempt to create new asso data when empty
			if (asso_data.is_empty()) {
				arma::mat::col_iterator cita = meas_data.begin_col(n);
				arma::mat::col_iterator cita_end = meas_data.end_col(n);

				arma::mat::col_iterator citr = meas_data.begin_col(n + 1);
				arma::mat::col_iterator citr_end = meas_data.end_col(n + 1);
				arma::mat p_re = arma::zeros(4 * maxdetect, 4);
				int i = 0;
				for (; cita != cita_end; cita++) {
					if ((*cita != 0)) {
						arma::mat re_create(1, 6);
						re_create << t_id(t_id.index_max()) + i + 1 << init_quality
							<< (*cita) << (*citr)
							<< (*cita) << (*citr) << arma::endr;
						asso_data.insert_rows(asso_data.n_rows, re_create);


						p_re(arma::span(i * 4, i * 4 + 3), arma::span(0, 3)) = p00;
						citr++;
						i++;
					}
				}
				p_re.resize(4 * maxdetect, 4);
				P = p_re;
			}*/



			//ipda.ipda(meas_data, cita, hz_store, id_store, sk_store, 0, st, q, asso_data);
			//std::cout <<" asso_data " << asso_data << std::endl;
			t_id(arma::span(0, asso_data.n_rows - 1), k / st) = asso_data.col(0);
			//std::cout << asso_data.col(0) << std::endl;
			//std::cout << t_id.col(k / st) << std::endl;
			//P.resize
			arma::mat P_temp = P;
			P_temp.resize(4 * asso_data.n_rows, 4);
			int past_size = 0;
			int add_new_size = 0;
			for (int i = 0; i < asso_data.n_rows; i++) {
				arma::uvec find_past_data = arma::find(t_id.col((k - 1) / st) == asso_data(i, 0));
				
				if (!find_past_data.is_empty()) {
					int find_past_data_at = find_past_data(0);
					past_size++;
					//std::cout << find_past_data_at << std::endl;
					P_temp(arma::span(i * 4, i * 4 + 3), arma::span(0, 3)) = P(arma::span(find_past_data_at * 4, find_past_data_at * 4 + 3), arma::span(0, 3));
				}
				else {
					add_new_size++;
					P_temp(arma::span(i * 4, i * 4 + 3), arma::span(0, 3)) = p00;
				}
				
				//extract points since ipda output z, need to convert to x
				arma::mat ar_tmp(1, 2);
				arma::mat xy_tmp(1, 2);
				ar_tmp << asso_data(i, 2) << asso_data(i, 3) << arma::endr;
				ra2xy2(ar_tmp, xy_tmp, radi, sigv);
				double px = xy_tmp(0);
				double py = xy_tmp(1);

				ar_tmp << asso_data(i, 4) << asso_data(i, 5) << arma::endr;
				ra2xy2(ar_tmp, xy_tmp, radi, sigv);
				double cx = xy_tmp(0);
				double cy = xy_tmp(1);
				phxk1k1(i * 4, 0) = px;
				phxk1k1(i * 4 + 2, 2) = py;
				chxk1k1(i * 4, 0) = cx;
				chxk1k1(i * 4 + 2, 2) = cy;

				double dx, dy, da, dr;
				double ca = asso_data(i, 4);
				double cr = asso_data(i, 5);
				double pa = asso_data(i, 2);
				double pr = asso_data(i, 3);
				arma::mat F(4, 4);
				arma::mat Q(4, 4);
				kfmodel.getdata(radi, k, st, sigv, F, Q,
					cx, cy, px, py, dx, dy, da, dr);
				//std::cout << da << "  " << asso_data(i, 4) - asso_data(i, 2) << std::endl;
				
				//std::cout << Q << std::endl;
			



				//start of kalman filter
				std::default_random_engine gene;
				gene.seed(std::chrono::system_clock::now().time_since_epoch().count());
				double b = exp(-pow(sigv, 2) / 2);
				arma::mat H = arma::zeros(2, 4);
				H(0, 0) = 1;
				H(1, 2) = 1;
				arma::mat Hx = arma::zeros(2, 4);

				arma::mat w = arma::zeros(4, 4);
				for (int i = 0; i<w.n_rows; i++) {
					for (int j = 0; j<w.n_cols; j++) {
						arma::mat rn = arma::randn(4, 4);
						if (Q(i, j) != 0) {
							std::normal_distribution<double> distw(0, sqrt(Q(i, j)));
							//std::cout << rn << std::endl;
							//w(i, j) = rn(i, j)*sqrt(Q(i, j));
							w(i, j) = distw(gene);
						}
						else {
							w(i, j) = 0;
						}
						
					}
				}
				arma::mat v = arma::zeros(2, 4);
				for (int i = 0; i < R.n_rows; i++) {

					std::normal_distribution<double> distv(0, sqrt(R(i, i)));
					v(i, i * 2) = distv(gene);
				}

				arma::mat xx = arma::zeros(4, 4);
				arma::mat x00 = arma::zeros(4, 4);

				arma::mat zk1 = arma::zeros(2, 4);
				arma::mat zk1_tmp = arma::zeros(4, 4);
				//arma::mat ra_tmp(1, 2);
				//arma::mat xy_tmp(1, 2);

				arma::mat hxkk = arma::zeros(4, 4);
				arma::mat pk1k1 = arma::zeros(4, 4);

				xx(0, 0) = cx;
				xx(2, 2) = cy;
				xx(1, 1) = dx;
				xx(3, 3) = dy;

				x00 = xx;
				xx = xx + w;

				zk1_tmp(0, 0) = ca;
				zk1_tmp(2, 2) = cr;
				zk1_tmp(1, 1) = da;
				zk1_tmp(3, 3) = dr;

				zk1 = H * zk1_tmp + v;

				hxkk = xx;
				//P = pk1k1; //P=p00;
				arma::mat Pkk = arma::zeros(4, 4);
				Pkk = P_temp(arma::span(i * 4, i * 4 + 3), arma::span(0, 3));

				arma::mat hxk1k = arma::zeros(4, 4);
				hxk1k = F * hxkk;
				
				//arma::mat xy_tmp(1, 2);
				//arma::mat ar_tmp(1, 2);
				arma::mat hxk1k_tmp = arma::zeros(4, 4);
				xy_tmp(0, 0) = hxk1k(0, 0);
				xy_tmp(0, 1) = hxk1k(2, 2);
				//xy_tmp << hxk1k(0, 0) << hxk1k(2, 2) << arma::endr;
				xy2ra2(xy_tmp, ar_tmp, radi, sigv);

				hxk1k_tmp(0, 0) = ar_tmp(0, 0);
				hxk1k_tmp(2, 2) = ar_tmp(0, 1);
				//hxk1k_tmp

				//calc da dr
				double x_tmp = (cx + dx) - radi;
				double y_tmp = -((cy + dy) - radi);
				double r_tmp = sqrt(pow(x_tmp, 2) + pow(y_tmp, 2));
				double a_tmp = atan2(y_tmp, x_tmp);
				double da_tmp = (a_tmp - ca) / st;
				double dr_tmp = (r_tmp - cr) / st;
				//std::cout << da_tmp << "da dr tmpppppp   " << dr_tmp << std::endl;
				hxk1k_tmp(1, 1) = da_tmp;
				hxk1k_tmp(3, 3) = dr_tmp;

				
				Hx << -(cy / (pow(cx, 2) + pow(cy, 2))) << 0 << cx / (pow(cx, 2) + pow(cy, 2)) << 0 << arma::endr
					<< cx / sqrt(pow(cx, 2) + pow(cy, 2)) << 0 << cy / sqrt(pow(cx, 2) + pow(cy, 2)) << 0 << arma::endr;
				arma::mat hzk1k_tmp = arma::zeros(2, 4);
				hzk1k_tmp = H * hxk1k_tmp;

				
				arma::mat vk1 = arma::zeros(2, 4);
				vk1 = zk1 - hzk1k_tmp;

				arma::mat pk1k = arma::zeros(4, 4);
				pk1k = F * Pkk*F.t() + Q;

				//arma::mat sk1_tmp = arma::zeros(2, 2);
				//arma::mat pk1k_xy_tmp(1, 2);
				//arma::mat pk1k_ar_tmp(1, 2);
				//pk1k_xy_tmp << pk1k(0, 0) << pk1k(2, 2) << arma::endr;
				//xy2ra2(pk1k_xy_tmp, pk1k_ar_tmp, radi, sigv);
				//arma::mat pk1k_tmp = arma::zeros(4, 4);
				//pk1k_tmp = pk1k;
				//pk1k_tmp(0, 0) = pk1k_ar_tmp(0, 0);
				//pk1k_tmp(2, 2) = pk1k_ar_tmp(0, 1);
				//sk1_tmp = R + H * pk1k_tmp*H.t();
				////std::cout << Hx << std::endl;
				//sk1_tmp = abs(sk1_tmp);
				////sk1_tmp << init_a_var + init_pos_var << init_vel_var << arma::endr
				//	//<< init_vel_var << init_r_var + init_pos_var << arma::endr;

				arma::mat sk12 = arma::zeros(2, 2);
				sk12 = R + Hx * pk1k*Hx.t();

				/*
				arma::mat convxy(1, 2);
				arma::mat convar(1, 2);
				convar << sk1_tmp(0, 0) << sk1_tmp(1, 1) << arma::endr;
				ra2xy2(convar, convxy, radi, sigv);
				arma::mat sk1_tmp_tmp = arma::zeros(2, 2);
				sk1_tmp_tmp(0, 0) = convxy(0, 0);
				sk1_tmp_tmp(1, 1) = convxy(0, 1);
				sk1_tmp_tmp = abs(sk1_tmp_tmp);
				std::cout << "sk1 tmp" << sk1_tmp << std::endl;
				std::cout << "sk1 tmp tmp" << sk12 << std::endl;*/
				arma::mat wk1 = arma::zeros(4, 2);
				wk1 = pk1k * Hx.t()*arma::inv(sk12);

				//pk1k1 = pk1k - wk1 * sk1_tmp*wk1.t();
				//pk1k1 = p00;
				pk1k1 = (arma::eye(4, 4) - wk1 * Hx)*pk1k;
				//std::cout << "pk1k1   " << pk1k1 << std::endl;
				arma::mat hxk1k1 = arma::zeros(4, 4);
				hxk1k1 = hxk1k + wk1 * vk1;

				if (dx == 0) {
					//while (!GetAsyncKeyState(VK_SPACE)) {}
				}
				


				chxk1k1(arma::span(i * 4, i * 4 + 3), arma::span(0, 3)) = hxk1k1;

				hzk1k(arma::span(i * 2, i * 2 + 1), arma::span(0, 3)) = hzk1k_tmp;
				
				
				P_temp(arma::span(i * 4, i * 4 + 3), arma::span(0, 3)) = pk1k1;


				//hzk1k(i * 2, 0) = asso_data(i, 4);
				//hzk1k(i * 2 + 1, 2) = asso_data(i, 5);
					
				//sk1_tmp << init_a_var + init_pos_var << 0 << arma::endr
				//	<< 0 << init_r_var + init_pos_var << arma::endr;
				
				sk1(arma::span(i * 2, i * 2 + 1), arma::span(0, 1)) = sk12;

				
				std::default_random_engine genu;
				genu.seed(std::chrono::system_clock::now().time_since_epoch().count());
				std::uniform_int_distribution<int> distu(0, 255);

				//add comfirmed tracks
				arma::mat com_data;
				if (asso_data(i, 1) >= com_th) {
					arma::mat com_now_id(1, 1);
					com_now_id(0, 0) = asso_data(i, 0);
					arma::mat com_now_color(1, 3);
					com_now_color << distu(genu) << distu(genu) << distu(genu) << arma::endr;
					
					
					arma::uvec find_com = arma::find(com_id == com_now_id(0, 0));
					if (find_com.is_empty()) {
						com_id.insert_rows(com_id.n_rows, com_now_id);
						com_color.insert_rows(com_color.n_rows, com_now_color);
						com_data.insert_rows(com_data.n_rows, asso_data.row(i));
					}
					//for display com tracks
					//std::cout << "k = " << k << std::endl;
					//std::cout << "GGGGGGGGGGGGGGGGGGGGGGGGGGGG" << std::endl;
					//while (!GetAsyncKeyState(VK_SPACE)) {}
				


					
				}
				//if comfirmed
				arma::uvec find_if_com = arma::find(com_id == asso_data(i, 0));
				if (!find_if_com.is_empty()) {
					cv::Point com_track(hxk1k1(0, 0), hxk1k1(2, 2));
					cv::Point com_c(cx, cy);
					cv::Point com_p(px, py);
					//cv::line(frame, com_track, com_track, cv::Scalar(0, 255, 0), 2);
					//cv::line(frame, com_p, com_track, cv::Scalar(com_color(find_if_com(0), 0), com_color(find_if_com(0), 1), com_color(find_if_com(0), 2)), 2);
					//while (!GetAsyncKeyState(VK_SPACE)) {}
					cv::line(frame, com_p, com_c, cv::Scalar(com_color(find_if_com(0), 0), com_color(find_if_com(0), 1), com_color(find_if_com(0), 2)), 2);
					
					arma::mat save_tmp(1, 7);
					save_tmp(0, 0) = k;
					save_tmp(0, arma::span(1, 6)) = asso_data.row(i);
					save_com.insert_rows(save_com.n_rows, save_tmp);

					char idname[100];
					sprintf_s(idname, "ID: %i", (int)asso_data(i,0));
					putText(frame, idname, com_c, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 255), 1);
				
					//std::cout << "save " << save_com << std::endl;
					//while (!GetAsyncKeyState(VK_SPACE)) {}
				
				
				}
				//if comfirmed dropped
				if (asso_data(i, 1) <= term_th) {
					arma::uvec find_term_com = arma::find(com_id == asso_data(i, 0));
					if (!find_term_com.is_empty()) {
						
						//asso_data.shed_row(find_term_com(0));
						com_id.shed_row(find_term_com(0));
						com_color.shed_row(find_term_com(0));
						//while (!GetAsyncKeyState(VK_SPACE)) {}
					}
					
				}
				//std::cout <<"com id "<< com_id << std::endl;
				
			}


			P = P_temp;
			hzk1k.elem(arma::find_nonfinite(hzk1k)).zeros();
			id_store = t_id;
			hz_store = hzk1k;
			sk_store = sk1;
		//}
			out_save = save_com;
		
	}
	//cv::imshow("Output", frame);
	//cv::waitKey(1);
	//while (!GetAsyncKeyState(VK_SPACE)) {}

}

