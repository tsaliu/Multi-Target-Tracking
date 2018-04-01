#include "IPDA.h"
#include "chi2inv.h"
#include "calcdis.h"
#include "likelihood.h"
#include "xyraconv2.h"
#include <cmath>


void IPDA::getpara(double l, double p_d) {
	lambda = l;
	pd = p_d;
}


void IPDA::ipda(arma::mat meas_data,
		arma::mat::col_iterator ccita,
		arma::mat hzk1k, arma::mat t_id, arma::mat sk1, int k, double st,
		arma::mat &q, arma::mat &data_output,
		int radi, double sigv, arma::mat &P){
	int n = k / st * 2;
	arma::mat curr_meas = arma::zeros(maxdetect * 2, 4);
	

	gamma = chi2inv(pg, deg);
	arma::uvec find_id_e = arma::find(t_id.col((k - 1) / st) != 0); // 0 need to change to k/st
	//std::cout << t_id.col((k-1) / st) << std::endl;
	int find_id_size = arma::max(find_id_e) + 1;		//existing num of tracks
	//std::cout << "HERE" << std::endl;
	arma::mat likeihood = arma::zeros(find_id_size, maxdetect);
	arma::mat beta = arma::zeros(find_id_size, maxdetect);

	arma::mat dk_store_tmp = arma::zeros(maxdetect, 1);

	int current_meas_count = 0;
	if (*ccita != 0) {					//check if have measurement, by check the first row element of column
		//std::cout << "Have measurement " << k << std::endl;
		for (int i = 0; i < find_id_size; i++) { //loop for every track
			//reset iterator
			//std::cout << meas_data << std::endl;
			arma::mat::col_iterator cita = meas_data.begin_col(n);
			arma::mat::col_iterator cita_end = meas_data.end_col(n);
			arma::mat::col_iterator citr = meas_data.begin_col(n + 1);
			arma::mat::col_iterator citr_end = meas_data.end_col(n + 1);
			
			current_meas_count = 0;
			double tot_like = 0;
			
			for (; cita != cita_end; cita++) {	//increment for each row
				if ((*cita != 0)) {				//eliminate the reserved memory (zeros)

					arma::mat zk1(2, 2);
					zk1 << (*cita) << 0 << 0 << 0 << arma::endr
						<< 0 << 0 << (*citr) << 0 << arma::endr;
					curr_meas(current_meas_count * 2, 0) = (*cita);
					curr_meas(current_meas_count * 2 + 1, 2) = (*citr);
					
					double dis = calcdis(zk1, hzk1k, sk1, i, radi, sigv);	//calc distance to apply gate
					//std::cout << dis << std::endl;
					if (dis <= pow(gamma, 2)) {			//apply gate
						//measurement that pass the gate
						//std::cout << "pass gate" << std::endl;
						
						double like = likeicalc(sk1, zk1, hzk1k, i, radi, sigv) / pg;
						tot_like = tot_like + like;
						//std::cout << (*cita) << "  " << (*citr) << std::endl;
						
						//ignore too small likeihood
						if (like <= 3e-5) {
							like = 0;
						}
						likeihood(i, current_meas_count) = like;
						//std::cout << "like  loop" << current_meas_count << std::endl;
					}
					else {
						//std::cout << "no pass" << std::endl;
					}
					current_meas_count++;
					citr++;
				}
			}
			dk_mid = pd * pg * tot_like / lambda;	//passed all gated measurement
			dk_store_tmp(i, 0) = pd * pg - dk_mid;

			//std::cout << "tot like " << tot_like << std::endl;
		}					
		
	}
	else {
		//since no meas dk=pd pg
		for (int i = 0; i < find_id_size; i++) {
			std::cout << "No measurement " << k  << std::endl;
			dk_store_tmp(i, 0) = pg * pd;
		}
		likeihood.fill(0);
		
	}
	dk = dk_store_tmp;		//this is delta_k
	dk.resize(find_id_size, 1);
	likeihood.resize(find_id_size, current_meas_count);
	beta.resize(find_id_size, current_meas_count);
	curr_meas.resize(current_meas_count * 2, 4);
	//likeihood.elem(arma::find_nonfinite(likeihood)).zeros();
	//std::cout << "like  " << likeihood << std::endl;
	//std::cout <<"dk "<< dk << std::endl;
	
	//track quality, beta prob
	for (int i = 0; i < find_id_size; i++) {
		arma::mat::col_iterator meas_check = meas_data.begin_col(n);
		
		if (*meas_check == 0) {
			beta(i, 0) = (1 - pd * pg) / (1 - dk(i, 0)); //dk change
		}
		else {
			for (int ii = 0; ii < current_meas_count; ii++) {
				beta(i, ii) = (pd*pg / lambda * likeihood(i, ii)) / (1 - dk(i, 0));//dkchange
				//beta(i, ii) = abs(beta(i, ii));
			}
		}
		

		//std::cout << "asso? " << beta.row(i).index_max() << std::endl;
		//asso_m(i, beta.row(i).index_max()) += 1;

		
		bool hist_check = false;	//reset check
		arma::mat unitmax(1, 1);
		arma::uvec find_hist = arma::find(t_id == t_id(i, (k-1)/st)); //need ot change to k/st
		if (arma::size(find_hist) == arma::size(unitmax)) {
			hist_check = true;		//check if the track is just initialize by checking track history
		}
		if (hist_check == true) {
			q(i, k / st) = init_quality;		//need to change to k/st
		}
		else {
			q(i, k / st) = (1 - dk(i, 0))*(q(arma::max(find_hist))) / (1 - dk(i, 0) * q(arma::max(find_hist))); //need to change to k/st
		}
		if (q(i, k / st) >= 0.99) {
			q(i, k / st) = q(i, k / st)*ps;
		}
		arma::uvec find_qz = arma::find(q.col(k / st) == 0);
		for (int j = 0; j < find_qz.n_rows; j++) {
			//q(find_qz(j), k / st) = init_quality;
		}
		
		//std::cout << q(i, k / st) << std::endl;
		//q(i, k / st) = (1 - dk(i, 0))*(q(i, k / st)) / (1 - dk(i, 0) * q(i, k / st)); //need to change to k/st
		//std::cout <<"hist check " <<hist_check << std::endl;
	}
	
	//std::cout <<"beta " << beta << std::endl;
	
	//std::cout << "currmeascout " << current_meas_count << std::endl;
	
	arma::mat asso_data = arma::zeros(std::max(current_meas_count, find_id_size), 6);
	/*
	arma::mat unasso_at = arma::zeros(current_meas_count, 1);
	for (int i = 0; i < current_meas_count; i++) {
		unasso_at(i, 0) = i;
	}
	*/
	//std::cout << unasso_at << std::endl;
	arma::mat beta_tmp = beta;
	beta.fill(0);
	//std::cout << beta_tmp << std::endl;

	arma::mat P_tmp = P;
	P_tmp.resize(4 * find_id_size, 4);
	for (int iii = 0; iii < find_id_size; iii++) {
		arma::mat weight_sum = arma::zeros(1, 2);
		arma::mat weight_ar = arma::zeros(1, 2);
		arma::mat p_tot = P;
		if (current_meas_count != 0) {
			for (int iiii = 0; iiii < current_meas_count; iiii++) {
				arma::mat curr_meas_tmp = arma::zeros(1, 2);
				arma::mat curr_st_tmp = arma::zeros(1, 2);
				curr_meas_tmp << curr_meas(iiii * 2, 0) << curr_meas(iiii * 2 + 1, 2) << arma::endr;
				ra2xy4(curr_meas_tmp, curr_st_tmp, radi, sigv);

				//curr_st_tmp(0, 0) = curr_st_tmp(0, 0) - 500;
				//curr_st_tmp(0, 1) = -(curr_st_tmp(0, 1) - 500);

				//std::cout << curr_st_tmp << std::endl;
				weight_sum(0, 0) = weight_sum(0, 0) + beta_tmp(iii, iiii)*curr_st_tmp(0, 0);
				weight_sum(0, 1) = weight_sum(0, 1) + beta_tmp(iii, iiii)*curr_st_tmp(0, 1);
				weight_ar(0, 0) = weight_ar(0, 0) + beta_tmp(iii, iiii)*curr_meas_tmp(0, 0);
				weight_ar(0, 1) = weight_ar(0, 1) + beta_tmp(iii, iiii)*curr_meas_tmp(0, 1);
				
			}
			//std::cout <<"sum "<< weight_sum << std::endl;
			//arma::mat weight_ar = arma::zeros(1, 2);

			//weight_sum(0, 0) = weight_sum(0, 0) + 500;
			//weight_sum(0, 1) = -weight_sum(0, 1) + 500;

			//xy2ra4(weight_sum, weight_ar, radi, sigv);

			/*double testfind = beta_tmp.row(iii).index_max();
			asso_data(iii, 4) = curr_meas(testfind * 2, 0);
			asso_data(iii, 5) = curr_meas(testfind * 2 + 1, 2);*/
			//asso points
			asso_data(iii, 4) = weight_ar(0, 0);
			asso_data(iii, 5) = weight_ar(0, 1);


			for (int iiii = 0; iiii < current_meas_count; iiii++) {
				arma::mat curr_meas_tmp = arma::zeros(1, 2);
				arma::mat curr_st_tmp = arma::zeros(1, 2);
				curr_meas_tmp << curr_meas(iiii * 2, 0) << curr_meas(iiii * 2 + 1, 2) << arma::endr;
				ra2xy4(curr_meas_tmp, curr_st_tmp, radi, sigv);

				arma::mat hxi = arma::zeros(4, 4);
				arma::mat hx = arma::zeros(4, 4);
				arma::mat pxi = arma::zeros(4, 4);

				
				hxi(0, 0) = curr_st_tmp(0, 0);
				hxi(2, 2) = curr_st_tmp(0, 1);
				hx(0, 0) = weight_sum(0, 0);
				hx(2, 2) = weight_sum(0, 1);
				
				/*hxi(0, 0) = curr_meas_tmp(0, 0);
				hxi(2, 2) = curr_meas_tmp(0, 1);
				hx(0, 0) = weight_ar(0, 0);
				hx(2, 2) = weight_ar(0, 1);*/

				//arma::mat p_mid = P_tmp(arma::span(iii * 4, iii * 4 + 3), arma::span(0, 3)) + (hxi - hx)*(hxi - hx).t();
				//p_tot(arma::span(iii * 4, iii * 4 + 3), arma::span(0, 3)) = p_tot(arma::span(iii * 4, iii * 4 + 3), arma::span(0, 3)) + beta_tmp(iii, iiii)*p_mid;
				arma::mat p_mid = P_tmp(arma::span(iii * 4, iii * 4 + 3), arma::span(0, 3)) + (hxi - hx)*(hxi - hx).t();
				p_tot(arma::span(iii * 4, iii * 4 + 3), arma::span(0, 3)) = p_tot(arma::span(iii * 4, iii * 4 + 3), arma::span(0, 3)) + beta_tmp(iii, iiii)*p_mid;
			}
			//std::cout << "p pp" << P(arma::span(iii * 4, iii * 4 + 3), arma::span(0, 3)) << std::endl;
			//std::cout << "p tot" << p_tot(arma::span(iii * 4, iii * 4 + 3), arma::span(0, 3)) << std::endl;
			P = p_tot; //P=P+p_tot;

		}
		else {
			asso_data(iii, 4) = asso_data(iii, 2);
			asso_data(iii, 5) = asso_data(iii, 3);
		}
		
		//id and quailty
		asso_data(iii, 0) = t_id(iii, (k - 1) / st);
		asso_data(iii, 1) = q(iii, k / st);
		//past track point
		asso_data(iii, 2) = hzk1k(iii * 2, 0);
		asso_data(iii, 3) = hzk1k(iii * 2 + 1, 2);
		
		


		/*
		int max_at = beta_tmp.index_max();
		
		if (beta_tmp(max_at) != 0) {
			int at_row = max_at % find_id_size;
			int at_col = max_at / find_id_size;

			//std::cout << "track  " << at_row << " asso with meas " << at_col << std::endl;
			//std::cout << beta_tmp(max_at) << std::endl;
			beta(at_row, at_col) = beta_tmp(beta_tmp.index_max());

			beta_tmp.row(at_row).fill(0);
			beta_tmp.col(at_col).fill(0);

			//asso points
			asso_data(at_row, 4) = curr_meas(at_col * 2, 0);
			asso_data(at_row, 5) = curr_meas(at_col * 2 + 1, 2);
			
			//unasso meas find where
			unasso_at(arma::find(unasso_at == at_col)).fill(0);
			
			
		}
		//id and quailty
		asso_data(iii, 0) = t_id(iii, (k-1)/st);
		asso_data(iii, 1) = q(iii, k/st);
		//past track point
		asso_data(iii, 2) = hzk1k(iii * 2, 0);
		asso_data(iii, 3) = hzk1k(iii * 2 + 1, 2);
		*/
	}

	

	//unasso meas
	int unasso_num = 0;
	arma::mat unasso_meas_at;
	//unasso_meas_at.clear();
	for (int iiii = 0; iiii < current_meas_count; iiii++) {
		//std::cout << arma::sum(beta_tmp.col(iiii)) << std::endl;
		arma::mat for_insert_tmp(1, 1);
		for_insert_tmp(0, 0) = iiii;

		if (arma::sum(beta_tmp.col(iiii)) == 0) {
			unasso_num += 1;
			unasso_meas_at.resize(unasso_num, 1);
			unasso_meas_at(unasso_num - 1, 0) = iiii;
			
		}
	}
	asso_data.resize(find_id_size + unasso_num, 6);
	for (int i = 0; i < unasso_num; i++) {
		asso_data(find_id_size + i, 0) = t_id(t_id.index_max()) + i + 1;
		asso_data(find_id_size + i, 2) = curr_meas(unasso_meas_at(i) * 2, 0);
		asso_data(find_id_size + i, 3) = curr_meas(unasso_meas_at(i) * 2 + 1, 2);
		asso_data(find_id_size + i, 4) = curr_meas(unasso_meas_at(i) * 2, 0);
		asso_data(find_id_size + i, 5) = curr_meas(unasso_meas_at(i) * 2 + 1, 2);
	}




	//std::cout << asso_data << std::endl;
	int del_track_size = find_id_size;
	arma::mat q_tmp = q.col(k / st);

	for (int i = 0; i < del_track_size; i++) {
		if (asso_data(i, 1) <= term_th) {
			
			asso_data.shed_row(i);
			q_tmp.shed_row(i);
			del_track_size -= 1;
		}
	}
	q_tmp.resize(maxdetect, 1);
	q.col(k / st) = q_tmp;
	

	//for (int iiii = 0; iiii < asso_data.n_rows; iiii++) {
		int at_prow = 1000;
		arma::uvec find_unasso_track = arma::find(asso_data(arma::span(0, asso_data.n_rows - 1), arma::span(0, 5)) == 0);
		for (int iiiii = 0; iiiii < find_unasso_track.n_rows; iiiii++) {
			int at_row = find_unasso_track(iiiii) % asso_data.n_rows;
			
			if (find_unasso_track(iiiii) >= 2 * asso_data.n_rows) {
				asso_data(find_unasso_track(iiiii)) = asso_data(find_unasso_track(iiiii) - 2 * asso_data.n_rows);
				
				if (at_row != at_prow) {
					asso_data(at_row, 1) = asso_data(at_row, 1)*ps;
					at_prow = at_row;
				}
			}
		}
	//}

	data_output = asso_data;
	//std::cout << asso_data << std::endl;
}


void IPDA::nn(arma::mat meas_data,
	arma::mat::col_iterator ccita,
	arma::mat hzk1k, arma::mat t_id, arma::mat sk1, int k, double st,
	arma::mat &q, arma::mat &data_output,
	int radi, double sigv, arma::mat &P) {
	int n = k / st * 2;
	arma::mat curr_meas = arma::zeros(maxdetect * 2, 4);


	gamma = chi2inv(pg, deg);
	arma::uvec find_id_e = arma::find(t_id.col((k - 1) / st) != 0); // 0 need to change to k/st
																	//std::cout << t_id.col((k-1) / st) << std::endl;
	int find_id_size = arma::max(find_id_e) + 1;		//existing num of tracks
														//std::cout << "HERE" << std::endl;
	arma::mat likeihood = arma::zeros(find_id_size, maxdetect);
	arma::mat beta = arma::zeros(find_id_size, maxdetect);

	arma::mat dk_store_tmp = arma::zeros(maxdetect, 1);

	int current_meas_count = 0;
	if (*ccita != 0) {					//check if have measurement, by check the first row element of column
										//std::cout << "Have measurement " << k << std::endl;
		for (int i = 0; i < find_id_size; i++) { //loop for every track
												 //reset iterator
												 //std::cout << meas_data << std::endl;
			arma::mat::col_iterator cita = meas_data.begin_col(n);
			arma::mat::col_iterator cita_end = meas_data.end_col(n);
			arma::mat::col_iterator citr = meas_data.begin_col(n + 1);
			arma::mat::col_iterator citr_end = meas_data.end_col(n + 1);

			current_meas_count = 0;
			double tot_like = 0;

			for (; cita != cita_end; cita++) {	//increment for each row
				if ((*cita != 0)) {				//eliminate the reserved memory (zeros)

					arma::mat zk1(2, 2);
					zk1 << (*cita) << 0 << 0 << 0 << arma::endr
						<< 0 << 0 << (*citr) << 0 << arma::endr;
					curr_meas(current_meas_count * 2, 0) = (*cita);
					curr_meas(current_meas_count * 2 + 1, 2) = (*citr);

					double dis = calcdis(zk1, hzk1k, sk1, i, radi, sigv);	//calc distance to apply gate
					//std::cout << dis << std::endl;
					if (dis <= pow(gamma, 2)) {			//apply gate
														//measurement that pass the gate
														//std::cout << "pass gate" << std::endl;

						double like = likeicalc(sk1, zk1, hzk1k, i, radi, sigv) / pg;
						tot_like = tot_like + like;
						//std::cout << (*cita) << "  " << (*citr) << std::endl;

						//ignore too small likeihood
						if (like <= 3e-5) {
							//like = 0;
						}
						likeihood(i, current_meas_count) = like;
						//std::cout << "like  loop" << current_meas_count << std::endl;
					}
					else {
						//std::cout << "no pass" << std::endl;
					}
					current_meas_count++;
					citr++;
				}
			}
			dk_mid = pd * pg * tot_like / lambda;	//passed all gated measurement
			dk_store_tmp(i, 0) = pd * pg - dk_mid;

			//std::cout << "tot like " << tot_like << std::endl;
		}

	}
	else {
		//since no meas dk=pd pg
		for (int i = 0; i < find_id_size; i++) {
			std::cout << "No measurement " << k << std::endl;
			dk_store_tmp(i, 0) = pg * pd;
		}
		likeihood.fill(0);

	}
	dk = dk_store_tmp;		//this is delta_k
	dk.resize(find_id_size, 1);
	likeihood.resize(find_id_size, current_meas_count);
	beta.resize(find_id_size, current_meas_count);
	curr_meas.resize(current_meas_count * 2, 4);
	//likeihood.elem(arma::find_nonfinite(likeihood)).zeros();
	//std::cout << "like  " << likeihood << std::endl;
	//std::cout <<"dk "<< dk << std::endl;

	//track quality, beta prob
	for (int i = 0; i < find_id_size; i++) {
		arma::mat::col_iterator meas_check = meas_data.begin_col(n);

		if (*meas_check == 0) {
			beta(i, 0) = (1 - pd * pg) / (1 - dk(i, 0)); //dk change
		}
		else {
			for (int ii = 0; ii < current_meas_count; ii++) {
				beta(i, ii) = (pd*pg / lambda * likeihood(i, ii)) / (1 - dk(i, 0));//dkchange
																				   //beta(i, ii) = abs(beta(i, ii));
			}
		}


		//std::cout << "asso? " << beta.row(i).index_max() << std::endl;
		//asso_m(i, beta.row(i).index_max()) += 1;


		bool hist_check = false;	//reset check
		arma::mat unitmax(1, 1);
		arma::uvec find_hist = arma::find(t_id == t_id(i, (k - 1) / st)); //need ot change to k/st
		if (arma::size(find_hist) == arma::size(unitmax)) {
			hist_check = true;		//check if the track is just initialize by checking track history
		}
		if (hist_check == true) {
			q(i, k / st) = init_quality;		//need to change to k/st
		}
		else {
			q(i, k / st) = (1 - dk(i, 0))*(q(arma::max(find_hist))) / (1 - dk(i, 0) * q(arma::max(find_hist))); //need to change to k/st
		}
		if (q(i, k / st) >= 0.99) {
			q(i, k / st) = q(i, k / st)*ps;
		}
		arma::uvec find_qz = arma::find(q.col(k / st) == 0);
		for (int j = 0; j < find_qz.n_rows; j++) {
		}

	
	}

	arma::mat asso_data = arma::zeros(std::max(current_meas_count, find_id_size), 6);
	
	arma::mat unasso_at = arma::zeros(current_meas_count, 1);
	for (int i = 0; i < current_meas_count; i++) {
		unasso_at(i, 0) = i;
	}
	
	//std::cout << unasso_at << std::endl;
	arma::mat beta_tmp = beta;
	beta.fill(0);
	//std::cout << beta_tmp << std::endl;


	for (int iii = 0; iii < find_id_size; iii++) {
		if (current_meas_count != 0) {
			for (int iiii = 0; iiii < current_meas_count; iiii++) {
				arma::mat curr_meas_tmp = arma::zeros(1, 2);
				arma::mat curr_st_tmp = arma::zeros(1, 2);
				curr_meas_tmp << curr_meas(iiii * 2, 0) << curr_meas(iiii * 2 + 1, 2) << arma::endr;
				ra2xy4(curr_meas_tmp, curr_st_tmp, radi, sigv);

				

			}
			

		}
		else {
			asso_data(iii, 4) = asso_data(iii, 2);
			asso_data(iii, 5) = asso_data(iii, 3);
		}


		
		int max_at = beta_tmp.index_max();

		if (beta_tmp(max_at) != 0) {
			int at_row = max_at % find_id_size;
			int at_col = max_at / find_id_size;

			//std::cout << "track  " << at_row << " asso with meas " << at_col << std::endl;
			//std::cout << beta_tmp(max_at) << std::endl;
			beta(at_row, at_col) = beta_tmp(beta_tmp.index_max());

			beta_tmp.row(at_row).fill(0);
			beta_tmp.col(at_col).fill(0);

			//asso points
			asso_data(at_row, 4) = curr_meas(at_col * 2, 0);
			asso_data(at_row, 5) = curr_meas(at_col * 2 + 1, 2);

			//unasso meas find where
			unasso_at(arma::find(unasso_at == at_col)).fill(0);


		}
		//id and quailty
		asso_data(iii, 0) = t_id(iii, (k-1)/st);
		asso_data(iii, 1) = q(iii, k/st);
		//past track point
		asso_data(iii, 2) = hzk1k(iii * 2, 0);
		asso_data(iii, 3) = hzk1k(iii * 2 + 1, 2);
		
	}


	//std::cout << unasso_at << std::endl;
	//unasso meas
	arma::uvec unasso_tot = arma::find(unasso_at != 0);
	int unasso_num = unasso_tot.n_rows;
	asso_data.resize(find_id_size + unasso_num, 6);
	for (int i = 0; i < unasso_num; i++) {
		asso_data(find_id_size + i, 0) = t_id(t_id.index_max()) + i + 1;
		asso_data(find_id_size + i, 2) = curr_meas(unasso_tot(i) * 2, 0);
		asso_data(find_id_size + i, 3) = curr_meas(unasso_tot(i) * 2 + 1, 2);
		asso_data(find_id_size + i, 4) = curr_meas(unasso_tot(i) * 2, 0);
		asso_data(find_id_size + i, 5) = curr_meas(unasso_tot(i) * 2 + 1, 2);
	}




	//std::cout << asso_data << std::endl;
	int del_track_size = find_id_size;
	arma::mat q_tmp = q.col(k / st);

	for (int i = 0; i < del_track_size; i++) {
		if (asso_data(i, 1) <= term_th) {

			asso_data.shed_row(i);
			q_tmp.shed_row(i);
			del_track_size -= 1;
		}
	}
	q_tmp.resize(maxdetect, 1);
	q.col(k / st) = q_tmp;


	//for (int iiii = 0; iiii < asso_data.n_rows; iiii++) {
	int at_prow = 1000;
	arma::uvec find_unasso_track = arma::find(asso_data(arma::span(0, asso_data.n_rows - 1), arma::span(0, 5)) == 0);
	for (int iiiii = 0; iiiii < find_unasso_track.n_rows; iiiii++) {
		int at_row = find_unasso_track(iiiii) % asso_data.n_rows;

		if (find_unasso_track(iiiii) >= 2 * asso_data.n_rows) {
			asso_data(find_unasso_track(iiiii)) = asso_data(find_unasso_track(iiiii) - 2 * asso_data.n_rows);

			if (at_row != at_prow) {
				asso_data(at_row, 1) = asso_data(at_row, 1)*ps;
				at_prow = at_row;
			}
		}
	}
	//}

	data_output = asso_data;
	//std::cout << asso_data << std::endl;
}