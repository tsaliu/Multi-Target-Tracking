#pragma once
#pragma once
#include <armadillo>
#include <time.h>

void randseed_set3() {
	arma::arma_rng::seed_type seed1 = arma::arma_rng::seed_type(0);
	arma::arma_rng::seed_type seed2 = arma::arma_rng::seed_type(0);
	arma::arma_rng::seed_type seed3 = arma::arma_rng::seed_type(0);
	arma::arma_rng::seed_type seed4 = arma::arma_rng::seed_type(0);

	double seconds;
	time_t timer = std::time(nullptr);
	struct tm y2k = { 0 };
	y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
	y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;
	seconds = difftime(timer, mktime(&y2k));
	time(&timer);
	seed3 = static_cast<arma::arma_rng::seed_type>(seconds);

	seed4 = static_cast<arma::arma_rng::seed_type>(std::time(NULL) & 0xFFFF);
	arma::arma_rng::seed_type seed5 = arma::arma_rng::seed_type(0);
	arma::arma_rng::set_seed(seed1 + seed2 + seed3 + seed4 + seed5);
}