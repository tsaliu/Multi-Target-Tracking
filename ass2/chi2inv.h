#pragma once

double chi2inv(double pg, int deg) {
	if (pg == 0.99 && deg == 1) {
		return 6.6348966;
	}
	if (pg == 0.99 && deg == 2) {
		return 9.21034037;
	}
	if (pg == 0.95 && deg == 1) {
		return 3.84145882;
	}
	if (pg == 0.95 && deg == 2) {
		return 5.99146454;
	}
}