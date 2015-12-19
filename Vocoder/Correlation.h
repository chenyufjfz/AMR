/*
 * Correlation.h
 *
 *  Created on: Mar 16, 2013
 *      Author: yuchen
 */

#ifndef CORRELATION_H_
#define CORRELATION_H_

#include <vector>
using namespace std;

template <class T>
class Correlation {
private:
	int degree;
	vector<T> r;
public:
	Correlation(int n) {
		degree = n;
		r.assign(2*n+1, 0);
	}
	vector<T> get_corr() {
		return r;
	}
	void clear_corr() {
		r.assign(2*degree+1, 0);
	}
	void relate(const vector<T> &x0, const vector<T> &x1) {

	}
};


#endif /* CORRELATION_H_ */
