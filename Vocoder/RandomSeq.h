/*
 * RandomSeq.h
 *
 *  Created on: Feb 15, 2013
 *      Author: yuchen
 */

#ifndef RANDOMSEQ_H_
#define RANDOMSEQ_H_

#include "Filter.h"
#include <cmath>
template<class T>
class RandomSeq {
	const double mod = 4294967296.0;
	const long mul = 663608941;
	const long add = 1023812313;

private:
	unsigned seed;
	T se;
	Filter<T> h;
public:
	RandomSeq() {
		seed = 1;
		se = 1;
	}
	RandomSeq(T d, const Filter<T> & f) {
		seed = 1;
		se =d;
		h = f;
	}
	/*
	 * need support for type_cast from double to T
	 */
	T gauss() {
		double x1,x2;

		seed *= mul;
		seed = (seed ^ 0x55555555) +add;
		x1 = (seed!=0) ? seed/mod : 1/mod;
		seed *= mul;
		seed = (seed ^ 0x33333333) +add;
		x2 = seed/mod;
		T norm_randn = static_cast<T>(sqrt(-2*log(x1))*cos(PI*x2*2));
		return norm_randn* se;
	}

	vector<T> getseq(int n) {
		vector<T> wn, ret;

		wn.resize(n, 0);
		for (int i=0; i<n; i++)
			wn[i] = gauss();
		h.process(wn, ret);
		return ret;
	}
};


#endif /* RANDOMSEQ_H_ */
