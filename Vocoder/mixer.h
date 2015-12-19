/*
 * mixer.h
 *
 *  Created on: Mar 16, 2013
 *      Author: yuchen
 */

#ifndef MIXER_H_
#define MIXER_H_

#include "public.h"

template <class T>
class Mixer {
private:
	T c0, c1, c2;
	int input;
public:
	Mixer(T a0) {
		c0 = a0;
		input = 1;
	}
	Mixer(T a0, T a1) {
		c0 = a0;
		c1 = a1;
		input = 2;
	}
	Mixer(T a0, T a1, T a2) {
		c0 = a0;
		c1 = a1;
		c2 = a2;
		input = 3;
	}
	void mix(const vector<T> & x, vector<T> & y) {
		ASSERT(input==1);
		y.resize(x.size());
		for (unsigned i=0; i<x.size(); i++)
			y[i] = x[i] * c0;
	}
	void mix(const vector<T> & x0, const vector<T> & x1, vector<T> & y) {
		ASSERT(input==2 && x0.size()==x1.size());
		y.resize(x0.size());
		for (unsigned i=0; i<x0.size(); i++)
			y[i] = x0[i] * c0 + x1[i] * c1;
	}
	void mix(const vector<T> & x0, const vector<T> & x1, const vector<T> & x2, vector<T> & y) {
		ASSERT(input==3 && x0.size()==x1.size() && x0.size()==x2.size());
		y.resize(x0.size());
		for (unsigned i=0; i<x0.size(); i++)
			y[i] = x0[i] * c0 + x1[i] * c1 + x2[i] * c2;
	}
};


#endif /* MIXER_H_ */
