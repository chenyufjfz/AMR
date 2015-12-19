/*
 * Resample.h
 *
 *  Created on: Feb 23, 2013
 *      Author: yuchen
 */

#ifndef RESAMPLE_H_
#define RESAMPLE_H_
#include <vector>
using namespace std;
#include "public.h"
#include "Process.h"

template<class T>
class Resample : public Process<T>{
public:
	enum ResampleType {
		COMBINE,
		CLASSIC
	};

private:
	vector<T> h;
	int up, down;
	unsigned h_size;
	unsigned long long prev_grain;
	ResampleType type;
	Resample<T> * r1, *r2;
	vector<T> buf_in;
public:
	/*
	 * hh is low pass filter's coefficient parameter
	 */
	Resample(T * hh, unsigned size, int u, int d) {
		h.assign(hh, hh+size);
		h_size = size;
		up = u;
		down = d;
		prev_grain =0;
		type = CLASSIC;
		r1 = NULL;
		r2 = NULL;
	}
	Resample(vector <T> hh, int u, int d) {
		h=hh;
		h_size = hh.size();
		up = u;
		down = d;
		prev_grain =0;
		type = CLASSIC;
		r1 = NULL;
		r2 = NULL;
	}
	Resample(Resample<T> *rs1, Resample<T> *rs2) {
		r1 = rs1;
		r2 = rs2;
		h.clear();
		h_size = 0;
		prev_grain =0;
		up = rs1->up * rs2->up;
		down = rs1->down * rs2->down;
		type = COMBINE;
	}
	~Resample() {

	}
	/*
	 * resample process date like following
	 * ^up --> lowpassfilter  --> |down
	 */
	void process(const vector<T> &x, vector<T> &y) {
		unsigned i, j, k, x_align, xx_size;
		unsigned long long prev_k;
		if (type==COMBINE) {
			vector<T> m;
			r1->process(x, m);
			r2->process(m, y);
			return;
		}
		if (type==CLASSIC) {
			vector <T> xx;
			if (buf_in.empty())
				buf_in.assign(h_size/(up*2), 0);
			xx = buf_in;
			xx.insert(xx.end(), x.begin(), x.end());
			xx_size =xx.size();	//use xx_size to save time for xx.size()
			y.assign(xx_size*up/down+20, 0); //+20 is for leave enough space
			prev_k = (prev_grain+up-1) /up;
			for (i=0; ; i++) {
				j = (up -prev_grain %up) %up; //up-(i*down)%up)%up, prev_grain=i*down
				x_align = (prev_grain+up-1) /up -prev_k; //(i*down-1)/up+1
				for (k=x_align; j<h_size && k<xx_size; j+=up,k++)
					y[i] += h[j] * xx[k];

				if (k >= xx_size) {
					buf_in.assign(xx.begin()+ x_align, xx.end());
					y.resize(i);
					return;
				} else
					prev_grain +=down;
			}
		}
	}

};

#endif /* RESAMPLE_H_ */
