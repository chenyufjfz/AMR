/*
 * FFT.h
 *
 *  Created on: Dec 4, 2015
 *      Author: yuchen
 */

#ifndef FFT_H_
#define FFT_H_
#include "public.h"
#include <vector>
#include <complex>
using namespace std;

class FFT {
protected:
	unsigned revert_coef[32];
	unsigned revert(unsigned idx, unsigned m) {
		int ret=0;
		for (unsigned i=0; i<m; i++)
			ret += (idx & (1<<i)) ? revert_coef[m-1-i] : 0;
		return ret;
	}

public:
	typedef complex<double> ComplexDouble;

	static void convert(const vector<int> &x, vector<ComplexDouble> &y)
	{
		y.resize(x.size());
		for (unsigned i=0; i<x.size(); i++)
			y[i] = ComplexDouble(x[i], 0);
	}

	static void convert(const vector<ComplexDouble> &x, vector<int> &y)
	{
		y.resize(x.size());
		for (unsigned i=0; i<x.size(); i++) {
			ASSERT(abs(x[i].imag()) <0.001 && abs(x[i].real())<(1<<30));
			y[i] = (int) x[i].real();
		}
	}

	static void convert(const vector<double> &x, vector<ComplexDouble> &y)
	{
		y.resize(x.size());
		for (unsigned i=0; i<x.size(); i++)
			y[i] = ComplexDouble(x[i], 0);
	}

	static void convert(const vector<ComplexDouble> &x, vector<double> &y)
	{
		y.resize(x.size());
		for (unsigned i=0; i<x.size(); i++) {
			ASSERT(abs(x[i].imag()) <0.001 && abs(x[i].real())<(1<<30));
			y[i] = x[i].real();
		}
	}


	FFT() {
		for (int i=0; i<32; i++)
			revert_coef[i] = 1 << i;
	}

	//x and f can't be same
	static void dft(vector<ComplexDouble> &x, vector<ComplexDouble> &f, bool invert=false)
	{
		int i, j, n;
		vector <ComplexDouble> w;
		double sqrt_n;
		n = x.size();
		w.resize(n);
		f.resize(n);
		sqrt_n = sqrt(n);
		for (i=0; i<n; i++)
			w[i] = invert ? polar(1.0, 2*PI*i/n) :polar(1.0, -2*PI*i/n);
		for (i=0; i<n; i++) {
			f[i] = ComplexDouble(0, 0);
			for (j=0; j<n; j++)
				f[i] += x[j] * w[(i*j) %n];
			f[i] = f[i] / sqrt_n;
		}
	}
	//x and y can't be same
	static void dft(vector<double> &x, vector<ComplexDouble> &y, bool invert=false)
	{
		vector<ComplexDouble> s;

		convert(x, s);
		dft(s, y, invert);
	}

	static void dft(vector<int> &x, vector<ComplexDouble> &y, bool invert=false)
	{
		vector<ComplexDouble> s;

		convert(x, s);
		dft(s, y, invert);
	}

	void fft2(vector<ComplexDouble> &x, vector<ComplexDouble> &y, bool invert=false)
	{
		int i,j,k,m,n,mask;
		vector <ComplexDouble> w, f1;

		n =x.size();
		for (m=0,i=n; i!=1; i=i>>1, m++)
			if (i%2!=0)
				throw MyException("fft2 input is not 2^m");
		y=x;
		w.resize(n/2);
		for (i=0; i<n/2; i++)
			w[i] = invert ? polar(1.0, 2*PI*i/n) :polar(1.0, -2*PI*i/n);
		for (i=m-1; i>=0; i--) {
			mask = 1<<i;
			for (j=0; j<n; j++)
				if (!(j & mask)) {
					ComplexDouble scale;
					k = j+mask;
					scale = w[revert(k & ~(mask*2-1), m)<<i] * y[k];
					y[k] = y[j] - scale;
					y[j] = y[j] + scale;
				}
		}
		f1 =y;
		double sqrt_n = sqrt(n);
		for (j=0; j<n; j++)
			y[revert(j,m)] = f1[j]/sqrt_n;

	}

	void fft2(vector<double> &x, vector<ComplexDouble> &y, bool invert=false)
	{
		vector<ComplexDouble> s;

		convert(x, s);
		fft2(s, y, invert);
	}

	void fft2(vector<int> &x, vector<ComplexDouble> &y, bool invert=false)
	{
		vector<ComplexDouble> s;

		convert(x, s);
		fft2(s, y, invert);
	}

	static double energy(vector<ComplexDouble> &x)
	{
		double ret =0;
		for (unsigned i=0; i<x.size(); i++)
			ret += norm(x[i]);
		return ret;
	}
	static double energy(vector<double> &x)
	{
		double ret =0;
		for (unsigned i=0; i<x.size(); i++)
			ret += x[i] * x[i];
		return ret;
	}
	static void self_check()
	{
		int i, j, n=1024;
		vector<ComplexDouble> s, f, f1;
		srand(2);
		s.resize(n);
		FFT myfft;
		for (i=0; i<100; i++) {
			for (j=0; j<n; j++)
				s[j] = ComplexDouble(rand()/1000, rand()/1000);
			myfft.fft2(s, f,false);
			FFT::dft(s,f1,false);
			for (j=0; j<n; j++)
				if (norm(f[j]-f1[j])>0.000001)
					cout<<"error\n";
		}
		cout<<"test finished";
	}

	//y and x can be same
	static void win_cos(vector<double> &x, vector <double> &y, int win_len)
	{
		int i, n;
		n = x.size();
		y.resize(n);
		if (win_len==1)
			y[0] = x[0] * 0.5;
		else
			for (i=0; i<win_len; i++)
				y[i] = x[i] * (1-cos(i *2 *PI /(win_len*2-2)))/2;
		for (i=win_len; i<n-win_len; i++)
			y[i] = x[i];
		if (win_len==1)
			y[n-1] = x[n-1] * 0.5;
		else
			for (i=n-win_len; i<n; i++)
				y[i] = x[i] * (1+cos((i-n+win_len) *2 *PI /(win_len*2-2)))/2;
	}

	//y and x can be same
	static void win_ladder(vector<double> &x, vector <double> &y, int win_len)
	{
		int i, n;
		n = x.size();
		y.resize(n);
		if (win_len==1)
			y[0] = x[0] * 0.5;
		else
			for (i=0; i<win_len; i++)
				y[i] = x[i] * i / (win_len-1);
		for (i=win_len; i<n-win_len; i++)
			y[i] = x[i];
		if (win_len==1)
			y[n-1] = x[n-1] * 0.5;
		else
			for (i=n-win_len; i<n; i++)
				y[i] = x[i] * (n-1-i) / (win_len-1);
	}
};

#endif /* FFT_H_ */
