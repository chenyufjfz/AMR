/*
 * filter.h
 *
 *  Created on: Feb 12, 2013
 *      Author: yuchen
 */

#ifndef FILTER_H_
#define FILTER_H_

#include "Polynomial.h"
#include "Process.h"

template <class T>
class Filter : public Process<T> {
private:
	Polynomial<T> numerator, denominator;
	vector<T> buf_in, buf_out;
	/*
	 * Pass equation: numerator(z) / denominator(z), where
	 * numerator(z) = numerator[0] + numerator[1]*z^(-1) + numerator[2]*z^(-2)...;
	 * denominator(z) = denominator[0] + denominator[1]*z^(-1) + denominator[2]*z^(-2)...;
	 */
public:
	enum FilterType {
		ALL_PASS,
		IIR,
		FIR,
	};
public:
	Filter() {
		numerator = Polynomial<T> (1);
		denominator = Polynomial<T> (1);
	}

	Filter(FilterType type, const vector<T> & para_zi, const vector<T> & para_mu) {
		switch (type) {
		case IIR:
			denominator = Polynomial<T> (para_mu);
			numerator = Polynomial<T> (para_zi);
			this->normalize();
			break;
		case ALL_PASS:
			throw MyException("Filter create error, not need two para for ALL_PASS");
			break;
		case FIR:
			throw MyException("Filter create error, not need two para for FIR");
			break;
		}
	}
	Filter(FilterType type, const Polynomial<T> & para_zi, const Polynomial<T> & para_mu) {
		switch (type) {
		case IIR:
			denominator = para_mu;
			numerator = para_zi;
			this->normalize();
			break;
		case ALL_PASS:
			throw MyException("Filter create error, not need two para for ALL_PASS");
			break;
		case FIR:
			throw MyException("Filter create error, not need two para for FIR");
			break;
		}
	}
	Filter(FilterType type, const vector<T> &para) {
		switch (type) {
		case ALL_PASS:
			denominator = Polynomial<T> (para);
			numerator = denominator.reverse();
			break;
		case FIR:
			numerator = Polynomial<T> (para);
			denominator = Polynomial<T> (1);
			break;
		case IIR:
			throw MyException("Filter create error, need two para for IIR");
			break;
		}
	}
	Filter(FilterType type, const Polynomial<T> &para) {
		switch (type) {
		case ALL_PASS:
			denominator = para;
			numerator = denominator.reverse();
			break;
		case FIR:
			numerator = para;
			denominator = Polynomial<T> (1);
			break;
		case IIR:
			throw MyException("Filter create error, need two para for IIR");
			break;
		}
	}
	//change filter parameter dynamically without clearing buf
	void dynamic_change(const vector<T> & para0, const vector<T> & para1)
	{
		ASSERT(para0.size() == numerator.size() && para1.size() == denominator.size());
		numerator = Polynomial<T> (para0);
		denominator = Polynomial<T> (para1);
		this->normalize();
	}

	//change filter parameter dynamically without clearing buf
	void dynamic_change(const Polynomial<T> & para0, const Polynomial<T> & para1)
	{
		ASSERT(para0.size() == numerator.size() && para1.size() == denominator.size());
		numerator = para0;
		denominator = para1;
		this->normalize();
	}

	~Filter() {

	}
	void normalize() {
		while (denominator[0]==0 && numerator[0]==0) {
			denominator >>= 1;
			numerator >>=1;
			if (denominator.size()==1)
				break;

		}
		if (denominator[0]==0 && denominator.size()==1)
			throw MyException("denominator is zero");
		buf_in.clear();
		buf_out.clear();
	}
	Polynomial<T> get_numerator() {
		return numerator;
	}
	Polynomial<T> get_denominator() {
		return denominator;
	}
	Filter <T> operator - () {
		Filter <T> temp;
		temp.numerator = -this->numerator;
		temp.denominator = this->denominator;
		this->normalize();
		return temp;
	}
	Filter <T>& operator *=(const T a) {
		numerator = numerator * a;
		this->normalize();
		return *this;
	}
	Filter <T>& operator /=(const T a) {
		denominator = denominator * a;
		this->normalize();
		return *this;
	}
	Filter <T>& operator +=(const Filter<T>& a) {
		numerator = numerator *a.denominator + denominator *a.numerator;
		denominator *= a.denominator;
		this->normalize();
		return *this;
	}
	Filter <T>& operator -=(const Filter<T>& a) {
		numerator = numerator *a.denominator - denominator *a.numerator;
		denominator *= a.denominator;
		this->normalize();
		return *this;
	}
	Filter <T>& operator *=(const Filter<T>& a) {
		numerator *= a.numerator;
		denominator *= a.denominator;
		this->normalize();
		return *this;
	}
	Filter <T>& operator /=(const Filter<T>& a) {
		numerator = numerator * a.denominator;
		denominator *= a.numerator;
		this->normalize();
		return *this;
	}
	void clear_buf() {
		buf_in.clear();
		buf_out.clear();
	}

	/*
	 * prefill zero for x and y
	 * xx:  0, 0,...0, x[0], x[1], x[2],...., x[n-1]
	 *      |<- n1 ->| |<-       n               ->|
	 * yy:  0, 0,...0, y[0], y[1], y[2],...., y[n-1]
	 *      |<- n2 ->| |<-       n               ->|
	 */
	void process(const vector<T> &x, vector<T> &y) {
		int n, n1, n2, i, j;
		vector <T> xx, yy;

		if (denominator[0] ==0)
			throw MyException("filter is not a casual filter");
		n = x.size();
		n1 = numerator.degree();
		n2 = denominator.degree();
		if (buf_in.empty()) {
			buf_in.assign(n1, 0);
			buf_out.assign(n2, 0);
		}
		xx = buf_in;
		xx.insert(xx.end(), x.begin(), x.end());
		yy = buf_out;
		yy.resize(n+n2);
		for (i=0; i<n; i++) {
			T e;
			e=0;
			for (j=0; j<=n1; j++)
				e += xx[i+j] * numerator[n1-j];
			for (j=0; j<n2; j++)
				e -= yy[i+j] * denominator[n2-j];
			yy[n2+i] = e / denominator[0];
		}
		y.assign(yy.begin()+n2, yy.end());
		buf_out.assign(yy.begin()+n, yy.end());
		ASSERT(buf_out.size()==(unsigned) n2);
		buf_in.assign(xx.begin()+n, xx.end());
		ASSERT(buf_in.size()==(unsigned) n1);
	}
};

template <class T>
Filter<T> operator + (const Filter<T>& a, const Filter<T>& b)
{
	Filter <T> temp;
	temp =a;
	temp +=b;
	return temp;
}
template <class T>
Filter<T> operator - (const Filter<T>& a, const Filter<T>& b)
{
	Filter <T> temp;
	temp =a;
	temp -=b;
	return temp;
}
template <class T>
Filter<T> operator * (const Filter<T>& a, const Filter<T>& b)
{
	Filter <T> temp;
	temp =a;
	temp *=b;
	return temp;
}
template <class T>
Filter<T> operator / (const Filter<T>& a, const Filter<T>& b)
{
	Filter <T> temp;
	temp =a;
	temp /=b;
	return temp;
}
template <class T>
Filter<T> operator * (const Filter<T>& a, const T& b)
{
	Filter <T> temp;
	temp =a;
	temp *= b;
	return temp;
}
template <class T>
Filter<T> operator * (const T & a, const Filter<T>& b)
{
	Filter <T> temp;
	temp =b;
	temp *= a;
	return temp;
}
template <class T>
Filter<T> operator / (const Filter<T>& a, const T& b)
{
	Filter <T> temp;
	temp =a;
	temp /= b;
	return temp;
}
template <class T>
Filter<T> operator / (const T & a, const Filter<T>& b)
{
	Filter <T> temp;
	temp =b;
	temp /= a;
	return temp;
}
template <class T>
ostream& operator << (ostream& ostrm, Filter<T>& m) {
	ostrm <<"numerator="<<m.get_numerator()<<'\n';
	ostrm <<"denominator="<<m.get_denominator()<<'\n';
	return ostrm;
}

template <class T>
class FilterN : public Process<T> {
private:
	int up;
	Polynomial <T> numerator, denominator;
	vector <T> buf_in, buf_out;
public:
	enum FilterNType {
		IIR,
		FIR,
	};
	FilterN(FilterNType type, int n, const vector<T> & para0, const vector<T> & para1) {
		switch (type) {
		case IIR:
			ASSERT(n>0);
			up = n;
			denominator = Polynomial<T> (para1);
			numerator = Polynomial<T> (para0);
			this->normalize();
			break;
		case FIR:
			throw MyException("Filter create error, not need two para for FIR");
			break;
		}
	}
	FilterN(FilterNType type, int n, const Polynomial<T> & para0, const Polynomial<T> & para1) {
		switch (type) {
		case IIR:
			ASSERT(n>0);
			up = n;
			denominator = para1;
			numerator = para0;
			this->normalize();
			break;
		case FIR:
			throw MyException("Filter create error, not need two para for FIR");
			break;
		}
	}
	FilterN(FilterNType type, int n, const vector<T> &para) {
		switch (type) {
		case FIR:
			ASSERT(n>0);
			up = n;
			numerator = Polynomial<T> (para);
			denominator = Polynomial<T> (1);
			break;
		case IIR:
			throw MyException("Filter create error, need two para for IIR");
			break;
		}
	}
	FilterN(FilterNType type, int n, const Polynomial<T> &para) {
		switch (type) {
		case FIR:
			ASSERT(n>0);
			up = n;
			numerator = para;
			denominator = Polynomial<T> (1);
			break;
		case IIR:
			throw MyException("Filter create error, need two para for IIR");
			break;
		}
	}
	FilterN(FilterNType type, int n, const Filter<T> &para) {
		up = n;
		numerator = para.get_numerator();
		denominator = para.get_denominator();
	}
	~FilterN() {

	}
	void normalize() {
		while (denominator[0]==0 && numerator[0]==0) {
			denominator >>= 1;
			numerator >>=1;
			if (denominator.size()==1)
				break;
		}
		if (denominator[0]==0 && denominator.size()==1)
			throw MyException("denominator is zero");
		buf_in.clear();
		buf_out.clear();
	}
	Polynomial<T> get_numerator() {
		return numerator;
	}
	Polynomial<T> get_denominator() {
		return denominator;
	}
	int get_up() {
		return up;
	}
	/*
	 * prefill zero for x and y
	 * xx:  0, 0,... 0, x[0], x[1], x[2],...., x[n-1]
	 *      |<-n1*up->| |<-       n               ->|
	 * yy:  0, 0,... 0, y[0], y[1], y[2],...., y[n-1]
	 *      |<-n2*up->| |<-       n               ->|
	 */
	void process(const vector<T> &x, vector<T> &y) {
		int n, n1, n2, i, j, k;
		vector <T> xx, yy;

		if (denominator[0] ==0)
			throw MyException("filter is not a casual filter");
		n = x.size();
		n1 = numerator.degree();
		n2 = denominator.degree();
		if (buf_in.empty()) {
			buf_in.assign(n1*up, 0);
			buf_out.assign(n2*up, 0);
		}
		xx = buf_in;
		xx.insert(xx.end(), x.begin(), x.end());
		yy = buf_out;
		yy.resize(n+n2*up);
		for (i=0; i<n; i++) {
			T e;
			e=0;
			for (j=0, k=0; j<=n1; j++, k+=up)
				e += xx[i+k] * numerator[n1-j];
			for (j=0, k=0; j<n2; j++, k+=up)
				e -= yy[i+k] * denominator[n2-j];
			yy[n2*up+i] = e / denominator[0];
		}
		y.assign(yy.begin()+n2*up, yy.end());
		buf_out.assign(yy.begin()+n, yy.end());
		ASSERT(buf_out.size()==(unsigned) n2*up);
		buf_in.assign(xx.begin()+n, xx.end());
		ASSERT(buf_in.size()==(unsigned) n1*up);
	}
};
#endif /* FILTER_H_ */
