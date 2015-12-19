/*
 * polynomial.h
 *
 *  Created on: Feb 12, 2013
 *      Author: yuchen
 */

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include "public.h"
using namespace std;

template <class T>
class Polynomial
{
private:
	vector<T> coef;
public:
   // construct:
	Polynomial() {
		coef.resize(1,0);
	}
	Polynomial(const T value) {
		coef.resize(1,value);
	}

	Polynomial(const vector<T> & c) {
		coef = c;
	}

	// access:
	unsigned size()const {
		return coef.size();
	}
	vector<T> get_vec() {
		return coef;
	}
	int degree()const {
		return coef.size() -1;
	}
	T& operator[](int i) {
		return coef[i];
	}
	T operator[](int i) const{
		return coef[i];
	}
	void normalize() {
		while (coef.back()==0 && coef.size()>1)
			coef.pop_back();
	}
	// operators:
	Polynomial <T>& operator +=(const T value) {
	   if (size()==0)
		   throw MyException("polynomial operator += empty polynomial");
	   coef[0] += value;
	   return *this;
	}

	Polynomial <T>& operator -=(const T value) {
	   if (size()==0)
		   throw MyException("polynomial operator += empty polynomial");
	   coef[0] -= value;
	   return *this;
	}

	Polynomial <T>& operator *=(const T value) {
	   if (size()==0)
		   throw MyException("polynomial operator += empty polynomial");
	   for (int i=0; i<size(); i++)
		   coef[i] *= value;
	   this->normalize();
	   return *this;
	}

	Polynomial <T> operator - () {
		Polynomial <T> temp;
		temp.coef = coef;
		for (int i=0; i<size(); i++)
			temp.coef[i] = -temp.coef[i];
		return temp;
	}

	Polynomial <T>& operator +=(const Polynomial<T>& a) {
		unsigned i;
		unsigned a_size=a.size();
		if (a_size > size())
		   coef.resize(a_size, 0);
		for (i=0; i<coef.size(); i++)
		   if (i<a_size)
			   coef[i] += a.coef[i];
		this->normalize();
		return *this;
	}

	Polynomial <T>& operator -=(const Polynomial<T>& a) {
		int i;
		int a_size=a.size();
		if (a_size > size())
		   coef.resize(a_size, 0);
		for (i=0; i<coef.size(); i++)
		   if (i<a_size)
			   coef[i] -= a.coef[i];
		this->normalize();
		return *this;
	}

	Polynomial <T>& operator *=(const Polynomial<T>& a) {
		vector<T> result;
		int a_size=a.size();

		if (a_size==0)
		   throw MyException("polynomial operator *= wrong input");
		result.assign(size() +a_size -1, 0);
		for (int i=0; i<(int)result.size(); i++)
		   for (int j=0; j<(int)size(); j++)
			   if (i-j<a_size && i-j>=0)
				   result[i] += coef[j] * a.coef[i-j];
		coef = result;
		this->normalize();
		return * this;
	}

	Polynomial <T> reverse() {
		Polynomial <T> temp;
		temp.coef.resize(this->size());
		reverse_copy(coef.begin(), coef.end(), temp.coef.begin());
		temp.normalize();
		return temp;
	}
	Polynomial <T>& operator >>=(int a) {
		if (a>=(int)coef.size())
			coef.assign(1, 0);
		else
			coef.erase(coef.begin(), coef.begin()+a);
		return * this;
	}
	Polynomial <T>& operator <<=(int a) {
		coef.insert(coef.begin(), a, 0);
		return * this;
	}
	Polynomial <T>& operator =(const vector<T> & a)
	{
		coef = a;
		return * this;
	}
};

template <class T>
Polynomial<T> operator + (const Polynomial<T>& a, const Polynomial<T>& b)
{
	Polynomial <T> temp;
	temp =a;
	temp +=b;
	return temp;
}
template <class T>
Polynomial<T> operator - (const Polynomial<T>& a, const Polynomial<T>& b)
{
	Polynomial <T> temp;
	temp =a;
	temp -=b;
	return temp;
}
template <class T>
Polynomial<T> operator * (const Polynomial<T>& a, const Polynomial<T>& b)
{
	Polynomial <T> temp;
	temp =a;
	temp *=b;
	return temp;
}

template <class T>
Polynomial<T> operator + (const Polynomial<T>& a, const T& b)
{
	Polynomial <T> temp;
	temp =a;
	temp += b;
	return temp;
}

template <class T>
Polynomial<T> operator - (const Polynomial<T>& a, const T& b)
{
	Polynomial <T> temp;
	temp =a;
	temp -= b;
	return temp;
}
template <class T>
Polynomial<T> operator * (const Polynomial<T>& a, const T& b)
{
	Polynomial <T> temp;
	temp =a;
	temp *= b;
	return temp;
}
template <class T>
Polynomial<T> operator + (const T& a, const Polynomial<T>& b)
{
	Polynomial <T> temp;
	temp =b;
	temp += a;
	return temp;
}
template <class T>
Polynomial<T> operator - (const T& a, const Polynomial<T>& b)
{
	Polynomial <T> temp;
	temp =-b;
	temp += a;
	return temp;
}
template <class T>
Polynomial<T> operator * (const T& a, const Polynomial<T>& b)
{
	Polynomial <T> temp;
	temp =b;
	temp *= a;
	return temp;
}
template <class T>
Polynomial<T> operator << (const Polynomial<T>& a, const T& b)
{
	Polynomial <T> temp;
	temp =a;
	temp <<= b;
	return temp;
}
template <class T>
Polynomial<T> operator >> (const Polynomial<T>& a, const T& b)
{
	Polynomial <T> temp;
	temp =a;
	temp >>= b;
	return temp;
}
template <class T>
ostream& operator << (ostream& ostrm, const Polynomial<T>& m) {
	ostrm<<'(';
	for (unsigned i=0; i<m.size(); i++) {
		ostrm<<m[i];
		if (i+1==m.size())
			ostrm<<')';
		else
			ostrm<<',';
	}
	return ostrm;
}
#endif /* POLYNOMIAL_H_ */
