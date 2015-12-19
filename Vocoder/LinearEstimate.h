/*
 * LinearEstimate.h
 *
 *  Created on: Feb 15, 2013
 *      Author: yuchen
 */

#ifndef LINEARESTIMATE_H_
#define LINEARESTIMATE_H_

#define _NO_NAMESPACE
#include "matrix.h"
#include "public.h"
#include "fitsutil.h"
#include <math.h>
#include "FFT.h"
#include "Resample.h"
#include "delay.h"
#include "mixer.h"
#define USE_LSP 1

template <class T>
class LinearEstimate {
public:
	vector<T> r;
public:
	enum AR_METHODS {
		AR_SELF_COR,
		AR_FORWARD_MMSE,
		AR_BACKWARD_MMSE,
		AR_BURG
	};

	static void r_az(const vector<T> &r, vector<T> &a) {
		int p=r.size();
		vector<T> new_a(p, 0);
		a = new_a;
		for (int i=0; i<p; i++) {
			for (int j=0; j<i; j++)
				new_a[j] = a[j] - r[i]*a[i-j-1];
			new_a[i] = -r[i];
			a = new_a;
		}
	}

	/*
	 * input x, p
	 * output a,
	 * methods:
	 * AR_SELF_COR:
	 *   R[0] + a[0]*R[1] + a[1]*R[2]... + a[p-1] * R[p]  =e^2
	 *   R[1] + a[0]*R[0] + a[1]*R[1]... + a[p-1]* R[p-1] =0
	 *   R[2] + a[0]*R[1] + a[1]*R[0]... + a[p-1]* R[p-2] =0
	 * AR_FORWARD_MMSE
	 *   x[i] = -a[0]*x[i-1]-a[1]*x[i-2]-a[2]*x[i-3]...-a[p-1]*x[i-p]
	 * AR_BACKWARD_MMSE
	 *   x[i] = -a[0]*x[i+1]-a[1]*x[i+2]-a[2]*x[i+3]...-a[p-1]*x[i+p]
	 * After ar estimation
	 *   x * (1+a[0]z^-1 + a[1]z^-2 + a[2]z^-3....) is white noise
	 */
	T ar_est(AR_METHODS methods, const vector<T> &x, vector<T> &a, int p) {
		int n = (int) x.size();
		int i,j;
		Matrix<T> xx(n-p,p);
		vector<T> y(n-p);
		vector<T> rx(p+1, 0), new_a(p,0);
		vector<T> f, g;
		T d, var;
		r.clear();
		r.resize(p,0);

		switch (methods) {
		case AR_SELF_COR:
			/*
			 * use Levinson_Durbin algorithm, see xian_dai_shu_zhi_xin_hao_chu_li P139
			 */
			for (i=0; i<=p; i++) {
				for (j=0; j<n-i; j++)
					rx[i] += x[j] * x[i+j];
				//chenyu add	rx[i] = rx[i] / (n-i);
			}
			a=new_a;
			var = rx[0];
			//cout << "var0="<<var<<"\n";
			for (i=0; i<p; i++) {
				d = rx[i+1];
				for (j=i; j>0; j--)
					d += a[i-j] * rx[j];
				r[i] = d / var;
				var = (1-r[i]*r[i]) * var;
				for (j=0; j<i; j++)
					new_a[j] =a[j] - r[i]*a[i-j-1];
				new_a[i] = -r[i];
				a=new_a;
				/*cout<<'r'<<i<<'='<<r[i]<<",var"<<i<<'='<<var/n<<'\n';
				for (j=0; j<=i; j++)
					cout<<'a'<<j<<'='<<a[j]<<',';
				cout<<'\n';*/
			}
			return var/n;
			break;
		case AR_FORWARD_MMSE:
			/*
			 * use correlation algorithm, see xian_dai_shu_zhi_xin_hao_chu_li P140
			 */
			for (i=0; i<n-p; i++) {
				y[i] = x[p+i];
				for (j=0; j<p; j++)
					xx(i,j) = -x[p+i-j-1];
			}
			return mmse_solve(xx, y, a)/(n-p);
			break;
		case AR_BACKWARD_MMSE:
			for (i=0; i<n-p; i++) {
				y[i] = x[i];
				for (j=0; j<p; j++)
					xx(i,j) = -x[i+j+1];
			}
			return mmse_solve(xx, y, a)/(n-p);
			break;
		case AR_BURG:
			/*
			 * use burg algorithm, see xian_dai_shu_zhi_xin_hao_chu_li P142
			 */
			T s;
			f = x;
			g = x;
			var = 0;
			for (i=0; i<n; i++)
				var += x[i] * x[i];
			new_a.assign(p+1, 0);
			new_a[0] = 1;
			a = new_a;
			for (i=0; i<p; i++) {
				d = 0;
				s = 0;
				for (j=i; j<n-1; j++) {
					d += f[j+1] *  g[j];
					s += f[j+1] * f[j+1] + g[j] * g[j];
				}
				r[i] = 2*d /s;
				for (j=n-1; j>i; j--) { //warning
					g[j] = g[j-1]-r[i] * f[j];
					f[j] = f[j] - r[i] *g[j-1];
				}
				for (j=1; j<=i+1; j++)
					new_a[j] =a[j] - r[i]*a[i+1-j];
				new_a[0] = 1;
				a = new_a;
				var = (1-r[i]*r[i]) * var;
				cout<<'r'<<i<<'='<<r[i]<<",var"<<i<<'='<<var/n<<'\n';
				for (j=1; j<=i+1; j++)
					cout<<'a'<<j<<'='<<a[j]<<',';
				cout<<'\n';
			}
			a.erase(a.begin());
			return var/n;
			break;
		default:
			return -1;
		}
	}

	/*
	 * y=x*a, give y, x, solve a
	 * (xT*y) = (xT*x) a
	 * return (y-x*a)T * (y-x*a)
	 */
	static T mmse_solve(const Matrix<T> &x, const vector<T> &y, vector<T> &a) {
		if (x.RowNo() < x.ColNo())
			throw MyException("X row is less than x column in mmse_solve");

		Matrix <T> x_x, y_x, yy(y), e, am;

		x_x = (~x) * x;
		y_x = (~x) * yy;
		am = x_x.Solve(y_x);
		a =am.Col(0);
		e = yy - x*am;
		return (~e * e)(0, 0);
	}

	//return cos(n+1)w/2 + b[0]cos((n-1)w/2) + b[1]cos((n-3)w/2)... +b[n/2-1]cos(w/2)
	static T poly_cos(const vector<T> &b, T w)
	{
		int n = b.size() *2;
		T s = cos((n+1)*w/2);
		for (int i=0; i<n/2; i++)
			s += b[i] * cos ((n-1-2*i)*w/2);
		return s;
	}
	//return sin(n+1)w/2 + b[0]sin((n-1)w/2) + b[1]sin((n-3)w/2)... +b[n/2-1]sin(w/2)
	static T poly_sin(const vector<T> &b, T w)
	{
		int n = b.size() *2;
		T s = sin((n+1)*w/2);
		for (int i=0; i<n/2; i++)
			s += b[i] * sin ((n-1-2*i)*w/2);
		return s;
	}
	/*
	 * A(z) = 1 + a[0]z^-1 + a[1]z^-2 + a[2]z^-3......+a[n-1]z^-n
	 * P(z) = 1 + (a[0]+a[n-1])z^-1 + (a[1]+a[n-2])z^-2....+(a[n-1]+a[0])z^-n + z^-(n+1) (z=-1 is one root)
	 *      = z^-(n+1)/2 { (z^(n+1)/2+z^-(n+1)/2) + (a[0]+a[n-1]) (z^(n-1)/2+z^-(n-1)/2)...}
	 *      = z^-(n+1)/2 {cos(n+1)w/2 + (a[0]+a[n-1])cos((n-1)w/2) + (a[1]+a[n-2])cos((n-3)w/2)...(a[n/2-1]+a[n/2])cos(w/2)} (n%2==0)
	 *
	 * Q(z) = 1 + (a[0]-a[n-1])z^-1 + (a[1]-a[n-2])z^-2....+(a[n-1]-a[0])z^-n - z^-(n+1) (z=1 is one root)
	 *      = z^-(n+1)/2 { (z^(n+1)/2-z^-(n+1)/2) + (a[0]-a[n-1]) (z^(n-1)/2-z^-(n-1)/2)...}
	 *      = z^-(n+1)/2 {sin(n+1)w/2 + (a[0]-a[n-1])sin((n-1)w/2) + (a[1]-a[n-2])sin((n-3)w/2)...(a[n/2-1]-a[n/2])sin(w/2)} (n%2==0)
	 * Return 0 if success
	 */
	static int az_lsp(const vector<T> &az, vector<T> &lsp, int M=60, T precision=0.0001)
	{
		ASSERT(az.size()%2==0);
		unsigned n=az.size();
		vector<T> p(n/2), q(n/2);
		vector<T> nlsp;
		T w0, w1, w, f0, f1, fw, wi, fi;

		nlsp.clear();
		for (unsigned i=0; i<p.size(); i++) {
			p[i] = az[i] + az[az.size()-1-i];
			q[i] = az[i] - az[az.size()-1-i];
		}

		w0 = 0;
		f0 = poly_cos(p, w0);
		for (int i=1; i<M; i++, w0=wi, f0=fi) {
			wi = PI *i/M;
			fi = poly_cos(p,wi);
			if (f0*fi<0) {
				w = w1 = wi;
				fw = f1 = fi;
				while (fabs(fw)>precision) {
					w = (w0+w1)/2;
					fw = poly_cos(p,w);
					if (fw*f0<0) {
						w1 = w;
						f1 = fw;
					} else {
						w0 = w;
						f0 = fw;
					}
				}
				nlsp.push_back(w);
				if (nlsp.size()==n/2)
					break;
			}
		}
		if (nlsp.size()!=n/2) {
			printf("\nerror lsp p root!=n/2:");
			for (unsigned j=0; j<az.size(); j++)
				printf("%f,",az[j]);
			return -1;
		}


		w0 = PI/M;
		f0 = poly_sin(q, w0);
		for (int i=2; i<=M; i++, w0=wi, f0=fi) {
			wi = PI *i/M;
			fi = poly_sin(q,wi);
			if (f0*fi<0) {
				w = w1 = wi;
				fw = f1 = fi;
				while (fabs(fw)>precision) {
					w = (w0+w1)/2;
					fw = poly_sin(q,w);
					if (fw*f0<0) {
						w1 = w;
						f1 = fw;
					} else {
						w0 = w;
						f0 = fw;
					}
				}
				nlsp.push_back(w);
				if (nlsp.size()==n)
					break;
			}
		}
		if (nlsp.size()!=n){
			cout<<"\nerror lsp q root!=n:";
			for (unsigned j=0; j<az.size(); j++)
				printf("%f,",az[j]);
			return -2;
		}
		lsp = nlsp;
		return 0;
	}

	/*
	 * P(z) = 1 + (a[0]+a[n-1])z^-1 + (a[1]+a[n-2])z^-2....+(a[n-1]+a[0])z^-n + z^-(n+1)
	 *      = (1+z^-1) {(1-2cos(lsp[0])z^-1+z^-2) * ....(1-2cos(lsp[n/2-1])z^-1+z^-2)} (n%2==0)
	 *
	 * Q(z) = 1 + (a[0]-a[n-1])z^-1 + (a[1]-a[n-2])z^-2....+(a[n-1]-a[0])z^-n - z^-(n+1) (z=1 is one root)
	 *      = (1-z^-1) {(1-2cos(lsp[n/2])z^-1+z^-2) * ....(1-2cos(lsp[n-1])z^-1+z^-2)} (n%2==0)
	 *
	 * A(z) = (P(z)+Q(z))/2
	 */
	static void lsp_az(const vector<T> &lsp, vector<T> &az)
	{
		Polynomial<T> p, q;
		vector<T> v;
		unsigned n = lsp.size();
		v.push_back(1.0);
		v.push_back(1.0);
		p = v;
		v[1]=-1.0;
		q = v;
		v.resize(3);
		for (unsigned i=0; i<n/2; i++) {
			v[0] =1;
			v[2] =1;
			v[1] =-cos(lsp[i])*2;
			p*=v;
		}
		for (unsigned i=n/2; i<n; i++) {
			v[0] =1;
			v[2] =1;
			v[1] =-cos(lsp[i])*2;
			q*=v;
		}
		p+=q;
		az.resize(n);
		for (unsigned i=1; i<=n; i++)
			az[i-1] =p[i]/2;
	}

};

template <class T>
class AMRProcess : public LinearEstimate <T>{
public:
	vector<T> buf;
	int n, m;
	double cos_freq;
	vector<T> w1;
	vector<T> w2;
	vector<T> A[4], a[4];
	vector<T> LSP[4], lsp[4];
	Filter<T> *f;
	bool tx;
public:
	AMRProcess(bool _tx=true)
	{
		tx = _tx;
		n =240;
		m =10;
		cos_freq = 1.7;
		buf.resize(80, 0);
		w1.resize(240, 0);
		w2.resize(240, 0);
		for (int i=0; i<n; i++) {
			w1[i] = (i<160) ? 0.54 -0.46*cos(PI*i/159) : 0.54+0.46*cos(PI*(n-160)/79);
			w2[i] = (i<232) ? 0.54 -0.46*cos(PI*2*i/(2*232-1)) : cos(PI*2*(i-232)/31);
		}
		for (int i=0; i<4; i++) {
			if (i==0) {
				A[i].resize(m+1, 1.0);
				a[i].resize(m+1, 1.0);
			} else {
				A[i].resize(m, 1.0);
				a[i].resize(m, 1.0);
			}
			LSP[i].resize(m, 0);
			lsp[i].resize(m, 0);
		}
		f = NULL;
	}

	~AMRProcess()
	{
		delete f;
	}

	void change_lsp(const vector<double> &L, vector<double> &l)
	{
		ASSERT(L.size()==l.size());
		if (0) {
			if (tx) {
				for (unsigned i=0; i<l.size()/2; i++) {
					l[i] = L[l.size()/2+i] -L[0];
					l[l.size()/2+i] = L[i+1] -L[0];
				}
				l[l.size()-1] = L[l.size()-1];
			} else {
				l[0] = L[l.size()-1] - L[l.size()/2-1];
				for (unsigned i=0; i<l.size()/2; i++) {
					l[l.size()/2+i] = L[i] +l[0];
					if (i!=0)
						l[i] = L[l.size()/2+i-1] +l[0];
				}
			}
		} else {
			l =L;
			reverse(l.begin(), l.end());
			for (unsigned j=0; j<l.size(); j++)
				l[j] = PI-l[j];
		}

	}

	/*
	 * in  => A1, A3, => LSP1, LSP3 => LSP0, LSP2 => A0, A2 => d --> out
	 *                     ||                                     /
	 *                   lsp1, lsp3 => lsp0, lsp2 => a0,a1,a2,a3 /
	 */
	void transform(const vector<int> & in, vector<int> & out)
	{
		vector<T> yy, xx, x1(n), x2(n);
		vector<T> old_LSP3, old_lsp3;

		if (f==NULL)
			f = new Filter<T>(Filter<T>::IIR, A[0], a[0]);

		ASSERT(in.size()==160);
		CCfits::fill(in, xx);
		xx.insert(xx.begin(), buf.begin(), buf.end());
		buf.assign(xx.begin()+160, xx.end());
		ASSERT(buf.size()==80 && xx.size()==(unsigned) n);
		old_LSP3 = LSP[3];
		if (FFT::energy(xx)>n*100) {
			for (int i=0; i<n; i++) {
				x1[i] = xx[i]*w1[i];
				x2[i] = xx[i]*w2[i];
			}
			ar_est(LinearEstimate<T>::AR_SELF_COR, x1, A[1], m);
#if USE_LSP
			if (az_lsp(A[1], LSP[1], 100)!=0)
				LSP[1] =old_LSP3;
#else
			LSP[1] = this->r;
#endif
			ar_est(LinearEstimate<T>::AR_SELF_COR, x2, A[3], m);
#if USE_LSP
			if (az_lsp(A[3], LSP[3], 100)!=0)
				LSP[3] =LSP[1];
#else
			LSP[3] = this->r;
#endif
		} else {
			LSP[1] =old_LSP3;
			LSP[3] =old_LSP3;
#if USE_LSP
			lsp_az(LSP[1], A[1]);
			lsp_az(LSP[3], A[3]);
#else
			r_az(LSP[1], A[1]);
			r_az(LSP[3], A[3]);
#endif
		}
		ASSERT(LSP[1].size()==(unsigned) m && LSP[3].size()==(unsigned) m);
		for (int i=0; i<m; i++) {
			LSP[0][i] = (LSP[1][i] + old_LSP3[i])/2;
			LSP[2][i] = (LSP[1][i] + LSP[3][i])/2;
		}
#if USE_LSP
		lsp_az(LSP[0], A[0]);
		lsp_az(LSP[2], A[2]);
#else
		r_az(LSP[0], A[0]);
		r_az(LSP[2], A[2]);
#endif
		old_lsp3 = lsp[3];
		change_lsp(LSP[1], lsp[1]);
		change_lsp(LSP[3], lsp[3]);
		for (int i=0; i<m; i++) {
			lsp[0][i] = (lsp[1][i] + old_lsp3[i])/2;
			lsp[2][i] = (lsp[1][i] + lsp[3][i])/2;
		}
		for (int i=0; i<4; i++) {
#if USE_LSP
			lsp_az(lsp[i], a[i]);
#else
			r_az(lsp[i], a[i]);
#endif
			ASSERT(A[i].size()==(unsigned) m && a[i].size()==(unsigned) m);
		}

		for (int i=0; i<4; i++) {
			vector <T> x, y;
			x.assign(xx.begin() +80+i*40, xx.begin()+i*40+120);
			A[i].insert(A[i].begin(), 1);
			a[i].insert(a[i].begin(), 1);
			f->dynamic_change(A[i], a[i]);
			f->process(x, y);
			yy.insert(yy.end(), y.begin(), y.end());
		}
		if (tx)
			for (unsigned i=0; i<yy.size(); i++)
				yy[i] = yy[i] *0.5;
		else
			for (unsigned i=0; i<yy.size(); i++)
				yy[i] = yy[i] * 2;
		ASSERT(yy.size() ==160);
		CCfits::fill(yy, out);
	}

	void transform2(const vector<int> & in, vector<int> & out)
	{
		vector<T> yy, xx;
		ASSERT(in.size()==160);
		vector<T> a(3), b(3);

		CCfits::fill(in, xx);
		cos_freq -= 0.05;
		if (cos_freq <1.26)
			cos_freq = 1.7;
		a[0] = 1;
		a[1] = cos_freq;
		a[2] = 0.95;
		b[0] = 1;
		b[1] = cos_freq*2/3;
		b[2] = 0.95;
		if (f==NULL)
			f = new Filter<T>(Filter<T>::IIR, b, a);
		if (tx)
			f->dynamic_change(b, a);
		else
			f->dynamic_change(a, b);
		f->process(xx, yy);
		if (tx)
		for (unsigned i=0; i<yy.size(); i++)
			yy[i] = yy[i] *0.5;
		ASSERT(yy.size() ==160);
		CCfits::fill(yy, out);
	}

	void transform3(const vector<int> & in, vector<int> & out)
	{
		vector<T> x1(n), xx;
		out.resize(in.size());
		for (unsigned i=0; i<in.size(); i++)
			out[i] = (i%2==0) ? in[i] : -in[i];

		CCfits::fill(in, xx);
		xx.insert(xx.begin(), buf.begin(), buf.end());
		buf.assign(xx.begin()+160, xx.end());
		ASSERT(buf.size()==80 && xx.size()==(unsigned) n);

		if (FFT::energy(xx)>n*100) {
			for (int i=0; i<n; i++)
				x1[i] = xx[i]*w1[i];
			ar_est(LinearEstimate<T>::AR_SELF_COR, x1, A[1], m);
			A[2] = A[1];
			A[3] = A[1];
			A[0] = A[1];
		}

		CCfits::fill(out, xx);
		xx.insert(xx.begin(), buf.begin(), buf.end());
		buf.assign(xx.begin()+160, xx.end());
		ASSERT(buf.size()==80 && xx.size()==(unsigned) n);

		if (FFT::energy(xx)>n*100) {
			for (int i=0; i<n; i++)
				x1[i] = xx[i]*w1[i];
			ar_est(LinearEstimate<T>::AR_SELF_COR, x1, a[1], m);
			a[2] = a[1];
			a[3] = a[1];
			a[0] = a[1];
		}
	}

};

template <class T>
class FreqSwap {
protected:
	Filter<T> * high;  //1800 ~ 1900 drop band
	Filter<T> * low; //1700 ~1800 drop band
	Delay<T> * delay1, *delay2;
	Mixer<T> * mix1, *mix2;
public:
	FreqSwap(vector <T> & high_coef, vector <T> & low_coef)
	{
		high = new Filter<T>(Filter<T>::FIR, high_coef);
		low = new Filter<T>(Filter<T>::FIR, low_coef);
		delay1 = new Delay<T> ((high_coef.size() -1)/2);
		delay2 = new Delay<T> ((low_coef.size()-2)/2);
		mix1 = new Mixer<T> (1, -1);
		mix2 = new Mixer<T> (2, 1);
	}

	~FreqSwap()
	{
		delete high;
		delete low;
		delete delay1;
		delete delay2;
		delete mix1;
		delete mix2;
	}

	void transform(const vector<int> & in, vector<int> & out)
	{
		vector<T> x1, x2, x, xd, x1_down(in.size()/2), x1_up(in.size()), x2d, x1_up_l, y;
		CCfits::fill(in, x);
		high->process(x, x1); //x1(z) = x(z) * f1(z)
		delay1->process(x, xd); //xd(z) = x(z) * z^-n
		mix1->mix(xd, x1, x2); //x2(z) = x(z) * (z^-n- f1(z))
		for (unsigned i=0; i<x1_down.size(); i++)
			x1_down[i] = (i%2==0) ? x1[i*2] : -x1[i*2];  //x1_down(w) = 0.5 x1((w+pi)/2)

		for (unsigned i=0; i<x1_up.size(); i++)
			x1_up[i] = (i%2==0) ? x1_down[i/2] : 0; //x1_up(w) = 0.5 x1(w+pi/2)
		low->process(x1_up, x1_up_l);
		delay2->process(x2, x2d);
		mix2->mix(x1_up_l, x2d, y);
		ASSERT(y.size()==in.size());
		CCfits::fill(y, out);
	}

	void transform2(const vector<int> & in, vector<int> & out)
	{		
		out.resize(in.size());
		for (unsigned i = 0; i<in.size(); i++)
			out[i] = (i % 2 == 0) ? in[i] : -in[i];
	}
};
#endif /* LINEARESTIMATE_H_ */
