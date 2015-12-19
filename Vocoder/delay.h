/*
 * delay.h
 *
 *  Created on: Mar 16, 2013
 *      Author: yuchen
 */

#ifndef DELAY_H_
#define DELAY_H_

#include "Process.h"

template <class T>
class Delay : public Process<T> {
private:
	vector<T> buf_in;
	int delay;
public:
	Delay(int d) {
		delay =d;
	}
	/*
		 * prefill zero for x and y
		 * xx:  0, 0,...0, x[0], x[1], x[2],...., x[n-1]
		 *      |<-delay->| |<-       n               ->|
		 */
	void process(const vector<T> &x, vector<T> &y) {
		vector <T> xx;
		if (buf_in.empty())
			buf_in.assign(delay, 0);
		xx = buf_in;
		xx.insert(xx.end(), x.begin(), x.end());
		y.assign(xx.begin(), xx.begin()+x.size());
		buf_in.assign(xx.begin()+x.size(), xx.end());
		ASSERT(buf_in.size()==(unsigned) delay);
	}
};


#endif /* DELAY_H_ */
