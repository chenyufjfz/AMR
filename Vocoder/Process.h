/*
 * Process.h
 *
 *  Created on: Mar 10, 2013
 *      Author: yuchen
 */

#ifndef PROCESS_H_
#define PROCESS_H_

#include <vector>
using namespace std;

template <class T>
class Process {
public:
	virtual void process(const vector<T> &x, vector<T> &y)=0;
	virtual ~Process() {}
};

template <class T>
class ProcessFlow : public Process<T> {
private:
	vector<Process<T> *> pset;
	vector<bool> del;
public:
	~ProcessFlow() {
		for (unsigned i=0; i<del.size(); i++)
			if (del[i])
				delete pset[i];
	}
	void append(Process<T> * p, bool d=false) {
		pset.push_back(p);
		del.push_back(d);
	}
	void process(const vector<T> &x, vector<T> &y) {
		vector<T> m1,m2;
		unsigned i;
		if (pset.size()==1) {
			pset[0]->process(x, y);
			return;
		}
		pset[0]->process(x,m1);
		for (i=0; i+2<pset.size(); i++) {
			if (i%2==0)
				pset[i+1]->process(m1,m2);
			else
				pset[i+1]->process(m2,m1);
		}
		if (i%2==0)
			pset[i+1]->process(m1,y);
		else
			pset[i+1]->process(m2,y);
	}
};
#endif /* PROCESS_H_ */
