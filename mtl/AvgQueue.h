/*****************************************************************************************
CTSat -- Copyright (c) 2020, Marc Hartung
                        Zuse Institute Berlin, Germany

Maple_LCM_Dist_Chrono -- Copyright (c) 2018, Vadim Ryvchin, Alexander Nadel

GlucoseNbSAT -- Copyright (c) 2016,Chu Min LI,Mao Luo and Fan Xiao
                           Huazhong University of science and technology, China
                           MIS, Univ. Picardie Jules Verne, France

MapleSAT -- Copyright (c) 2016, Jia Hui Liang, Vijay Ganesh

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef SOURCES_MTL_AVGQUEUE_H_
#define SOURCES_MTL_AVGQUEUE_H_

#include "mtl/Vec.h"
#include <cassert>

namespace ctsat {
template<typename T>
class AvgQueue {
	int max_sz, q_sz;
	int ptr;
	int64_t sum;
	vec<T> q;
public:
	AvgQueue(int sz) :
			max_sz(sz), q_sz(0), ptr(0), sum(0) {
		assert(sz > 0);
		q.growTo(sz);
	}
	inline bool full() const {
		return q_sz == max_sz;
	}
#ifdef INT_QUEUE_AVG
	inline T avg() const
	{
		assert(full());
		return sum / max_sz;
	}
#else
	inline double avg() const {
		assert(full());
		return sum / (double) max_sz;
	}
#endif
	inline void clear() {
		sum = 0;
		q_sz = 0;
		ptr = 0;
	}
	void push(T e) {
		if (q_sz < max_sz)
			q_sz++;
		else
			sum -= q[ptr];
		sum += e;
		q[ptr++] = e;
		if (ptr == max_sz)
			ptr = 0;
	}
};
}

#endif /* SOURCES_MTL_AVGQUEUE_H_ */
