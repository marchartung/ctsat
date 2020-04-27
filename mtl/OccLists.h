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
#ifndef SOURCES_MTL_OCCLISTS_H_
#define SOURCES_MTL_OCCLISTS_H_

#include "mtl/Vec.h"

namespace ctsat {

//=================================================================================================
// OccLists -- a class for maintaining occurence lists with lazy deletion:

template<class Idx, class Vec, class Deleted>
class OccLists {
	vec<Vec> occs;
	vec<char> dirty;
	vec<Idx> dirties;
	Deleted deleted;

public:
	OccLists(const Deleted& d) :
			deleted(d) {
	}

	void init(const Idx& idx) {
		occs.growTo(toInt(idx) + 1);
		dirty.growTo(toInt(idx) + 1, 0);
	}
	// Vec&  operator[](const Idx& idx){ return occs[toInt(idx)]; }
   Vec const & operator[](const Idx& idx) const {
      return occs[toInt(idx)];
   }
	Vec& operator[](const Idx& idx) {
		return occs[toInt(idx)];
	}
	Vec& lookup(const Idx& idx) {
		if (dirty[toInt(idx)])
			clean(idx);
		return occs[toInt(idx)];
	}

	void cleanAll();
	void clean(const Idx& idx);
	void smudge(const Idx& idx) {
		if (dirty[toInt(idx)] == 0) {
			dirty[toInt(idx)] = 1;
			dirties.push(idx);
		}
	}

	void clear(bool free = true) {
		occs.clear(free);
		dirty.clear(free);
		dirties.clear(free);
	}
};

template<class Idx, class Vec, class Deleted>
void OccLists<Idx, Vec, Deleted>::cleanAll() {
	for (int i = 0; i < dirties.size(); i++)
		// Dirties may contain duplicates so check here if a variable is already cleaned:
		if (dirty[toInt(dirties[i])])
			clean(dirties[i]);
	dirties.clear();
}

template<class Idx, class Vec, class Deleted>
void OccLists<Idx, Vec, Deleted>::clean(const Idx& idx) {
	Vec& vec = occs[toInt(idx)];
	int i, j;
	for (i = j = 0; i < vec.size(); i++)
		if (!deleted(vec[i]))
			vec[j++] = vec[i];
	vec.shrink(i - j);
	dirty[toInt(idx)] = 0;
}

}

#endif /* SOURCES_MTL_OCCLISTS_H_ */
