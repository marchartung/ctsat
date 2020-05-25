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
#ifndef DATABASE_RESOLUTIONTRACKER_H_
#define DATABASE_RESOLUTIONTRACKER_H_

#include "mtl/Vec.h"
#include "mtl/Alg.h"

#include <iostream>

namespace ctsat {
template<bool active, typename Lit, typename Clause>
class ResolutionTracker {
public:
	ResolutionTracker();

	void addResolvent(Clause const & c);

	template<typename LitVec, typename Level>
	void checkResolution(LitVec const & c, Level const & level);

	void startTracking();

private:
	vec<Lit> reso;
	vec<Clause const *> resolvents;

	bool isResolved(Lit const & l) const {
		assert(active);
		bool foundNeg = false, foundPos = false;
		for (int i = 0; i < resolvents.size(); ++i) {
			Clause const & r = *resolvents[i];
			for (int j = 0; j < r.size(); ++j) {
				foundPos |= r[j] == l;
				foundNeg |= r[j] == ~l;
				if (foundPos && foundNeg)
					return true;
			}
		}
		return false;
	}

	template<typename LitVec>
	void printClause(LitVec const & c) {
		assert(active);
		std::cout << "[";
		for (int i = 0; i < c.size(); ++i)
			std::cout << ((c[i].sign()) ? "-" : "") << c[i].var() << " ";
		std::cout << "]" << std::endl;
	}

};

template<bool active, typename Lit, typename Clause>
inline ResolutionTracker<active, Lit, Clause>::ResolutionTracker() {
}

template<bool active, typename Lit, typename Clause>
inline void ResolutionTracker<active, Lit, Clause>::addResolvent(
		const Clause& c) {
	if (active)
		resolvents.push(&c);
}

template<bool active, typename Lit, typename Clause>
template<typename LitVec, typename Level>
inline void ResolutionTracker<active, Lit, Clause>::checkResolution(
		const LitVec& c, Level const & level) {
	if (active) {
		reso.clear();
		bool res = true;
		for (int i = 0; i < resolvents.size() && res; ++i) {
			Clause const & r = *resolvents[i];
			for (int j = 0; j < r.size(); ++j) {
				if (level(r[j]) == 0)
					continue;
				if (isResolved(r[j]) == find(c, r[j])) {
					std::cout << r[j].var() << " is resolved: "
							<< isResolved(r[j]) << " is in resolution: "
							<< find(c, r[j]) << std::endl;
					res = false;
					break;
				}
			}
		}
		if (!res) {
			std::cout << "Wrong resolution:\n";
			for (int i = 0; i < resolvents.size(); ++i) {
				std::cout << i << ":";
				printClause(*resolvents[i]);
			}

			std::cout << "res:";
			printClause(c);
			assert(false);
		}
	}
}

template<bool active, typename Lit, typename Clause>
inline void ResolutionTracker<active, Lit, Clause>::startTracking() {
	if (active)
		resolvents.clear();
}
}

#endif /* DATABASE_RESOLUTIONTRACKER_H_ */
