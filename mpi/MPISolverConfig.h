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

#ifndef MPI_MPISOLVERCONFIG_H_
#define MPI_MPISOLVERCONFIG_H_

#include "utils/CPUBind.h"
#include "initial/Inputs.h"

namespace ctsat
{
class MPISolverConfig
{
public:
   MPISolverConfig(int const rank);

   int getMaxNumThreads() const;
   int getNumThreads() const;

private:
   int maxThreads;
   int usedThreads;
};


MPISolverConfig::MPISolverConfig(int const rank)
: maxThreads(Inputs::nThreads),
  usedThreads(Inputs::nThreads)
{
   if(Inputs::mpiAutoThreads)
   {
      // we map relative to the number of available threads
      // using this sequence against the rank: 0 means all threads, 1 means one thread, < 1 means maxThread/"<1"
      int rankMap[] = {0,2,0,2,0,2,0,2,0,1,0,2,0,2,0,2,0,2,0,1};

      NumaAwareSet const & numa = NumaAwareSet::instance;
      maxThreads = numa.getNumCores();
      int const rrank = rank%20;
      usedThreads = (rankMap[rrank] < 2) ? ((rankMap[rrank] == 1) ? numa.getNumNumaNodes() : maxThreads) : maxThreads/rankMap[rrank];
   }
   assert(usedThreads > 0);
}

int MPISolverConfig::getMaxNumThreads() const
{
   return maxThreads;
}
int MPISolverConfig::getNumThreads() const
{
   return usedThreads;
}
}




#endif /* MPI_MPISOLVERCONFIG_H_ */
