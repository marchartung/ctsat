/*****************************************************************************************[Main.cc]
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

#ifndef SOURCES_PARALLEL_PARALLELSTATISTIC_H_
#define SOURCES_PARALLEL_PARALLELSTATISTIC_H_

#include "core/Statistic.h"

namespace CTSat
{
struct PStatistic
{
   uint64_t nAdded;
   uint64_t starts;
   uint64_t conflicts;
   uint64_t decisions;
   uint64_t propagations;
   uint64_t spropagations;
   uint64_t watchedLearnts;

   uint64_t nPromoted;
   uint64_t nReceivedClauses;
   uint64_t nSendClauses;
   uint64_t nHoldBackImported;

   PStatistic()
         : nAdded(0),
           starts(0),
           conflicts(0),
           decisions(0),
           propagations(0),
           spropagations(0),
           watchedLearnts(0),

           nPromoted(0),
           nReceivedClauses(0),
           nSendClauses(0),
           nHoldBackImported(0)
   {
   }
   void reset()
   {
      *this = PStatistic();
   }

   void add(PStatistic const & stat)
   {
      nAdded += stat.nAdded;
      starts += stat.starts;
      conflicts += stat.conflicts;
      decisions += stat.decisions;
      propagations += stat.propagations;
      spropagations += stat.spropagations;
      watchedLearnts += stat.watchedLearnts;

      nPromoted += stat.nPromoted;
      nReceivedClauses += stat.nReceivedClauses;
      nSendClauses += stat.nSendClauses;
      nHoldBackImported += stat.nHoldBackImported;
   }

   void add(Statistic const & stat)
   {
      ++nAdded;
      starts += stat.starts;
      conflicts += stat.conflicts;
      decisions += stat.decisions;
      propagations += stat.propagations;
      spropagations += stat.s_propagations;
      watchedLearnts += stat.nWatchedLearnts;

      nPromoted += stat.nPromoted;
      nReceivedClauses += stat.nReceivedClauses;
      nSendClauses += stat.nSendClauses;
      nHoldBackImported += stat.nHoldBackImported;
   }

   void print()
   {
      uint64_t const added = (nAdded > 0) ? nAdded : 1;

      printf("c ############################## per Solver ############################\n");
      printf("c rest:%-12" PRIu64" confl:%-12" PRIu64" dec:%-12" PRIu64" prop:%-12" PRIu64"\n",
             starts / added, conflicts / added, decisions / added, propagations / added);
      printf("c learnts:%-12" PRIu64 " vivi load:%2.2f%%\n", watchedLearnts / added,
             static_cast<double>(spropagations * 100) / (spropagations + std::max(propagations,1ul)));

      printf("c prom:%-12" PRIu64" receive:%-12" PRIu64" send:%-12" PRIu64"holdIm:%-12" PRIu64"\n",
             nPromoted / added, nReceivedClauses / added, nSendClauses / added,
             nHoldBackImported / added);
   }
};

}

#endif /* SOURCES_PARALLEL_PARALLELSTATISTIC_H_ */
