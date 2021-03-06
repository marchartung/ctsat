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

#ifndef SOURCES_CORE_STATISTIC_H_
#define SOURCES_CORE_STATISTIC_H_

#include <cstdint>
#include <cinttypes>
#include <cstdio>
#include "utils/System.h"

namespace ctsat
{
class Statistic
{
 public:

   Statistic()
         : restarts(0),
           conflicts(0),
           decisions(0),
           propagations(0),
           clauses_literals(0),
           nWatchedClauses(0),
           learnts_literals(0),
           nWatchedLearnts(0),
           nAdditionalLearnt(0),
           nPromoted(0),
           nReceivedClauses(0),
           nSendClauses(0),
           nLostClauses(0),
           nHoldBackImported(0),
           chrono_backtrack(0),
           non_chrono_backtrack(0),
           level_backtracked(0),
           s_propagations(0),
           simpDB_props(0),
           simpDB_assigns(0),
           global_lbd_sum(0)
   {
   }

   uint64_t restarts;

   uint64_t conflicts;

   uint64_t decisions;
   uint64_t propagations;

   uint64_t clauses_literals;
   uint64_t nWatchedClauses;
   uint64_t learnts_literals;
   uint64_t nWatchedLearnts;

   uint64_t nAdditionalLearnt;

   uint64_t nPromoted;
   uint64_t nReceivedClauses;
   uint64_t nSendClauses;
   uint64_t nLostClauses;
   uint64_t nHoldBackImported;

   uint64_t chrono_backtrack;
   uint64_t non_chrono_backtrack;
   uint64_t level_backtracked;

   uint64_t s_propagations;

   int64_t simpDB_props;  // Remaining number of propagations that must be made before next execution of 'simplify()'.
   int simpDB_assigns;  // Number of top-level assignments since last execution of 'simplify()'.

   float global_lbd_sum;

   void print() const
   {
      printf("c rest:%-12" PRIu64" confl:%-12" PRIu64" dec:%-12" PRIu64" prop:%-12" PRIu64"\n",
             restarts, conflicts, decisions, propagations);
      printf(
            "c learnts:%-12" PRIu64 " vivi load:%2.2f%%   mem: %5.2fMB\n",
            nWatchedLearnts,
            static_cast<double>(s_propagations * 100)
               / (s_propagations + std::max(propagations, 1ul)),
            memUsedPeak());
      printf("c ######################################################################\n");
   }

   void printFinal() const
   {
      double const cpu_time = cpuTime();
      printf("c restarts              : %-12" PRIu64"\n", restarts);
      printf("c conflicts             : %-12" PRIu64"   (%.0f /sec)\n", conflicts,
             conflicts / cpu_time);
      printf("c decisions             : %-12" PRIu64" (%.0f /sec)\n", decisions,
             decisions / cpu_time);
      printf("c propagations          : %-12" PRIu64"   (%.0f /sec)\n", propagations,
             propagations / cpu_time);
      printf("c level backtracked     : %-12" PRIu64"   (NCB %0.f%% , CB %0.f%%)\n",
             level_backtracked,
             (non_chrono_backtrack * 100) / (double) (non_chrono_backtrack + chrono_backtrack),
             (chrono_backtrack * 100) / (double) (non_chrono_backtrack + chrono_backtrack));
      printf("c additional learnt     : %-12" PRIu64"\n", nAdditionalLearnt);

      double const mem_used = memUsedPeak();
      if (mem_used != 0)
         printf("c Memory used           : %.2f MB\n", mem_used);
      printf("c CPU time              : %g s\n", cpu_time);
      if (nLostClauses > 0)
         std::cout
            << "c Warning: Lost "
            << nLostClauses
            << " clauses due to full buffer. Increase buffer sizes."
            << std::endl;
   }

};
}
#endif /* SOURCES_CORE_STATISTIC_H_ */
