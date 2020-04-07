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

#ifndef SOURCES_MPI_MPISTATISTC_H_
#define SOURCES_MPI_MPISTATISTC_H_

#include "database/BasicTypes.h"
#include "parallel/ParallelStatistic.h"
#include "utils/Timer.h"

namespace CTSat
{
struct MPIStatistc : public PStatistic
{

   MPIStatistc(int const rank)
         : abort(false),
           res(lbool::Undef()),
           rank(rank)
   {

   }
   void add(MPIStatistc const & h)
   {
      this->abort |= h.abort;
      if (!h.res.isUndef())
      {
         assert(res == h.res || res.isUndef());
         rank = (res.isUndef() || rank > h.rank) ? h.rank : rank;
         res = h.res;
      }
      PStatistic::add(h);
   }

   bool abort;
   lbool res;
   int rank;

};

struct MPICommunicationStat
{
   int const nRanks;
   uint64_t const bufferSize;
   uint64_t numRounds;
   uint64_t numBytesReceived;
   uint64_t numBytesSend;
   double sumFillState;
   double maxFillState;
   uint64_t clausesSend;
   uint64_t clausesRecv;
   uint64_t sumLbdSend;
   uint64_t sumLbdRecv;
   Timer t;

   MPICommunicationStat(int const nRanks, uint64_t const bufferSize)
         : nRanks(nRanks),
           bufferSize(bufferSize),
           numRounds(0),
           numBytesReceived(0),
           numBytesSend(0),
           sumFillState(0.0),
           maxFillState(0.0),
           clausesSend(0),
           clausesRecv(0),
           sumLbdSend(0),
           sumLbdRecv(0)
   {
   }

   void addComm(
                bool const isSend,
                uint64_t const nBytes,
                uint64_t const nClauses,
                uint64_t const sLbd)
   {
      if (isSend)
      {
         ++numRounds;
         numBytesSend += nBytes;
         clausesSend += nClauses;
         sumLbdSend += sLbd;

      } else
      {
         numBytesReceived += nBytes;
         clausesRecv += nClauses;
         sumLbdRecv += sLbd;
      }
      double const fillState = static_cast<double>(nBytes + sizeof(MPIStatistc)) / bufferSize;
      sumFillState += fillState / nRanks;
      maxFillState = (maxFillState < fillState) ? fillState : maxFillState;
   }

   void print() const
   {

      printf("c ##############################  per Node  ############################\n");
      printf("c rounds: %-6" PRIu64 " cl send: %-12" PRIu64 " recv: %-12" PRIu64 "\n", numRounds,
             clausesSend, clausesRecv);

      printf("c buffer usage: %4.2f%% max: %4.2f%%\n",
             100.0 * sumFillState / (numRounds * nRanks), 100.0 * maxFillState);
      printf("c bandwidth: %4.2f MB/s interval: %2.2f rounds/s\n",static_cast<double>(numBytesSend + numBytesReceived)
             / (1024.0 * 1024.0 * t.getPassedTime()),static_cast<double>(numRounds) / t.getPassedTime());
      printf("c ######################################################################\n\n");
   }
};

}

#endif /* SOURCES_MPI_MPISTATISTC_H_ */
