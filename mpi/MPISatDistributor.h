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

#ifndef MPI_MPISATDISTRIBUTOR_H_
#define MPI_MPISATDISTRIBUTOR_H_

#include <tuple>
#include <algorithm>
#include <mpi.h>
#include <cstring>
#include "mpi/MPIHelper.h"
#include "database/BasicTypes.h"
#include "mtl/Vec.h"

namespace ctsat
{

// HELPER:


inline int32_t calcResultInt(int const rank, lbool const res)
{
   assert(rank >= 0);
   if (res.isUndef())
      return std::numeric_limits<int32_t>::min();
   else if (res.isFalse())
      return -(rank + 1);
   else
      return rank + 1;
}

inline int calcRankFromResult(int32_t const res)
{
   if (res == std::numeric_limits<int32_t>::min())
      return -1;
   else
      return std::abs(res) - 1;
}

inline constexpr lbool calcLboolFromResult(int32_t const res)
{
   return
         (res == std::numeric_limits<int32_t>::min()) ?
               lbool::Undef() : ((res < 0) ? lbool::False() : lbool::True());
}

class MPISatDistributor
{
 public:

   MPISatDistributor();
   ~MPISatDistributor();

   void sendInstance(void const * data, uint64_t const nBytes);
   std::tuple<void const *, uint64_t> receiveInstance();
   void freeInstance();

   void receiveModel(vec<lbool> & model);
   void sendModel(vec<lbool> const & model);

   bool trySetRootResult(lbool const res);
   lbool getRootResult();
   int getWinner();

 private:
   int rank;
   int32_t * resultWinMem;
   uint8_t * instanceBuffer;
   uint8_t * modelBuffer;
   MPI_Request instanceReq;
   MPI_Win resultWin;

   int32_t pollResult();
};

inline MPISatDistributor::MPISatDistributor()
      : rank(-1),
        resultWinMem(NULL),
        instanceBuffer(NULL),
        modelBuffer(NULL),
        instanceReq(MPI_Request()),
        resultWin(MPI_Win())

{
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int nBytes = 0;
   if(rank == 0)
   {
      nBytes = sizeof(int32_t);
      resultWinMem = allocMpiMem<int32_t>(nBytes);
      *resultWinMem = calcResultInt(0, lbool::Undef());
      MPI_Win_create(resultWinMem,nBytes,1,MPI_INFO_NULL,MPI_COMM_WORLD,&resultWin);
//      MPI_Win_sync(resultWin);

   }
   else
      MPI_Win_create(resultWinMem,0,1,MPI_INFO_NULL,MPI_COMM_WORLD,&resultWin);
}

MPISatDistributor::~MPISatDistributor()
{
   MPI_Win_free(&resultWin);
   freeMpiMem(resultWinMem);
   freeMpiMem(instanceBuffer);
}

bool MPISatDistributor::trySetRootResult(lbool const result)
{
   assert(!result.isUndef());
   int32_t des = calcResultInt(rank, result), exp = calcResultInt(0, lbool::Undef()), res;
   MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, resultWin);
   MPI_Compare_and_swap(&des, &exp, &res, MPI_INT32_T, 0, 0, resultWin);
   MPI_Win_unlock(0, resultWin);

   if (res == exp)
      assert(pollResult() == des);

   return res == exp;
}

int32_t MPISatDistributor::pollResult()
{
   int32_t res = calcResultInt(0, lbool::Undef());
   while (true)
   {
      MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, resultWin);
      if (rank == 0)
         res = *resultWinMem;
      else
         MPI_Get(&res, 1, MPI_INT32_T, 0, 0, 1, MPI_INT32_T, resultWin);
      MPI_Win_unlock(0, resultWin);

      if (res != calcResultInt(0, lbool::Undef()))
         break;
      else
         usleep(100);
   }
   return res;
}

lbool MPISatDistributor::getRootResult()
{
   return calcLboolFromResult(pollResult());
}
int MPISatDistributor::getWinner()
{
   return calcRankFromResult(pollResult());
}

inline void MPISatDistributor::sendModel(vec<lbool> const & model)
{
   static_assert(sizeof(lbool) == 1, "lbool size has changed");
   assert(rank != 0);
   uint64_t const nBytes = model.size() * sizeof(lbool);
   modelBuffer = allocMpiMem<uint8_t>(nBytes);

   lbool * const m = reinterpret_cast<lbool*>(modelBuffer);
   for (int i = 0; i < model.size(); ++i)
   {
      m[i] = model[i];
   }
   memcpy(modelBuffer,model.getData(),model.size());
   MPI_Send(modelBuffer, nBytes, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
   freeMpiMem(modelBuffer);
}

inline void MPISatDistributor::receiveModel(vec<lbool> & model)
{
   int const winnerRank = getWinner();
   static_assert(sizeof(lbool) == 1, "lbool size has changed");
   assert(rank == 0 && rank != std::abs(winnerRank));

   uint64_t const nBytes = model.size() * sizeof(uint8_t);
   modelBuffer = allocMpiMem<uint8_t>(nBytes);
   MPI_Status status;
   MPI_Recv(modelBuffer, nBytes, MPI_BYTE, winnerRank, 0, MPI_COMM_WORLD, &status);
   int recvSz;
   MPI_Get_count(&status, MPI_BYTE, &recvSz);
   lbool const * m = reinterpret_cast<lbool const *>(modelBuffer);
   assert(model.size() >= recvSz);
   for (int i = 0; i < recvSz; ++i)
   {
      model[i] = m[i];
   }
   freeMpiMem(modelBuffer);
}
inline void MPISatDistributor::sendInstance(const void* data, const uint64_t nBytes)
{
   int waitFlag = 0;
// send size:
   uint64_t bytesSend = nBytes;
   MPI_Bcast(&bytesSend, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);  // local completition should be fast
   if (nBytes > 0)
   {
      // send instance
      instanceBuffer = allocMpiMem<uint8_t>(nBytes);
      memcpy(instanceBuffer, data, nBytes);
      MPI_Ibcast(instanceBuffer, nBytes, MPI_BYTE, 0, MPI_COMM_WORLD, &instanceReq);
      MPI_Test(&instanceReq, &waitFlag, MPI_STATUS_IGNORE);
   }
}

inline std::tuple<const void*, uint64_t> MPISatDistributor::receiveInstance()
{
// receive number of bytes to allocate:
   uint64_t nBytes = 0;
   MPI_Bcast(&nBytes, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

   if (nBytes > 0)
   {
      // receive instance:
      instanceBuffer = allocMpiMem<uint8_t>(nBytes);
      MPI_Ibcast(instanceBuffer, nBytes, MPI_BYTE, 0, MPI_COMM_WORLD, &instanceReq);
      MPI_Wait(&instanceReq, MPI_STATUS_IGNORE);
   }
   return std::tie(instanceBuffer, nBytes);
}

inline void MPISatDistributor::freeInstance()
{
   if (rank == 0)
   {
      MPI_Wait(&instanceReq, MPI_STATUS_IGNORE);
   }
   freeMpiMem(instanceBuffer);
}

}

#endif /* MPI_MPISATDISTRIBUTOR_H_ */
