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

#ifndef SOURCES_PARALLEL_MPIBROADCASTCONNECTOR_H_
#define SOURCES_PARALLEL_MPIBROADCASTCONNECTOR_H_

#include <cstdint>
#include <tuple>
#include <algorithm>
#include <mpi/mpi.h> // FIXME
#include <unistd.h>

#include "parallel/AtomicConnector.h"
#include "mpi/MPIStatistc.h"
#include "mpi/MPIHelper.h"

namespace ctsat
{

class MPIBroadcastConnector : public AtomicConnector
{
   MPIBroadcastConnector(MPIBroadcastConnector const &) = delete;
   MPIBroadcastConnector & operator=(MPIBroadcastConnector const &) = delete;
   MPIBroadcastConnector(MPIBroadcastConnector &&) = delete;
   MPIBroadcastConnector & operator=(MPIBroadcastConnector &&) = delete;

 public:

   struct size_type
   {
      AtomicConnector::size_type localPos;
      AtomicConnector::size_type globalPos;
      size_type(AtomicConnector::size_type const s)
            : localPos(s),
              globalPos(AtomicConnector::startPos())
      {
      }

      size_type()
            : size_type(AtomicConnector::startPos())
      {
      }

      size_type(size_type const& in)
            : localPos(in.localPos),
              globalPos(in.globalPos)
      {
      }

      size_type(size_type && in)
            : localPos(in.localPos),
              globalPos(in.globalPos)
      {
      }

      size_type & operator=(size_type const & in)
      {
         localPos = in.localPos;
         globalPos = in.globalPos;
         return *this;
      }

      size_type & operator=(size_type && in)
      {
         localPos = in.localPos;
         globalPos = in.globalPos;
         return *this;
      }
   };

   unsigned getUniqueId()
   {
      return NoConnector::getUniqueId() + rank * maxThreads;
   }

   MPIBroadcastConnector(
                         uint64_t bytesPerThread,
                         uint64_t bytesPerRank,
                         unsigned const nThreads,
                         unsigned const maxThreads,
                         MPI_Comm const comm = MPI_COMM_WORLD);
   ~MPIBroadcastConnector();

   uint64_t getSendBufferSize() const
   {
      return sendBufferSize - sizeof(MPIStatistc);
   }

   bool isAllowedToSend() const;

   void progress();

   MPIStatistc getAccumulatedHeader();

   void printState();

   bool isValid(size_type const & pos) const;
   size_type next(size_type const & pos);

   template <typename ClauseType>
   inline ClauseType const & get(size_type const & pos) const
   {
      assert(isValid(pos));
      bool const isLocal = AtomicConnector::ara.isValid(pos.localPos);
      AtomicRingAllocator::size_type rpos = (isLocal) ? pos.localPos : pos.globalPos;
      AtomicRingAllocator const & a = (isLocal) ? AtomicConnector::ara : mpiRecvClauses;
      return a.get<ClauseType>(rpos);
   }

   template <typename ClauseType, typename ... Args>
   inline bool exchange(uint64_t const & nBytes, Args ... args)
   {
      return AtomicConnector::exchange<ClauseType, Args...>(nBytes, args...);
   }

   inline uint64_t nFreeSendBytes() const
   {
      return sendBufferSize - sendBufferPos;
   }

   inline bool addToSend(void const * ptr, uint64_t const & nBytes)
   {
      assert(isAllowedToSend());
      if (nBytes + sendBufferPos > sendBufferSize)
         return false;
      memcpy(sendBuffer + sendBufferPos, ptr, nBytes);
      sendBufferPos += nBytes;
      return true;
   }

   uint8_t const * getRecvBuffer(int const rank) const
   {
      assert(isAllowedToSend());

      uint8_t const * res = recvBuffer + sizeof(MPIStatistc) + rank * sendBufferSize;
      return res;
   }

   bool shouldImport(size_type const & pos) const;

   AtomicRingAllocator & getMPIExchangeBuffer()
   {
      return mpiRecvClauses;
   }

   AtomicRingAllocator & getLocalExchangeBuffer()
   {
      return AtomicConnector::ara;
   }

   void send(const MPIStatistc& header);

   void abort(MPIStatistc& header);
   int getNumRanks() const
   {
      return nRanks;
   }
   int getRank() const
   {
      return rank;
   }
 protected:
   bool openSend;
   MPI_Comm const comm;
   const unsigned maxThreads;
   int rank;
   int nRanks;
   uint64_t sendBufferPos;
   uint64_t const sendBufferSize;

   AtomicRingAllocator mpiRecvClauses;
   MPIStatistc accHeader;
   MPI_Request sendRequest;
   uint8_t * sendBuffer;
   uint8_t * recvBuffer;

   void send(MPIStatistc const & header, uint64_t const nBytes);

   bool tryFinishSend();
   void finishSend();
   void updateHeader();

};

inline bool MPIBroadcastConnector::shouldImport(size_type const & pos) const
{
   return AtomicConnector::shouldImport(pos.localPos)
      || mpiRecvClauses.bytesToEnd(pos.globalPos) > 0.2 * mpiRecvClauses.capacity();
}

void MPIBroadcastConnector::updateHeader()
{
   assert(isAllowedToSend());
   accHeader.reset();
   for (int i = 0; i < nRanks; ++i)
   {
      MPIStatistc const & h = *reinterpret_cast<MPIStatistc const*>(recvBuffer + sendBufferSize * i);
      accHeader.add(h);
   }
}

inline bool MPIBroadcastConnector::isValid(const size_type& pos) const
{
   return AtomicConnector::ara.isValid(pos.localPos) || mpiRecvClauses.isValid(pos.globalPos);
}

inline MPIBroadcastConnector::size_type MPIBroadcastConnector::next(const size_type& pos)
{
   assert(isValid(pos));
   size_type res(pos);
   if (AtomicConnector::ara.isValid(pos.localPos))
      res.localPos = AtomicConnector::ara.getNextPos(pos.localPos);
   else
      res.globalPos = mpiRecvClauses.getNextPos(pos.globalPos);
   return res;
}

inline void MPIBroadcastConnector::progress()
{
   if (openSend)
      tryFinishSend();
}

bool MPIBroadcastConnector::tryFinishSend()
{
   if (openSend)
   {
      int flag;
      MPI_Test(&sendRequest, &flag, MPI_STATUS_IGNORE);
      if (flag)
      {
         openSend = false;
         updateHeader();
         if (accHeader.abort)
            NoConnector::abort();
      } else
         NoConnector::sleep();
   }
   return !openSend;
}

inline void MPIBroadcastConnector::finishSend()
{
   while (!tryFinishSend())
      NoConnector::sleep();
}

inline bool MPIBroadcastConnector::isAllowedToSend() const
{
   return !openSend && !accHeader.abort;
}

void MPIBroadcastConnector::abort(MPIStatistc& header)
{
   if (openSend)
      finishSend();
   if (!accHeader.abort)
   {
      header.abort = true;
      send(header);
   }
}

void MPIBroadcastConnector::send(const MPIStatistc& header)
{
   send(header, sendBufferPos);
   sendBufferPos = sizeof(MPIStatistc);
}
void MPIBroadcastConnector::send(MPIStatistc const & header, uint64_t const nBytes)
{
   assert(isAllowedToSend());
   assert(nBytes <= sendBufferSize);
   int flag;
   void const * hData = &header;
   memcpy(sendBuffer, hData, sizeof(MPIStatistc));
   MPI_Iallgather(sendBuffer, nBytes, MPI_BYTE, recvBuffer, sendBufferSize,
   MPI_BYTE,
                  comm, &sendRequest);
   MPI_Test(&sendRequest, &flag, MPI_STATUS_IGNORE);
   openSend = true;

}

inline MPIBroadcastConnector::MPIBroadcastConnector(
                                                    uint64_t bytesPerThread,
                                                    uint64_t bytesPerRank,
                                                    unsigned const nThreads,
                                                    unsigned const maxThreads,
                                                    MPI_Comm const comm)
      : AtomicConnector(bytesPerThread, nThreads),
        openSend(false),
        comm(comm),
        maxThreads(maxThreads),
        rank(-1),
        nRanks(-1),
        sendBufferPos(sizeof(MPIStatistc)),
        // we only use eights of the buffer size, to give solver a chance to read from mpiRecvClauses
        sendBufferSize(bytesPerRank),
        mpiRecvClauses(),
        accHeader(),
        sendRequest(),
        sendBuffer(NULL),
        recvBuffer(NULL)
{
   static bool initFirst = true;
   assert(initFirst);
   initFirst = false;

   MPI_Comm_size(comm, &nRanks);
   MPI_Comm_rank(comm, &rank);
   if (bytesPerThread > 0)
      mpiRecvClauses.growTo(bytesPerRank * 3 * nRanks, nThreads);
   sendBuffer = allocMpiMem<uint8_t>(sendBufferSize);
   recvBuffer = allocMpiMem<uint8_t>(sendBufferSize * nRanks);
}

inline MPIBroadcastConnector::~MPIBroadcastConnector()
{
   if (openSend)
      finishSend();
   freeMpiMem(sendBuffer);
   freeMpiMem(recvBuffer);
}

inline MPIStatistc MPIBroadcastConnector::getAccumulatedHeader()
{
   return accHeader;
}

}
#endif /* SOURCES_PARALLEL_MPIBROADCASTCONNECTOR_H_ */
