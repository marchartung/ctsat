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

#ifndef SOURCES_PARALLEL_MPIBROADCASTCONNECTOR_H_
#define SOURCES_PARALLEL_MPIBROADCASTCONNECTOR_H_

#include <cstdint>
#include <tuple>
#include <algorithm>
#include <mpi/mpi.h> // FIXME
#include <unistd.h>

#include "parallel/AtomicConnector.h"
#include "initial/Inputs.h"
#include "MPIStatistc.h"

namespace CTSat
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
              globalPos(0)
      {
      }

      size_type()
            : size_type(0)
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
      return NoConnector::getUniqueId() + rank * Inputs::nThreads;
   }

   MPIBroadcastConnector();
   ~MPIBroadcastConnector();

   bool isRoot() const;
   int getRank() const;
   int getNumRanks() const;
   uint64_t getSendBufferSize() const
   {
      return sendBufferSize - sizeof(MPIStatistc);
   }

   bool isAllowedToSend() const;

   void sendInstance(void const * data, uint64_t const nBytes);
   std::tuple<void const *, uint64_t> receiveInstance();
   void freeInstance();

   void progress();

   MPIStatistc getAccumulatedHeader();

   lbool receiveResult(int const nVars);

   void sendModel();

   void sendResult(MPIStatistc const & header);

   void printState();

   bool isValid(size_type const & pos) const;
   size_type next(size_type const & pos) const;

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
   inline size_type exchange(uint64_t const & nBytes, Args ... args)
   {
      AtomicConnector::size_type res = AtomicConnector::exchange<ClauseType, Args...>(nBytes,
                                                                                      args...);
      return size_type(res);
   }

   inline uint64_t nFreeSendBytes() const
   {
      return sendBufferSize-sendBufferPos;
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

 protected:
   bool openSend;
   MPI_Comm const comm;
   int rank;
   int nRanks;
   int numTimerIncrease;
   uint64_t sendBufferPos;
   uint64_t const sendBufferSize;

   AtomicRingAllocator mpiRecvClauses;
   MPIStatistc accHeader;
   MPI_Request instanceReq;
   MPI_Request sendRequest;
   uint8_t * instanceBuffer;
   uint8_t * sendBuffer;
   uint8_t * recvBuffer;
   uint8_t * modelBuffer;

   template <typename T>
   T * allocMpiMem(uint64_t const nBytes)
   {
      void * res = NULL;
      if (MPI_Alloc_mem(nBytes, MPI_INFO_NULL, &res))
      {
         std::cout << "Error: Could not allocate memory\n";
         MPI_Finalize();
         exit(-1);  // no destructor call on exit?!?
      }
      return reinterpret_cast<T*>(res);
   }

   void freeMpiMem(uint8_t * & p)
   {
      if (p != NULL)
      {
         MPI_Free_mem(p);
         p = NULL;
      }
   }
   void send(MPIStatistc const & header, uint64_t const nBytes);

   bool tryFinishSend();
   void finishSend();
   void updateHeader();

}
;

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

inline MPIBroadcastConnector::size_type MPIBroadcastConnector::next(const size_type& pos) const
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
         if (!accHeader.res.isUndef())
            NoConnector::setFinished(accHeader.res);
         else if (accHeader.abort)
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
   return !openSend && accHeader.res.isUndef() && !accHeader.abort;
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

void MPIBroadcastConnector::sendResult(MPIStatistc const & header)
{
   assert(!header.res.isUndef() || header.abort);
   finishSend();
   if (isAllowedToSend())
   {
      send(header, sizeof(MPIStatistc));
      finishSend();
   }
}

inline lbool MPIBroadcastConnector::receiveResult(int const nVars)
{
   assert(!accHeader.res.isUndef() || accHeader.abort);
   assert(isRoot());
   assert(getResult() == accHeader.res);
   if(accHeader.abort)
      return lbool::Undef();
   if (accHeader.res.isTrue() && accHeader.rank != rank)
   {
      uint64_t const nBytes = nVars * sizeof(lbool);
      modelBuffer = allocMpiMem<uint8_t>(nBytes);
      MPI_Recv(modelBuffer, nBytes, MPI_BYTE, accHeader.rank, 0, comm, MPI_STATUS_IGNORE);

      lbool const * model = reinterpret_cast<lbool const *>(modelBuffer);
      NoConnector::model.clear();
      NoConnector::model.growTo(nVars);
      for (int i = 0; i < nVars; ++i)
         NoConnector::model[i] = model[i];
      freeMpiMem(modelBuffer);

      modelCommited = true;
   }
   return accHeader.res;
}

inline void MPIBroadcastConnector::sendModel()
{
   finishSend();
   assert(accHeader.res.isTrue());
   if (rank == accHeader.rank && rank != 0)
   {
      vec<lbool> &model = getModel();
      uint64_t const nBytes = model.size() * sizeof(lbool);
      modelBuffer = allocMpiMem<uint8_t>(nBytes);

      lbool * m = reinterpret_cast<lbool*>(modelBuffer);
      for (int i = 0; i < model.size(); ++i)
         m[i] = model[i];
      MPI_Send(modelBuffer, nBytes, MPI_BYTE, 0, 0, comm);
      freeMpiMem(modelBuffer);
   }
}

inline MPIBroadcastConnector::MPIBroadcastConnector()
      : openSend(false),
        comm(MPI_COMM_WORLD),
        rank(-1),
        nRanks(-1),
        sendBufferPos(sizeof(MPIStatistc)),
        // we only use eights of the buffer size, to give solver a chance to read from mpiRecvClauses
        sendBufferSize(Inputs::mpiMbBufferSize * (128ull * 1024ull)),
        mpiRecvClauses(),
        accHeader(rank),
        sendRequest(),
        instanceBuffer(NULL),
        sendBuffer(NULL),
        recvBuffer(NULL),
        modelBuffer(NULL)
{
   static bool initFirst = true;
   assert(initFirst);
   initFirst = false;

   MPI_Init(Inputs::argc, Inputs::argv);
   MPI_Comm_size(comm, &nRanks);
   MPI_Comm_rank(comm, &rank);
   mpiRecvClauses.growTo(Inputs::mpiMbBufferSize * 1024ull * 1024ull * nRanks);
   sendBuffer = allocMpiMem<uint8_t>(sendBufferSize);
   recvBuffer = allocMpiMem<uint8_t>(sendBufferSize * nRanks);
}

inline MPIBroadcastConnector::~MPIBroadcastConnector()
{
   if (openSend)
      finishSend();
   freeMpiMem(instanceBuffer);
   freeMpiMem(sendBuffer);
   freeMpiMem(modelBuffer);
   MPI_Finalize();
}

inline bool MPIBroadcastConnector::isRoot() const
{
   return rank == 0;
}

inline int MPIBroadcastConnector::getRank() const
{
   return rank;
}

inline int MPIBroadcastConnector::getNumRanks() const
{
   return nRanks;
}

inline void MPIBroadcastConnector::sendInstance(const void* data, const uint64_t nBytes)
{
   int waitFlag = 0;
// send size:
   uint64_t bytesSend = nBytes;
   MPI_Bcast(&bytesSend, 1, MPI_UINT64_T, 0, comm);  // local completition should be fast
   if (nBytes > 0)
   {
      // send instance
      instanceBuffer = allocMpiMem<uint8_t>(nBytes);
      memcpy(instanceBuffer, data, nBytes);
      MPI_Ibcast(instanceBuffer, nBytes, MPI_BYTE, 0, comm, &instanceReq);
      MPI_Test(&instanceReq, &waitFlag, MPI_STATUS_IGNORE);
   }
}

inline std::tuple<const void*, uint64_t> MPIBroadcastConnector::receiveInstance()
{
// receive number of bytes to allocate:
   uint64_t nBytes = 0;
   MPI_Bcast(&nBytes, 1, MPI_UINT64_T, 0, comm);

   if (nBytes > 0)
   {
      // receive instance:
      instanceBuffer = allocMpiMem<uint8_t>(nBytes);
      MPI_Ibcast(instanceBuffer, nBytes, MPI_BYTE, 0, comm, &instanceReq);
      MPI_Wait(&instanceReq, MPI_STATUS_IGNORE);
   }
   return std::tie(instanceBuffer, nBytes);
}

inline void MPIBroadcastConnector::freeInstance()
{
   if (isRoot())
   {
      MPI_Wait(&instanceReq, MPI_STATUS_IGNORE);
   }
   freeMpiMem(instanceBuffer);
}

inline MPIStatistc MPIBroadcastConnector::getAccumulatedHeader()
{
   return accHeader;
}

}
#endif /* SOURCES_PARALLEL_MPIBROADCASTCONNECTOR_H_ */
