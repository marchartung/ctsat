/*
 * StealQueue.h
 *
 *  Created on: 24.04.2020
 *      Author: hartung
 */

#ifndef MPI_ATOMICSTEALRESULTQUEUE_H_
#define MPI_ATOMICSTEALRESULTQUEUE_H_

#include <atomic>
#include <limits>
#include <unistd.h>

#include "mtl/Vec.h"

namespace ctsat
{

template <typename StealType, typename ResultType>
class AtomicStealResultQueue
{
 public:
   typedef int size_type;
   AtomicStealResultQueue();

   static size_type npos();

   // master functions;
   int size() const;

   void reset();
   void add(ResultType const & elem);
   void openQueue();
   void wait();

// worker functions
   size_type trySteal();
   StealType const & getStolen(size_type const id);
   void setResult(size_type const id, ResultType const & in);

 private:
   struct Bounds
   {
      size_type stealPos;
      size_type nElems;

      void set(size_type const startPos, size_type const n)
      {
         stealPos = startPos;
         nElems = n;
      }
   };
   std::atomic<size_type> numResults;
   std::atomic<Bounds> bounds;
   vec<StealType> queue;
   vec<ResultType> result;
};

template <typename StealType, typename ResultType>
AtomicStealResultQueue<StealType, ResultType>::size_type AtomicStealResultQueue<StealType, ResultType>::npos()
{
   return std::numeric_limits<size_type>::max();
}

template <typename StealType, typename ResultType>
inline AtomicStealResultQueue<StealType, ResultType>::AtomicStealResultQueue()
      : numResults(0)
{
   bounds.init(0, 0);
}

template <typename StealType, typename ResultType>
inline int AtomicStealResultQueue<StealType, ResultType>::size() const
{
   return queue.size();
}


template <typename StealType, typename ResultType>
inline void AtomicStealResultQueue<StealType, ResultType>::reset()
{
   queue.clear();
   result.clear();
}

template <typename StealType, typename ResultType>
inline void AtomicStealResultQueue<StealType, ResultType>::add(const ResultType& elem)
{
   queue.push(elem);
}

template <typename StealType, typename ResultType>
inline void AtomicStealResultQueue<StealType, ResultType>::openQueue()
{
   result.growTo(size());
   numResults = 0;
}

template <typename StealType, typename ResultType>
inline void AtomicStealResultQueue<StealType, ResultType>::wait()
{
   while (numResults < bounds.nElems)
      sleep(100);
}

template <typename StealType, typename ResultType>
inline AtomicStealResultQueue<StealType, ResultType>::size_type AtomicStealResultQueue<StealType,
      ResultType>::trySteal()
{
   Bounds exp, desir;
   do
   {
      exp = bounds;
      desir = exp;
      ++desir.stealPos;
   } while (!bounds.compare_exchange_weak(exp, desir) || desir.stealPos < desir.nElems);

   return (desir.stealPos < desir.nElems) ? desir.stealPos : npos();
}

template <typename StealType, typename ResultType>
inline const StealType& AtomicStealResultQueue<StealType, ResultType>::getStolen(const size_type id)
{
   assert(id < size());
   return queue[id];
}

template <typename StealType, typename ResultType>
inline void AtomicStealResultQueue<StealType, ResultType>::setResult(
                                                                     const size_type id,
                                                                     const ResultType& in)
{
   assert(id < size());
   result[id] = in;
   assert(numResults < bounds.nElems);
   numResults.fetch_add(1);  // important, flushes result write out
}
}

#endif /* MPI_ATOMICSTEALRESULTQUEUE_H_ */
