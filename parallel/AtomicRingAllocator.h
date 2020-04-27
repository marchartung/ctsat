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

#ifndef SOURCES_PARALLEL_ATOMICRINGALLOCATOR_H_
#define SOURCES_PARALLEL_ATOMICRINGALLOCATOR_H_

#include <cstdint>
#include <cstdlib>
#include <tuple>
#include <limits>
#include <atomic>
#include <cstring>
#include <cassert>
#include <iostream>

#include "utils/Exceptions.h"
#include "mtl/XAlloc.h"

namespace ctsat
{

constexpr uint64_t makeMultiple(uint64_t const align, uint64_t const n)
{
   return n + n % align;
}

// does not guard buffer overflow
class AtomicRingAllocator
{
   AtomicRingAllocator(AtomicRingAllocator const &) = delete;
   AtomicRingAllocator(AtomicRingAllocator &&) = delete;
   AtomicRingAllocator & operator=(AtomicRingAllocator const &) = delete;
   AtomicRingAllocator & operator=(AtomicRingAllocator &&) = delete;

 public:

   typedef uint32_t size_type;
   typedef int32_t value_type;

   static_assert(ATOMIC_INT_LOCK_FREE > 0, "This implementation will be too slow");
   static_assert(ATOMIC_LONG_LOCK_FREE > 0, "This is a todo, where a lock version should be written");

   static size_type startPos()
   {
      return 1;
   }

   static size_type offset()
   {
      return 1;
   }

   AtomicRingAllocator(uint64_t const nBytes, unsigned const nThreads);
   AtomicRingAllocator();

   template <typename T, typename ... Args>
   bool allocConstruct(uint64_t const nBytes, Args ... args);
   template <typename T>
   T const & get(size_type const pos) const;

   bool isValid(size_type const pos) const;

   // only call when isValid is true:
   size_type getNextPos(size_type const pos);

   uint64_t bytesToEnd(size_type const pos) const;
   uint64_t capacity() const;

   size_type getEndPos();

   void growTo(uint64_t const & nBytes, unsigned const nThreads);

 private:

   // blocks either the touch to pos = 0 or the touch to pos >= secondBlockedPos
   // steppedOn can only be called once per pos

   /*
    * two states: [0 ############ secondBlockedPos ----------]
    * when  blockedPos == 0: --- is not allowed to be touched
    * when blockedPos = secondBlockedPos: ### is not allowed to be touched
    * TODO the naming is twisted, rename
    */
   struct BlockerPos
   {
      BlockerPos(BlockerPos&&) = delete;
      BlockerPos(BlockerPos const&) = delete;

      bool isIn(size_type const & pos, size_type const blockPos) const
      {
         if (blockPos == 0)
            return pos >= secondBlockedPos;
         else
            return pos < secondBlockedPos;
      }

      size_type nThreads;
      size_type secondBlockedPos;
      std::atomic<size_type> blockedPos;
      std::atomic<size_type> nOverStepped;

      BlockerPos(size_type const nThreads, size_type const secondBlockedPos)
            : nThreads(nThreads),
              secondBlockedPos(secondBlockedPos),
              blockedPos(0),
              nOverStepped(nThreads)
      {
      }

      BlockerPos()
            : nThreads(0),
              secondBlockedPos(0),
              blockedPos(0),
              nOverStepped(0)
      {
      }

      BlockerPos & operator=(BlockerPos const & in)
      {
         nThreads = in.nThreads;
         secondBlockedPos = in.secondBlockedPos;
         blockedPos = in.blockedPos.load();
         nOverStepped = in.nOverStepped.load();
         return *this;
      }

      BlockerPos & operator=(BlockerPos && in)
      {
         nThreads = in.nThreads;
         secondBlockedPos = in.secondBlockedPos;
         blockedPos = in.blockedPos.load();
         nOverStepped = in.nOverStepped.load();
         return *this;
      }

      void stepped(size_type const from, size_type const to)
      {
         assert(secondBlockedPos != 0);
         if (to == 0 || (from < secondBlockedPos && to >= secondBlockedPos))
            ++nOverStepped;
      }
      // not thread safe
      bool safeToTouch(size_type const pos)
      {
         assert(secondBlockedPos != 0);
         size_type const lockPos = blockedPos;
         bool isBlocked = isIn(pos, lockPos);
         if (isBlocked)
         {
            if (nOverStepped >= nThreads)
            {
               nOverStepped.fetch_sub(nThreads);
               isBlocked = false;
               blockedPos = (blockedPos == 0) ? secondBlockedPos : 0;
            }
         }
         return !isBlocked;
      }
   };

   static const uint64_t alignment = sizeof(value_type);
   size_type cap;
   std::atomic<bool> writeLocked;
   std::atomic<size_type> writeEndPos;
   std::atomic<value_type> * data;
   BlockerPos blocker;

   size_type getPosSafe(size_type const pos) const;
};

AtomicRingAllocator::size_type AtomicRingAllocator::getEndPos()
{
   return writeEndPos;
}

bool AtomicRingAllocator::isValid(size_type const pos) const
{
   return data[getPosSafe(pos)] != -1;
}

inline AtomicRingAllocator::size_type AtomicRingAllocator::getPosSafe(size_type const pos) const
{
   assert(pos < cap);
   size_type res = pos;
   if (data[pos].load(std::memory_order_acquire) == 0)
      res = 0;
   return res;
}

AtomicRingAllocator::size_type AtomicRingAllocator::getNextPos(size_type const pos)
{
   assert(isValid(pos));

   size_type res = getPosSafe(pos) + data[pos].load(std::memory_order_relaxed);
   assert(pos != res);
   blocker.stepped(pos, res);
   return res;
}
inline AtomicRingAllocator::AtomicRingAllocator(const uint64_t nBytes, unsigned const nThreads)
      : cap(0),
        writeLocked(false),
        writeEndPos(0),
        data(nullptr)
{
   if (nBytes > 0)
      growTo(nBytes, nThreads);
}

inline AtomicRingAllocator::AtomicRingAllocator()
      : cap(0),
        writeLocked(false),
        writeEndPos(0),
        data(nullptr)
{
}

// not thread safe
void AtomicRingAllocator::growTo(uint64_t const & nBytes, unsigned const nThreads)
{
   assert(data == nullptr);
   assert(nBytes > 500 * sizeof(value_type) && "Buffer size to small");
   cap = makeMultiple(alignment, nBytes) / alignment;
   blocker = BlockerPos(nThreads, cap / 2);
   data = reinterpret_cast<std::atomic<value_type>*>(xaligned_alloc(sizeof(std::atomic<value_type>),
                                                                    capacity()));
   writeEndPos = startPos();
   data[startPos()] = -1;
   assert(data[0].is_lock_free());
}

template <typename T>
bool compareTyped(T const * a, T const * b, unsigned const num)
{
   for (unsigned i = 0; i < num; ++i)
      if (a[i] != b[i])
      {
         assert(false);
         return false;
      }
   return true;
}

template <typename T, typename ... Args>
inline bool AtomicRingAllocator::allocConstruct(uint64_t const nBytes, Args ... args)
{
   size_type const addSize = makeMultiple(alignment, nBytes) / alignment + offset();
   assert(addSize < cap / 2);
   if (addSize >= cap / 2)
   {
      std::cout << "AtomicRingAllocator: Buffer size too small for allocation" << std::endl;
      throw OutOfMemoryException();
   }
   bool expected;
   do
   {
      expected = false;
   } while (!writeLocked.compare_exchange_weak(expected, true));  // TODO maybe switch to something without thread starvation
   assert(writeLocked);

   size_type start = writeEndPos;
   assert(data[start] == -1);
   if (!blocker.safeToTouch(start + addSize))
   {
      writeLocked.store(false, std::memory_order_seq_cst);
      return false;
   }
   if (cap <= start + addSize)
   {
      data[0].store(-1, std::memory_order_relaxed);
      data[start].store(0, std::memory_order_release);
      start = 0;
   }
   data[start + addSize].store(-1, std::memory_order_relaxed);
   writeEndPos.store(start + addSize, std::memory_order_seq_cst);
//   std::cout << "write: [" << start << "," << start+addSize << "]\n";
   writeLocked.store(false, std::memory_order_seq_cst);

   void * pos = reinterpret_cast<void*>(data + (start + offset()));
   if (pos != new (pos) T(args...))
      assert(false && "The current value_type is not the systems alignment");

   data[start].store(addSize, std::memory_order_release);
   return true;
}

template <typename T>
inline T const & AtomicRingAllocator::get(const size_type pos) const
{
   return *reinterpret_cast<T*>(data + (getPosSafe(pos) + offset()));
}

inline uint64_t AtomicRingAllocator::bytesToEnd(const size_type pos) const
{
   size_type res = writeEndPos - pos;
   if (writeEndPos < pos)
      res = cap - pos + writeEndPos;
   return static_cast<uint64_t>(res) * alignment;
}

inline uint64_t AtomicRingAllocator::capacity() const
{
   return static_cast<uint64_t>(cap) * alignment;
}
}

#endif /* SOURCES_PARALLEL_ATOMICRINGALLOCATOR_H_ */
