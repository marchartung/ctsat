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

#ifndef SOURCES_PARALLEL_ATOMICRINGALLOCATOR_H_
#define SOURCES_PARALLEL_ATOMICRINGALLOCATOR_H_

#include <cstdint>
#include <cstdlib>
#include <tuple>
#include <limits>
#include <atomic>
#include <cstring>
#include <cassert>

#include "mtl/XAlloc.h"

namespace CTSat
{

constexpr uint64_t makeMultiple(uint64_t const align, uint64_t const n)
{
   return n + n % align;
}

// does not guard buffer overflow
class AtomicRingAllocator
{
 public:

   typedef uint32_t size_type;
   typedef int32_t value_type;

   static_assert(ATOMIC_INT_LOCK_FREE > 0, "This implementation will be too slow");
   static_assert(ATOMIC_LONG_LOCK_FREE > 0, "This is a todo, where a lock version should be written");

   static size_type offset()
   {
      return 1;
   }

   AtomicRingAllocator(uint64_t const nBytes = 0);

   template <typename T, typename ... Args>
   size_type allocConstruct(uint64_t const nBytes, Args ... args);
   template <typename T>
   T const & get(size_type const pos) const;

   bool isValid(size_type const pos) const;

   // only call when isValid is true:
   size_type getNextPos(size_type const pos) const;

   uint64_t bytesToEnd(size_type const pos) const;
   uint64_t capacity() const;

   size_type getEndPos();

   void growTo(uint64_t const & nBytes);

 private:
   static const uint64_t alignment = sizeof(value_type);
   size_type cap;
   std::atomic<bool> writeLocked;
   std::atomic<size_type> writeEndPos;
   std::atomic<value_type> * data;

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

AtomicRingAllocator::size_type AtomicRingAllocator::getNextPos(size_type const pos) const
{
   assert(isValid(pos));

   size_type res = getPosSafe(pos) + data[pos].load(std::memory_order_relaxed);
   assert(pos != res);
   return res;
}
inline AtomicRingAllocator::AtomicRingAllocator(const uint64_t nBytes)
      : cap(0),
        writeLocked(false),
        writeEndPos(0),
        data(nullptr)
{
   if (nBytes > 0)
   {
      growTo(nBytes);
   }
}

// not thread safe
void AtomicRingAllocator::growTo(uint64_t const & nBytes)
{
   assert(data == nullptr);
   assert(nBytes > 500 * sizeof(value_type) && "Buffer size to small");
   cap = makeMultiple(alignment, nBytes) / alignment;
   data = reinterpret_cast<std::atomic<value_type>*>(xaligned_alloc(sizeof(std::atomic<value_type>),
                                                                    capacity()));
   data[0] = -1;
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
inline AtomicRingAllocator::size_type AtomicRingAllocator::allocConstruct(
                                                                          uint64_t const nBytes,
                                                                          Args ... args)
{
   size_type const addSize = makeMultiple(alignment, nBytes) / alignment + offset();
   assert(addSize <= cap);
   bool expected;
   do
   {
      expected = false;
   } while (!writeLocked.compare_exchange_weak(expected, true));  // TODO maybe switch to something without thread starvation
   assert(writeLocked);

   size_type start = writeEndPos;
   assert(data[start] == -1);
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
   return start;
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
