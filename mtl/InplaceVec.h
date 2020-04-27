/*******************************************************************************************[Vec.h]
 Copyright (c) 2003-2007, Niklas Een, Niklas Sorensson
 Copyright (c) 2007-2010, Niklas Sorensson

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef Concusat_InplaceVec_h
#define Concusat_InplaceVec_h

#include <assert.h>
#include <new>

#include <cstdint>
#include "mtl/XAlloc.h"
#include <string.h>
#include <algorithm>
#include <iostream>

namespace ctsat
{

//=================================================================================================
// Automatically resizable arrays
//
// NOTE! Don't use this InplaceVector on datatypes that cannot be re-located in memory (with realloc)

template <class T>
struct InplaceStorage
{
   int sz;
   int cap;
   T data[1];

   static int headSize()
   {
      return sizeof(int) * 2;
   }
};

template <class T>
class InplaceVec
{
   typedef InplaceStorage<T> VType;
   VType* _data;

   // Don't allow copying (error prone):
   InplaceVec<T>& operator =(const InplaceVec<T>& other) = delete;

   //static inline void nextCap(int& cap){ cap += ((cap >> 1) + 2) & ~1; }
   static inline void nextCap(int& cap)
   {
      cap += ((cap >> 1) + 2) & ~1;
   }

   InplaceVec(InplaceVec<T> const &) = delete;
   InplaceVec<T> & operator=(InplaceVec<T> const &) = delete;

 public:
   // Constructors:
   InplaceVec()
         : _data(nullptr)
   {
   }
   // Helpers for calculating next capacity:
   static inline int imax(int x, int y)
   {
      int mask = (y - x) >> (sizeof(int) * 8 - 1);
      return (x & mask) + (y & (~mask));
   }

   InplaceVec(InplaceVec<T> && in)
         : InplaceVec<T>()
   {
      std::swap(_data, in._data);
   }

   explicit InplaceVec(int size)
         : InplaceVec<T>()
   {
      growTo(size);
   }

   InplaceVec(int size, const T& pad)
         : InplaceVec<T>()
   {
      growTo(size, pad);
   }

   ~InplaceVec()
   {
      clear(true);
   }

   bool isNull() const
   {
      return _data == nullptr;
   }

   InplaceVec<T> & operator=(InplaceVec<T> && in)
   {
      std::swap(_data, in._data);
      return *this;
   }

   // Pointer to first element:
   operator T*(void)
   {
      return data();
   }

   inline void * plainData()
   {
      return _data;
   }

   inline T* data()
   {
      return _data->data;
   }
   inline const T* data() const
   {
      return _data->data;
   }

   inline uint64_t dynBytes() const
   {
      return (_data == nullptr) ? 0 : sizeof(T) * (_data->cap);
   }

   // Size operations:
   inline int size(void) const
   {
      return (_data == nullptr) ? 0 : _data->sz;
   }
   inline void shrink(int nelems)
   {
      if (nelems > 0)
      {
      assert(nelems <= size());
      for (int i = 0; i < nelems; i++)
         --_data->sz, _data->data[_data->sz].~T();
      }
   }
   inline void shrink_(int nelems)
   {
      if (nelems > 0)
      {
         assert(nelems <= size());
         _data->sz -= nelems;
      }
   }
   inline int capacity(void) const
   {
      return (_data == nullptr) ? 0 : _data->cap;
   }

   inline void resize_(int nelems)
   {
      assert(nelems <= capacity());
      if (size() != nelems)
         _data->sz = nelems;
   }

   template <typename VecType>
   InplaceVec<T> & intersect(const VecType & c)
   {
      int i = 0, j = 0;
      for (; i < size(); ++i)
      {
         for (int k = 0; k < c.size(); ++k)
            if (c[k] == _data->data[i])
            {
               _data->data[j++] = _data->data[i];
               break;
            }
      }
      shrink_(i - j);
      return *this;
   }

   void capacity(int min_cap);
   void growTo(int size);
   void lazy_growTo(int size);
   void growTo(int size, const T& pad);
   void clear(bool dealloc = false);
   void unordered_remove(const int & idx)
   {
      assert(size() > idx);
      _data->data[idx] = _data->data[--_data->sz];
   }
   void erase(const int & startIdx, const int & endIdx)
   {
      assert(endIdx <= size());
      if (startIdx == 0 && endIdx == size())
         clear(false);
      else
      {
         int i = startIdx, j = endIdx;
         for (; j < size(); ++i, ++j)
         {
            _data->data[i] = _data->data[j];
            _data->data[j].~T();
         }
         shrink_(endIdx - startIdx);
      }
   }

   inline void erase(const InplaceVec<int> & idxs)
   {
      int i = idxs.size();
      while (i-- > 0)
      {
         assert(i == 0 || idxs[i] > idxs[i - 1]);  // idxs must be sorted increasing
         this->unordered_remove(idxs[i]);
      }
   }

   bool contains(const T & elem)
   {
      bool res = false;
      for (int i = 0; i < size(); ++i)
         if (elem == _data->data[i])
         {
            res = true;
            break;
         }
      return res;
   }

   // Stack interface:
   void push(void)
   {
      if (_data == nullptr || _data->sz == _data->cap)
         capacity(_data->sz + 1);
      new (&_data->data[_data->sz]) T();
      ++_data->sz;
   }
   void push(const InplaceVec<T> & in)
   {
      if (size() == capacity())
         capacity(size() + in.size());
      for (int i = size(); i < size() + in.size(); ++i)
         _data->data[i] = in[i - size()];
      _data->sz += in.size();
   }
   inline void push(const T& elem)
   {
      if (size() == capacity())
         capacity(size() + 1);
      push_(elem);
   }

   void push(T&& elem)
   {
      if (size() == capacity())
         capacity(size() + 1);
      push_(std::move(elem));
   }

   void push_construct(const T& elem)
   {
      if (size() == capacity())
         capacity(size() + 1);
      new (&_data->data[size()]) T(elem);
      _data->sz++;
   }
   void push_(const T& elem)
   {
      assert(size() < capacity());
      _data->data[_data->sz++] = elem;
   }
   void push_(T&& elem)
   {
      assert(size() < capacity());
      _data->data[_data->sz++] = std::move(elem);
   }
   void pop(void)
   {
      assert(size() > 0);
      --_data->sz, _data->data[_data->sz].~T();
   }
   // NOTE: it seems possible that overflow can happen in the 'sz+1' expression of 'push()', but
   // in fact it can not since it requires that 'cap' is equal to INT_MAX. This in turn can not
   // happen given the way capacities are calculated (below). Essentially, all capacities are
   // even, but INT_MAX is odd.

   inline const T& last(void) const
   {
      return _data->data[size() - 1];
   }
   inline T& last(void)
   {
      return _data->data[size() - 1];
   }

   // Vector interface:
   inline const T& operator [](int index) const
   {
      return _data->data[index];
   }
   inline T& operator [](int index)
   {
      return _data->data[index];
   }

   // Duplicatation (preferred instead):
   void copyTo(InplaceVec<T>& copy) const
   {
      copy.clear();
      copy.lazy_growTo(size());
      for (int i = 0; i < size(); i++)
         copy[i] = _data->data[i];
   }

   void moveTo(InplaceVec<T>& dest)
   {
      dest.clear(true);
      dest._data = _data;
      _data = nullptr;
   }
   void memCopyTo(InplaceVec<T>& copy) const
   {
      copy.capacity(capacity());
      copy._data->sz = size();
      memcpy(copy._data, _data, sizeof(T) * size());
   }

   void swap(InplaceVec<T> & swappy)
   {
      std::swap(_data, swappy._data);
   }

   bool hasUniqueElements() const
   {
      bool res = true;
      for (int i = 0; i < size(); ++i)
         for (int j = i + 1; j < size(); ++j)
            if (_data->data[i] == _data->data[j])
            {
               res = false;
               break;
            }
      return res;
   }

   void shrinkClever()
   {

   }

};

//template <class T>
//void InplaceVec<T>::capacity(int min_cap)
//{
//   if (cap >= min_cap)
//      return;
//   int add = imax((min_cap - cap + 1) & ~1, ((cap >> 1) + 2) & ~1);  // NOTE: grow by approximately 3/2
//   cap += add;
//   _data = reinterpret_cast<T*>(xrealloc(_data, cap * sizeof(T)));
//}

template <class T>
void InplaceVec<T>::capacity(int min_cap)
{
   if (capacity() >= min_cap)
      return;
   int sz = size(), cap = capacity(), add = imax((min_cap - cap + 1) & ~1, ((cap >> 1) + 2) & ~1);  // NOTE: grow by approximately 3/2
   int newMem = VType::headSize() + (cap += add) * sizeof(T);
   if (add > INT_MAX - capacity()
      || (((_data = (VType*) Glucose::xrealloc(_data, newMem)) == NULL) && errno == ENOMEM))
      throw OutOfMemoryException();
   _data->sz = sz;
   _data->cap = cap;
}

template <class T>
void InplaceVec<T>::growTo(int s, const T& pad)
{
   if (size() >= s)
      return;
   capacity(s);
   for (int i = size(); i < s; i++)
      _data->data[i] = pad;
   _data->sz = s;
}

template <class T>
void InplaceVec<T>::growTo(int s)
{
   if (size() >= s)
      return;
   capacity(s);
   for (int i = size(); i < s; i++)
      new (&_data[i]) T();
   _data->sz = s;
}

template <class T>
void InplaceVec<T>::lazy_growTo(int s)
{
   if (size() >= s)
      return;
   capacity(s);
   _data->sz = s;
}

template <class T>
void InplaceVec<T>::clear(bool dealloc)
{
   if (_data != NULL)
   {
      for (int i = 0; i < size(); i++)
         _data->data[i].~T();
      _data->sz = 0;
      if (dealloc)
         Glucose::xfree(_data), _data = NULL;
   }
}

//=================================================================================================
}

#endif
