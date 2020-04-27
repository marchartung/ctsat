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

#ifndef MTL_DRATPRINT_H_
#define MTL_DRATPRINT_H_

#include <cstdio>
#include <string>
#include <iostream>
#include "mtl/Vec.h"

namespace ctsat
{

template <typename Lit>
class DratPrint
{
   template<typename E>
   DratPrint(DratPrint<E> const &) = delete;
   template<typename E>
   DratPrint& operator=(DratPrint<E> const &) = delete;
 public:
   DratPrint(const bool shouldPrint, const std::string & proofFileName)
         : proof((shouldPrint) ? fopen(proofFileName.c_str(), "wb") : NULL)
   {
   }

   DratPrint()
         : DratPrint<Lit>(false, "")
   {
   }

   template<typename E>
   DratPrint(DratPrint<E>&& in)
   : DratPrint()
   {
      *this = std::move(in);
   }
   template<typename E>
   DratPrint& operator=(DratPrint<E>&& in)
   {
      std::swap(proof,in.proof);
      return *this;
   }

   ~DratPrint()
   {
      if (isActive())
      {
         fclose(proof);
         proof = NULL;
      }
   }

   constexpr bool isActive() const
   {
      return proof != NULL;
   }

   template <typename VecType>
   void addClause(const VecType & c)
   {
      addClause('a', c);
   }

   template <typename VecType>
   void addClauseExcludeLit(const VecType & c, const Lit l)
   {
      if (proof != NULL)
      {
         write('a');
         for (int i = 0; i < c.size(); ++i)
            if (c[i] != l)
               writeLit(c[i]);
         write(0);
      }
   }

   void flush()
   {
      // FIXME maybe use a buffer
   }

   void addEmptyClause()
   {
      addClause('a', vec<Lit>());
   }
//-11479 17850
   template <typename VecType>
   void removeClause(const VecType & c)
   {
      addClause('d', c);
   }
   template <typename VecType>
   void addClause(const int prefix, const VecType & c)
   {
      if (isActive())
      {
         write(prefix);
         for (int i = 0; i < c.size(); ++i)
            writeLit(c[i]);
         write(0);
      }
   }
   inline void writeLit(const Lit in)
   {
      unsigned l = in.toInt() + 2;
      assert(l > 0u && in.toInt() + 2 > 0);
      while (l > 127u)
      {
         write(128u + (l & 127u));
         l >>= 7u;
      }
      write(l);
   }

 private:
   FILE * proof;

   inline void write(const unsigned & c)
   {
      assert(proof != NULL);
      putc_unlocked((unsigned char) c, proof);
   }
};
}

#endif /* MTL_DRATPRINT_H_ */

