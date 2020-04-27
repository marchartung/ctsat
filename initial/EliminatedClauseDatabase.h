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

#ifndef SOURCES_SIMP_ELIMINATEDCLAUSEDATABASE_H_
#define SOURCES_SIMP_ELIMINATEDCLAUSEDATABASE_H_

#include <cstdint>
#include "mtl/Vec.h"
#include "database/MinisatAllocatorDb.h"

namespace ctsat
{

class EliminatedClauseDatabase
{
   typedef ClauseAllocator::Clause Clause;
   typedef ClauseAllocator::Lit Lit;
   typedef ClauseAllocator::lbool lbool;
 public:

   void addElimClause(Var const v, Clause const & c);
   void addElimUnit(Lit const & l);

   void printModel(vec<lbool> & model) const;
   void extendModel(vec<lbool> & model) const;

   double mBytes() const;

 private:
   vec<uint32_t> elimclauses;

   void mkElimClause(vec<uint32_t>& elimclauses, Lit const x) const;
   void mkElimClause(vec<uint32_t>& elimclauses, Var const v, Clause& c) const;
};

inline double EliminatedClauseDatabase::mBytes() const
{
   return static_cast<double>(elimclauses.size() * sizeof(uint32_t)) / (1024.0 * 1024.0);
}

inline void EliminatedClauseDatabase::addElimUnit(const Lit& l)
{
   elimclauses.push(l.toInt());
   elimclauses.push(1);
}

class EliminatedClauseDatabase1
{
   typedef ClauseAllocator::Clause Clause;
   typedef ClauseAllocator::Lit Lit;
   typedef ClauseAllocator::lbool lbool;
 public:

   void addElimClause(Var const v, Clause const & c)
   {
      if(v == 284)
         std::cout << "hoh\n";
      if (v != curV)
      {
         assert(curV == var_Undef);
         curV = v;
         if (vElimClauses.size() <= v)
            vElimClauses.growTo(v + 1);
         vElimClauses[v].pos = elimclauses.size();
         order.push(v);
      }
      elimclauses.push(c.size() - 1);
      for (int i = 0; i < c.size(); ++i)
         if (c[i].var() != v)
            elimclauses.push(c[i].toInt());
   }
   void addElimUnit(Lit const & l)
   {
      assert(l.var() == curV);
      vElimClauses[l.var()].backUpSign = !l.sign();
      curV = var_Undef;
   }
   void extendModel(vec<lbool> & model) const
   {
      int lastPos = elimclauses.size();
      for (int i = order.size() - 1; i >= 0; ++i)
      {
         Var const v = order[i];
         assert(model[v].isUndef());
         const bool sign = vElimClauses[v].backUpSign;
         int pos = vElimClauses[v].pos;
         while (pos < lastPos)
         {
            bool sat = false;
            int const sz = elimclauses[pos++];
            assert(pos + sz <= lastPos);
            for (int k = 0; k < sz; ++k)
            {
               Lit const l = Lit::toLit(elimclauses[pos + k]);
               assert(!model[l.var()].isUndef());
               if (model[l.var()].sign() == l.sign())
               {
                  sat = true;
                  break;
               }
            }
            if(!sat)
            {
               assert(model[v].isUndef() || model[v] == lbool(sign));
               model[v] = lbool(sign);
            }
         }

         if(model[v].isUndef())
            model[v] = lbool(!sign);
         lastPos = vElimClauses[v].pos;
      }

   }

   void printModel(vec<lbool> & model) const
   {
      printf("v ");
      for (int i = 0; i < model.size(); i++)
         if (!model[i].isUndef())
            printf("%s%s%d", (i == 0) ? "" : " ", (model[i].isTrue()) ? "" : "-", i + 1);
      printf(" 0\n");
   }

   double mBytes() const;

 private:
   Var curV;
   struct ElimVar
   {
      bool backUpSign;
      int pos;
   };
   vec<Var> order;
   vec<ElimVar> vElimClauses;
   vec<uint32_t> elimclauses;
};

}

#endif /* SOURCES_SIMP_ELIMINATEDCLAUSEDATABASE_H_ */
