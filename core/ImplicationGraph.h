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

#ifndef SOURCES_CORE_IMPLICATIONGRAPH_H_
#define SOURCES_CORE_IMPLICATIONGRAPH_H_

#include "mtl/Vec.h"
#include <cassert>
#include <iostream>
#include "mtl/Sort.h"

namespace CTSat
{

template <typename Database>
class ImplicationGraph
{
 public:

   typedef typename Database::CRef CRef;
   typedef typename Database::Var Var;
   typedef typename Database::Lit Lit;
   typedef typename Database::lbool lbool;
   typedef typename Database::Clause Clause;

   ImplicationGraph(Database & ca);

   void newVar();

   lbool value(Var x) const;       // The current value of a variable.
   lbool value(Lit p) const;       // The current value of a literal.

   CRef reason(Var const v) const;
   CRef & reason(Var const v);
   int level(Var const x) const;
   uint32_t abstractLevel(Var x) const;  // Used to represent an abstraction of sets of decision levels.

   int levelEnd(int const lvl) const;

   int decisionLevel() const;  // Gives the current decisionlevel.
   void newDecisionLevel();                                          // Begins a new decision level.

   int nAssigns() const;
   int nVars() const;

   void assign(Lit const l, CRef const from, int const lvl);
   void assign(Lit const l, CRef const from);
   void assign(Lit const l);
   void unassign(Var const v);

   // returns true if clause is satisfied
   bool removeRedundant(vec<Lit> & c)
   {
      // Check if clause is satisfied and remove false/duplicate literals:
      sort(c);
      Lit p;
      int i = 0, j = 0;
      for ( p = Lit::Undef(); i < c.size(); i++)
         if (value(c[i]).isTrue() || c[i] == ~p)
            return true;
         else if (!value(c[i]).isFalse() && c[i] != p)
            c[j++] = p = c[i];
      c.shrink(i - j);
      return false;
   }

   // returns true if clause is satisfied
   template<typename ClauseType>
   bool removeSetLits(ClauseType & c)
   {
      // Check if clause is satisfied and remove false/duplicate literal
      int i = 0, j = 0;
      for (; i < c.size(); ++i)
      {
         lbool const val = value(c[i]);
         if (val.isTrue())
            return true;
         else if(!val.isFalse())
            c[j++] = c[i];
      }
      c.shrink(i - j);
      return false;
   }

   void shrink(int const n);

   lbool operator[](Var const idx) const;
   lbool operator[](Lit const l) const;

   Lit getTrailLit(int const idx) const;

   void backtrack(int const lvl);

   void setSeen(Var const v);
   void unsetSeen(Var const v);
   bool isSeen(Var const v) const;

   int nToClear() const;
   void markSeenToClear(Lit const & l);
   void clearSeen(int const keepN = 0);

   void initSeen2();
   void setSeen2(Var const v);
   void unsetSeen2(Var const v);
   bool isSeen2(Var const v) const;

   template <typename LitVec>
   void printClause(const LitVec & c)
   {
      std::cout << "[";
      for (int i = 0; i < c.size(); ++i)
      {
         lbool val = value(c[i]);
         std::cout
            << std::string(((c[i].sign()) ? "-" : ""))
            << c[i].var()
            << "="
            << std::string(
                  ((val == lbool::False()) ? "false" : ((val == lbool::True()) ? "true" : "undef")))
            << "("
            << level(c[i].var())
            << std::string((reason(c[i].var()) != Database::npos()) ? "p" : "d")
            << ") ";
      }
      std::cout << "]\n";
   }
   template <class V> int computeLBD(const V& c)
   {
      int lbd = 0;

      initSeen2();
      for (int i = 0; i < c.size(); i++)
      {
         int l = level(c[i].var());
         if (l != 0 && !isSeen2(l))
         {
            setSeen2(l);
            lbd++;
         }
      }
      return lbd;
   }

   bool satisfied(const Clause& c) const;
   bool locked(const Clause& c) const;
   bool locked(CRef const & cr) const;

 private:
   Database & ca;

   struct VarData
   {
      CRef reason;
      int level;
   };
   static inline VarData mkVarData(CRef cr, int l)
   {
      VarData d = { cr, l };
      return d;
   }

   uint64_t counter;  // Simple counter for marking purpose with 'seen2'.
   vec<Lit> analyze_toclear;
   vec<char> seen;
   vec<uint64_t> seen2;  // Mostly for efficient LBD computation. 'seen2[i]' will indicate if decision level or variable 'i' has been seen.

   vec<lbool> assigns;          // The current assignments.
   vec<Lit> trail;      // Assignment stack; stores all assigments made in the order they were made.
   vec<int> trail_lim;        // Separator indices for different decision levels in 'trail'.
   vec<VarData> vardata;          // Stores reason and level for each variable.

};

template <typename Database> inline uint32_t ImplicationGraph<Database>::abstractLevel(Var x) const
{
   return 1 << (level(x) & 31);
}

template <typename Database>
inline bool ImplicationGraph<Database>::satisfied(const ImplicationGraph<Database>::Clause& c) const
{
   for (int i = 0; i < c.size(); i++)
      if (value(c[i]).isTrue())
         return true;
   return false;
}


template <typename Database>
inline bool ImplicationGraph<Database>::locked(CRef const & cr) const
{
   return locked(ca[cr]);
}

template <typename Database>
inline bool ImplicationGraph<Database>::locked(const ImplicationGraph<Database>::Clause& c) const
{
   int i = c.size() != 2 ? 0 : (value(c[0]).isTrue() ? 0 : 1);
   return value(c[i]).isTrue()
      && reason(c[i].var()) != Database::npos()
      && ca.lea(reason(c[i].var())) == &c;
}

template <typename Database>
inline int ImplicationGraph<Database>::nToClear() const
{
   return analyze_toclear.size();
}

template <typename Database>
inline void ImplicationGraph<Database>::markSeenToClear(Lit const & l)
{
   analyze_toclear.push(l);
}
template <typename Database>
inline void ImplicationGraph<Database>::clearSeen(int const keepN)
{
   for (int i = keepN; i < analyze_toclear.size(); ++i)
      unsetSeen(analyze_toclear[i].var());
   analyze_toclear.shrink(analyze_toclear.size() - keepN);
//   if (keepN == 0)
//      for (int i = 0; i < seen.size(); ++i)
//         assert(!seen[i] || level(i) == 0);
}

template <typename Database>
inline ImplicationGraph<Database>::ImplicationGraph(Database & ca)
      : ca(ca),
        counter(0)
{

}

template <typename Database>
inline void ImplicationGraph<Database>::setSeen(Var const v)
{
   seen[v] = 1;
}
template <typename Database>
inline void ImplicationGraph<Database>::unsetSeen(Var const v)
{
   seen[v] = 0;
}
template <typename Database>
inline bool ImplicationGraph<Database>::isSeen(Var const v) const
{
   return seen[v];
}

template <typename Database>
inline void ImplicationGraph<Database>::initSeen2()
{
   ++counter;
}

template <typename Database>
inline void ImplicationGraph<Database>::setSeen2(Var const v)
{
   seen2[v] = counter;
}
template <typename Database>
inline void ImplicationGraph<Database>::unsetSeen2(Var const v)
{
   seen2[v] = counter - 1;
}
template <typename Database>
inline bool ImplicationGraph<Database>::isSeen2(Var const v) const
{
   return seen2[v] == counter;
}

template <typename Database>
inline void ImplicationGraph<Database>::backtrack(int const lvl)
{
   trail.shrink(trail.size() - trail_lim[lvl]);
   trail_lim.shrink(trail_lim.size() - lvl);
}

template <typename Database>
inline void ImplicationGraph<Database>::unassign(Var const v)
{
   assigns[v] = lbool::Undef();
}

template <typename Database>
inline void ImplicationGraph<Database>::shrink(int const n)
{
   trail.shrink(n);
}

template <typename Database>
inline int ImplicationGraph<Database>::levelEnd(int const lvl) const
{
   return trail_lim[lvl];
}

template <typename Database>
inline typename ImplicationGraph<Database>::Lit ImplicationGraph<Database>::getTrailLit(
                                                                                        int const idx) const
{
   return trail[idx];
}

template <typename Database>
inline void ImplicationGraph<Database>::newVar()
{
   seen.push(0);
   seen2.push(0);
   assigns.push(lbool::Undef());
   vardata.push(mkVarData(Database::npos(), 0));
   trail.capacity(vardata.size() + 1);
}

template <typename Database>
inline typename ImplicationGraph<Database>::lbool ImplicationGraph<Database>::value(Var x) const
{
   return assigns[x];
}
template <typename Database>
inline typename ImplicationGraph<Database>::lbool ImplicationGraph<Database>::value(Lit p) const
{
   return p.value(assigns[p.var()]);
}

template <typename Database>
inline typename ImplicationGraph<Database>::lbool ImplicationGraph<Database>::operator[](
                                                                                         Var const v) const
{
   return value(v);
}

template <typename Database>
inline typename ImplicationGraph<Database>::lbool ImplicationGraph<Database>::operator[](
                                                                                         Lit const l) const
{
   return value(l);
}

template <typename Database>
inline void ImplicationGraph<Database>::assign(Lit const l, CRef const from, int const lvl)
{
   Var const x = l.var();
   vardata[x].level = lvl;
   assign(l, from);
}

template <typename Database>
inline void ImplicationGraph<Database>::assign(Lit const l, CRef const from)
{
   Var const x = l.var();
   vardata[x].reason = from;
   assign(l);
}

template <typename Database>
inline void ImplicationGraph<Database>::assign(Lit const l)
{
   Var const x = l.var();
   assigns[x] = lbool(!l.sign());
   trail.push(l);
}

template <typename Database>
inline typename ImplicationGraph<Database>::CRef ImplicationGraph<Database>::reason(Var x) const
{
   return vardata[x].reason;
}
template <typename Database>
inline typename ImplicationGraph<Database>::CRef & ImplicationGraph<Database>::reason(Var x)
{
   return vardata[x].reason;
}
template <typename Database>
inline int ImplicationGraph<Database>::level(Var x) const
{
   return vardata[x].level;
}
template <typename Database>
inline int ImplicationGraph<Database>::nAssigns() const
{
   return trail.size();
}

template <typename Database>
inline int ImplicationGraph<Database>::nVars() const
{
   return vardata.size();
}

template <typename Database>
inline void ImplicationGraph<Database>::newDecisionLevel()
{
   trail_lim.push(trail.size());
}

template <typename Database>
inline int ImplicationGraph<Database>::decisionLevel() const
{
   return trail_lim.size();
}

} /* namespace Minisat */

#endif /* SOURCES_CORE_IMPLICATIONGRAPH_H_ */
