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

#ifndef SOURCES_CORE_IMPLICATIONGRAPH_H_
#define SOURCES_CORE_IMPLICATIONGRAPH_H_

#include "mtl/Vec.h"
#include <cassert>
#include <iostream>
#include "mtl/Sort.h"

namespace ctsat
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

   lbool value(Var const x) const;       // The current value of a variable.
   lbool value(Lit const p) const;       // The current value of a literal.

   CRef reason(Var const v) const;
   CRef & reason(Var const v);

   CRef reason(Lit const l) const;
   CRef & reason(Lit const l);

   int level(Var const x) const;
   int level(Lit const l) const;

   template <typename AccessType>
   uint32_t abstractLevel(AccessType const x) const;  // Used to represent an abstraction of sets of decision levels.

   int levelEnd(int const lvl) const;

   int decisionLevel() const;  // Gives the current decisionlevel.
   void newDecisionLevel();                                          // Begins a new decision level.

   int nAssigns() const;
   int nVars() const;

//   void checkTrail() const
//   {
//      int i = trail.size();
//      while(--i >= trail_lim[0])
//         assert(reason(trail[i]) != Database::npos() || level(trail[i]) == 0 || i == trail_lim[level(trail[i])-1]);
//   }

   void assign(Lit const l, CRef const from, int const lvl);
   void assign(Lit const l, CRef const from);
   void assign(Lit const l);
   void unassign(Var const v);

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
   void markSeen2(LitVec const & c, int const startIdx = 0);

   template <typename LitVec>
   bool removeNotSeen2(LitVec & c, int const startIdx = 0);

   template <typename LitVec>
   void printClause(const LitVec & c) const;

   template <class V>
   int computeLBD(const V& c);

   bool satisfied(const Clause& c) const;
   bool locked(const Clause& c) const;
   bool locked(CRef const & cr) const;
   // returns true if clause is satisfied
   bool removeRedundant(vec<Lit> & c) const;
   // returns true if clause is satisfied
   template <typename ClauseType>
   bool removeSetLits(ClauseType & c) const;

   template <typename ClauseVisited, typename LitVisited, typename RunFurther>
   int visitImplications(
                         int index,
                         CRef cref,
                         const int lvl,
                         ClauseVisited const & clauseVisited,
                         LitVisited const & litVisited,
                         RunFurther const & runFurther,
                         bool const unmarkSeen = true);

   bool isImplied(Lit const l) const;
   template <typename LitVec>
   void minimizeImplied(LitVec & c, int const startIdx = 1) const;

   bool isImpliedRecursive(
                           Lit const l,
                           uint32_t const abstract_levels = std::numeric_limits<uint32_t>::max());
   template <typename LitVec>
   void minimizeImpliedRecursive(LitVec & c, bool const full = false, int const startIdx = 1);

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

   vec<Lit> analyze_stack;

   vec<lbool> assigns;          // The current assignments.
   vec<Lit> trail;   // Assignment stack; stores all assigments made in the order they were made.
   vec<int> trail_lim;        // Separator indices for different decision levels in 'trail'.
   vec<VarData> vardata;          // Stores reason and level for each variable.

}
;

template <typename Database>
template <typename AccessType>
inline uint32_t ImplicationGraph<Database>::abstractLevel(AccessType const x) const
{
   return 1 << (level(x) & 31);
}

template <typename Database>
template <typename LitVec>
inline void ImplicationGraph<Database>::minimizeImplied(LitVec & c, int const startIdx) const
{
   int i = startIdx, j = startIdx;

   for (i = j = 1; i < c.size(); i++)
      if (!isImplied(c[i]))
         c[j++] = c[i];
   c.shrink(i - j);
}

template <typename Database>
template <typename LitVec>
inline void ImplicationGraph<Database>::minimizeImpliedRecursive(
                                                          LitVec & c,
                                                          bool const full,
                                                          int const startIdx)
{
   int i = startIdx, j = startIdx;
   uint32_t abstract_level = 0;
   if (full)
      abstract_level = std::numeric_limits<uint32_t>::max();
   else
      for (i = startIdx; i < c.size(); i++)
         abstract_level |= abstractLevel(c[i]);  // (maintain an abstraction of levels involved in conflict)

   for (i = j = startIdx; i < c.size(); i++)
      if (!isImpliedRecursive(c[i], abstract_level))
         c[j++] = c[i];
   c.shrink(i - j);
}

template <typename Database>
template <typename ClauseVisited, typename LitVisited, typename RunFurther>
inline int ImplicationGraph<Database>::visitImplications(
                                                  int index,
                                                  CRef cref,
                                                  const int lvl,
                                                  ClauseVisited const & clauseVisited,
                                                  LitVisited const & litVisited,
                                                  RunFurther const & runFurther,
                                                  bool const unmarkSeen)
{
   Lit p = Lit::Undef();
   // Generate conflict clause:
   //
   assert(lvl == level((ca[cref][0].var())));
   do
   {
      assert(cref != Database::npos());  // (otherwise should be UIP)
      void * debugP = &ca[cref];
      if (clauseVisited(p, cref))
      {
         Clause & c = ca[cref];
         assert(debugP==&c); // FIXME
         // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
         if (p != Lit::Undef() && c.size() == 2 && value(c[0]).isFalse())
         {
            assert(value(c[1]).isTrue());
            Lit const tmp = c[0];
            c[0] = c[1], c[1] = tmp;
         }

         for (int j = (p == Lit::Undef()) ? 0 : 1; j < c.size(); j++)
         {
            Var const v = c[j].var();

            if (!isSeen(v) && level(v) > 0)
            {
               setSeen(v);
               if (litVisited(c[j]))
                  markSeenToClear(c[j]);
            }
         }
      }
      // Select next clause to look at:
      do
      {
         while (!isSeen(getTrailLit(index--).var()))
            ;
         p = getTrailLit(index + 1);
      } while (level((p.var())) < lvl);

      Var const pv = p.var();
      cref = reason(pv);
      if(unmarkSeen)
         unsetSeen(pv);
   } while (runFurther(p));
   return index + 1;
}

template <typename Database>
template <typename LitVec>
void ImplicationGraph<Database>::printClause(const LitVec & c) const
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

template <typename Database>
inline bool ImplicationGraph<Database>::isImplied(Lit const l) const
{
   assert(isSeen(l.var()));
   bool res = false;
   CRef const ref = reason(l.var());
   if (ref != Database::npos())
   {
      int i = 0;
      Clause& c = ca[ref];
      for (; i < c.size(); ++i)
      {
         Var const v = c[i].var();
         if (!isSeen(v) && level(v) > 0)
            break;
      }
      res = i == c.size();
   }
   return res;
}

template <typename Database>
inline bool ImplicationGraph<Database>::isImpliedRecursive(
                                                           Lit const l,
                                                           uint32_t const abstract_levels)
{
   if (reason(l.var()) == Database::npos())
      return false;
   analyze_stack.clear();
   analyze_stack.push(l);
   int const top = nToClear();
   while (analyze_stack.size() > 0)
   {
      assert(isSeen(analyze_stack.last().var()));
      assert(reason((analyze_stack.last().var())) != Database::npos());
      Clause& c = ca[reason((analyze_stack.last().var()))];
      analyze_stack.pop();

      // Special handling for binary clauses like in 'analyze()'.
      if (c.size() == 2 && value(c[0]).isFalse())
      {
         assert(value(c[1]).isTrue());
         Lit const tmp = c[0];
         c[0] = c[1], c[1] = tmp;
      }

      for (int i = 1; i < c.size(); i++)
      {
         Lit const p = c[i];
         Var const v = c[i].var();
         if (!isSeen(v) && level(v) > 0)
         {
            if (reason(v) != Database::npos() && (abstractLevel(v) & abstract_levels) != 0)
            {
               setSeen(v);
               markSeenToClear(p);
               analyze_stack.push(p);

            } else
            {
               clearSeen(top);
               return false;
            }
         }
      }
   }
   return true;
}

template <typename Database>
template <class V>
int ImplicationGraph<Database>::computeLBD(const V& c)
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
template <typename LitVec>
void ImplicationGraph<Database>::markSeen2(LitVec const & c, int const startIdx)
{
   initSeen2();
   for (int i = startIdx; i < c.size(); ++i)
      setSeen2(c[i].var());
}
template <typename Database>
template <typename LitVec>
bool ImplicationGraph<Database>::removeNotSeen2(LitVec & c, int const startIdx)
{
   int i = startIdx, j = startIdx;
   for (; i < c.size(); ++i)
   {
      c[j] = c[i];
      j += isSeen2(c[i].var());
   }
   c.shrink(i - j);
   return i != j;
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
bool ImplicationGraph<Database>::removeRedundant(vec<Lit> & c) const
{
   // Check if clause is satisfied and remove false/duplicate literals:
   sort(c);
   Lit p;
   int i = 0, j = 0;
   for (p = Lit::Undef(); i < c.size(); i++)
      if (value(c[i]).isTrue() || c[i] == ~p)
         return true;
      else if (!value(c[i]).isFalse() && c[i] != p)
         c[j++] = p = c[i];
   c.shrink(i - j);
   return false;
}

template <typename Database>
template <typename ClauseType>
bool ImplicationGraph<Database>::removeSetLits(ClauseType & c) const
{
   // Check if clause is satisfied and remove false/duplicate literal
   int i = 0, j = 0;
   for (; i < c.size(); ++i)
   {
      lbool const val = value(c[i]);
      if (val.isTrue())
         return true;
      else if (!val.isFalse())
         c[j++] = c[i];
   }
   c.shrink(i - j);
   return false;
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
inline typename ImplicationGraph<Database>::lbool ImplicationGraph<Database>::value(
                                                                                    Var const x) const
{
   return assigns[x];
}
template <typename Database>
inline typename ImplicationGraph<Database>::lbool ImplicationGraph<Database>::value(
                                                                                    Lit const p) const
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
inline typename ImplicationGraph<Database>::CRef ImplicationGraph<Database>::reason(
                                                                                    Var const x) const
{
   return vardata[x].reason;
}
template <typename Database>
inline typename ImplicationGraph<Database>::CRef & ImplicationGraph<Database>::reason(Var const x)
{
   return vardata[x].reason;
}
template <typename Database>
inline typename ImplicationGraph<Database>::CRef ImplicationGraph<Database>::reason(
                                                                                    Lit const l) const
{
   return reason(l.var());
}
template <typename Database>
inline typename ImplicationGraph<Database>::CRef & ImplicationGraph<Database>::reason(Lit const l)
{
   return reason(l.var());
}

template <typename Database>
inline int ImplicationGraph<Database>::level(Var const x) const
{
   return vardata[x].level;
}
template <typename Database>
inline int ImplicationGraph<Database>::level(Lit const l) const
{
   return level(l.var());
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
