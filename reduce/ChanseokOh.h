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

#ifndef SOURCES_REDUCE_CHANSEOKOH_H_
#define SOURCES_REDUCE_CHANSEOKOH_H_

#include "core/Statistic.h"
#include "initial/SolverConfig.h"
#include "core/ImplicationGraph.h"

namespace ctsat
{

template <typename Database>
class ChanseokOhReduce
{
   typedef typename Database::Lit Lit;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;

 public:
   ChanseokOhReduce(
                    SolverConfig const & config,
                    Statistic & stat,
                    Database & db,
                    ImplicationGraph<Database> & ig);

   // returns true, if clauses were deleted
   template <typename RemoveFunctor>
   bool run(RemoveFunctor const & remove);

   template <typename RemoveFunctor>
   void removeSatisfied(RemoveFunctor const & remove);

   // applies the given functor to each clause. When the functor returns true, the clause will be removed
   template <typename ApplyFunctor>
   void applyRemoveGoodClauses(ApplyFunctor const & func);

   void addClause(CRef const ref);

   void clauseImproved(CRef const ref);
   void clauseUsedInConflict(CRef const ref);
   void adjustOnConflict();

   void notifyConflictResolved();

   int nClauses() const;

   void relocAll(Database& to);

 public:

   // Don't change the actual numbers.
   static const int LOCAL = 0;
   static const int TIER2 = 2;
   static const int CORE = 3;

   Statistic & stat;
   Database & ca;
   ImplicationGraph<Database> & ig;

   int core_lbd_cut;
   double clause_decay;
   double cla_inc;   // Amount to bump next clause with.
   uint64_t next_T2_reduce, next_L_reduce;

   vec<CRef> learnts_core;
   vec<CRef> learnts_tier2;
   vec<CRef> learnts_local;

   struct reduceDB_lt
   {
      Database & ca;
      reduceDB_lt(Database& ca_)
            : ca(ca_)
      {
      }
      bool operator ()(CRef x, CRef y) const
      {
         return ca[x].activity() < ca[y].activity();
      }
   };
   void claDecayActivity();  // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
   void claBumpActivity(Clause& c);   // Increase a clause with the current 'bump' value.

   template <typename RemoveFunctor>
   void reduceDB(RemoveFunctor const & remove);
   void reduceDB_Tier2();

   template <typename RemoveFunctor>
   void safeRemoveSatisfied(RemoveFunctor const & remove, vec<CRef>& cs, unsigned valid_mark);
};

template <typename Database>
template <typename RemoveFunctor>
void ChanseokOhReduce<Database>::reduceDB(RemoveFunctor const & remove)   // Reduce the set of learnt clauses.
{
   int i, j;
   //if (local_learnts_dirty) cleanLearnts(learnts_local, LOCAL);
   //local_learnts_dirty = false;

   sort(learnts_local, reduceDB_lt(ca));

   int limit = learnts_local.size() / 2;
   for (i = j = 0; i < learnts_local.size(); i++)
   {
      Clause& c = ca[learnts_local[i]];
      if (c.mark() == LOCAL)
         if (c.removable() && !ig.locked(c) && i < limit)
            remove(learnts_local[i]);
         else
         {
            if (!c.removable())
               limit++;
            c.removable(true);
            learnts_local[j++] = learnts_local[i];
         }
   }
   learnts_local.shrink(i - j);
}
template <typename Database>
void ChanseokOhReduce<Database>::reduceDB_Tier2()
{
   int i, j;
   const uint32_t nConfl = stat.conflicts;  // TODO long running props might have overflows
   for (i = j = 0; i < learnts_tier2.size(); i++)
   {
      Clause& c = ca[learnts_tier2[i]];
      if (c.mark() == TIER2)
         if (!ig.locked(c) && c.touched() + 30000 < nConfl)
         {
            learnts_local.push(learnts_tier2[i]);
            c.mark(LOCAL);
            //c.removable(true);
            c.activity() = 0;
            claBumpActivity(c);
         } else
            learnts_tier2[j++] = learnts_tier2[i];
   }
   learnts_tier2.shrink(i - j);
}

template <typename Database>
template <typename RemoveFunctor>
void ChanseokOhReduce<Database>::safeRemoveSatisfied(RemoveFunctor const & remove, vec<CRef>& cs, unsigned valid_mark)
{
   int i, j;
   for (i = j = 0; i < cs.size(); i++)
   {
      Clause& c = ca[cs[i]];
      if (c.mark() == valid_mark)
         if (ig.satisfied(c))
            remove(cs[i]);
         else
            cs[j++] = cs[i];
   }
   cs.shrink(i - j);
}

// returns true, if clauses were deleted
template <typename Database>
template <typename RemoveFunctor>
bool ChanseokOhReduce<Database>::run(RemoveFunctor const & remove)
{
   bool res = false;
   if (stat.conflicts >= next_T2_reduce)
   {
      next_T2_reduce = stat.conflicts + 10000;
      reduceDB_Tier2();
   }
   if (stat.conflicts >= next_L_reduce)
   {
      next_L_reduce = stat.conflicts + 15000;
      reduceDB(remove);
      res = true;
   }
   return res;
}

template <typename Database>
template <typename RemoveFunctor>
void ChanseokOhReduce<Database>::removeSatisfied(RemoveFunctor const & remove)
{
   safeRemoveSatisfied(remove, learnts_core, CORE);  // Should clean core first.
   safeRemoveSatisfied(remove, learnts_tier2, TIER2);
   safeRemoveSatisfied(remove, learnts_local, LOCAL);
}

// applies the given functor to each clause. When the functor returns true, the clause will be removed
template <typename Database>
template <typename ApplyFunctor>
void ChanseokOhReduce<Database>::applyRemoveGoodClauses(ApplyFunctor const & func)
{
   int i = 0, j = 0;
   for (; i < learnts_core.size(); ++i)
      if (!func(learnts_core[i]))
         learnts_core[j++] = learnts_core[i];
   learnts_core.shrink(i - j);
   i = 0;
   j = 0;
   for (; i < learnts_tier2.size(); ++i)
      if (!func(learnts_tier2[i]))
         learnts_tier2[j++] = learnts_tier2[i];
   learnts_tier2.shrink(i - j);
}

template <typename Database>
ChanseokOhReduce<Database>::ChanseokOhReduce(
                                             SolverConfig const & config,
                                             Statistic & stat,
                                             Database & db,
                                             ImplicationGraph<Database> & ig)
      : stat(stat),
        ca(db),
        ig(ig),
        core_lbd_cut(3),
        clause_decay(config.clause_decay),
        cla_inc(1),
        next_T2_reduce(config.next_T2_reduce),
        next_L_reduce(config.next_L_reduce),
        learnts_core(),
        learnts_tier2(),
        learnts_local()
{

}

template <typename Database>
void ChanseokOhReduce<Database>::addClause(CRef const ref)
{

   Clause & c = ca[ref];
   int const lbd = c.lbd();
   if (lbd <= core_lbd_cut)
   {
      learnts_core.push(ref);
      c.mark(CORE);
   } else if (lbd <= 6)
   {
      learnts_tier2.push(ref);
      c.mark(TIER2);
      c.touched() = stat.conflicts;
   } else
   {
      learnts_local.push(ref);
      claBumpActivity(c);
   }
}

template <typename Database>
inline void ChanseokOhReduce<Database>::notifyConflictResolved()
{
   claDecayActivity();
}

template <typename Database>
inline void ChanseokOhReduce<Database>::adjustOnConflict()
{
   if (stat.conflicts == 100000 && learnts_core.size() < 100)
      core_lbd_cut = 5;

}

template <typename Database>
inline void ChanseokOhReduce<Database>::clauseUsedInConflict(CRef const ref)
{
   Clause & c = ca[ref];
   // Update LBD if improved.
   if (c.learnt() && c.mark() != CORE)
   {
      int const lbd = ig.computeLBD(c);
      if (lbd < c.lbd())
      {
         if (c.lbd() <= 30)
            c.removable(false);  // Protect once from reduction.
         c.set_lbd(lbd);
         if (lbd <= core_lbd_cut)
         {
            learnts_core.push(ref);
            c.mark(CORE);
         } else if (lbd <= 6 && c.mark() == LOCAL)
         {
            // Bug: 'cr' may already be in 'learnts_tier2', e.g., if 'cr' was demoted from TIER2
            // to LOCAL previously and if that 'cr' is not cleaned from 'learnts_tier2' yet.
            learnts_tier2.push(ref);
            c.mark(TIER2);
         }
      }

      if (c.mark() == TIER2)
         c.touched() = stat.conflicts;
      else if (c.mark() == LOCAL)
         claBumpActivity(c);
   }
}

template <typename Database>
inline void ChanseokOhReduce<Database>::claDecayActivity()
{
   cla_inc *= (1 / clause_decay);
}

template <typename Database>
inline void ChanseokOhReduce<Database>::clauseImproved(const CRef ref)
{
   Clause & c = ca[ref];
   int const lbd = c.lbd();

   assert((c.mark() == CORE || c.mark() == TIER2) && "uncomment the else if");

   if (c.mark() != CORE && lbd <= core_lbd_cut)
   {
      learnts_core.push(ref);
      c.mark(CORE);
   }
//   else if (c.mark() != TIER2 && lbd <= 6)
//   {
//      learnts_tier2.push(ref);
//      c.mark(TIER2);
//      c.touched() = stat.conflicts;
//   }
}

template <typename Database>
inline void ChanseokOhReduce<Database>::claBumpActivity(Clause& c)
{
   if ((c.activity() += cla_inc) > 1e20)
   {
      // Rescale:
      for (int i = 0; i < learnts_local.size(); i++)
         ca[learnts_local[i]].activity() *= 1e-20;
      cla_inc *= 1e-20;
   }
}

template <typename Database>
inline int ChanseokOhReduce<Database>::nClauses() const
{
   return learnts_core.size() + learnts_tier2.size() + learnts_local.size();
}

template <typename Database>
inline void ChanseokOhReduce<Database>::relocAll(Database& to)
{

   // All learnt:
   //
   for (int i = 0; i < learnts_core.size(); i++)
      ca.reloc(learnts_core[i], to);
   for (int i = 0; i < learnts_tier2.size(); i++)
      ca.reloc(learnts_tier2[i], to);
   for (int i = 0; i < learnts_local.size(); i++)
      ca.reloc(learnts_local[i], to);
}

}

#endif /* SOURCES_REDUCE_CHANSEOKOH_H_ */
