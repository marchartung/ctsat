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

#ifndef SOURCES_REDUCE_GlucoseRestart_H_
#define SOURCES_REDUCE_GlucoseRestart_H_

#include "core/Statistic.h"
#include "initial/SolverConfig.h"
#include "core/ImplicationGraph.h"

namespace ctsat
{

template <typename Database>
class GlucoseReduce
{
   typedef typename Database::Lit Lit;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;

 public:
   GlucoseReduce(
                 SolverConfig const & config,
                 Statistic & stat,
                 Database & db,
                 ImplicationGraph<Database> & ig);

   // returns true, if clauses were deleted
   template <typename RemoveFunctor>
   inline bool run(RemoveFunctor const & remove)
   {
      bool res = false;
      if (stat.conflicts >= nRestarts * nbclausesbeforereduce)
      {

         if (learnts.size() > 0)
         {
            nRestarts = (stat.conflicts / nbclausesbeforereduce) + 1;
            reduceDB(remove);
            nbclausesbeforereduce += incReduceDB;
            res = true;
         }
      }
      return res;
   }

   template <typename RemoveFunctor>
   void removeSatisfied(RemoveFunctor const & remove)
   {
      safeRemoveSatisfied(remove, learnts);  // Should clean core first.
   }

   // applies the given functor to each clause. When the functor returns true, the clause will be removed
   template <typename ApplyFunctor>
   void applyRemoveGoodClauses(ApplyFunctor const & func)
   {
      int const limit = learnts.size() / 2;
      int i = limit, j = limit;
      sort(learnts, reduceDB_lt(ca));
      for (; i < learnts.size(); ++i)
         if (!func(learnts[i]))
            learnts[j++] = learnts[i];
      learnts.shrink(i - j);
   }

   void addClause(CRef const ref);

   void clauseImproved(CRef const ref);
   void clauseUsedInConflict(CRef const ref);
   void adjustOnConflict();

   void notifyConflictResolved();

   int nClauses() const;

   void relocAll(Database& to);

 public:

   Statistic & stat;
   Database & ca;
   ImplicationGraph<Database> & ig;

   int maxProtectableLbd;
   uint64_t nbclausesbeforereduce;
   uint64_t incReduceDB;
   uint64_t specialIncReduceDB;
   uint64_t nRestarts;
   double clause_decay;
   double cla_inc;   // Amount to bump next clause with.

   vec<CRef> learnts;
   struct reduceDB_lt
   {
      Database const & ca;
      reduceDB_lt(Database const& ca_)
            : ca(ca_)
      {
      }
      bool operator ()(CRef const & x, CRef const & y) const
      {
         Clause const & a = ca[x], &b = ca[y];
         if (a.size() == 2 || b.size() == 2)
            return a.size() > b.size();
         if (a.lbd() != b.lbd())
            return a.lbd() > b.lbd();
         return ca[x].activity() < ca[y].activity();
      }
   };
   void claDecayActivity();  // Decay all clauses with the specified factor. Implemented by increasing the 'bump' value instead.
   void claBumpActivity(Clause& c);   // Increase a clause with the current 'bump' value.

   template <typename RemoveFunctor>
   void reduceDB(RemoveFunctor const & remove);

   template <typename RemoveFunctor>
   void safeRemoveSatisfied(RemoveFunctor const & remove, vec<CRef>& cs)
   {
      int i, j;
      for (i = j = 0; i < cs.size(); i++)
      {
         Clause& c = ca[cs[i]];
         if (ig.satisfied(c))
            remove(cs[i]);
         else
            cs[j++] = cs[i];
      }
      cs.shrink(i - j);
   }

};

template <typename Database>
GlucoseReduce<Database>::GlucoseReduce(
                                       SolverConfig const & config,
                                       Statistic & stat,
                                       Database & db,
                                       ImplicationGraph<Database> & ig)
      : stat(stat),
        ca(db),
        ig(ig),
        maxProtectableLbd(config.maxProtectableLbd),
        nbclausesbeforereduce(config.firstReduceDB),
        incReduceDB(config.incReduceDB),
        specialIncReduceDB(config.specialIncReduceDB),
        nRestarts(1),
        clause_decay(config.clause_decay),
        cla_inc(1),
        learnts()
{

}

template <typename Database>
template <typename RemoveFunctor>
void GlucoseReduce<Database>::reduceDB(RemoveFunctor const & remove)   // Reduce the set of learnt clauses.
{
   int i, j;
   sort(learnts, reduceDB_lt(ca));

   // We have a lot of "good" clauses, it is difficult to compare them. Keep more !
   if (ca[learnts[learnts.size() / 2]].lbd() <= 3)
      nbclausesbeforereduce += specialIncReduceDB;
   // Useless :-)
   if (ca[learnts.last()].lbd() <= 5)
      nbclausesbeforereduce += specialIncReduceDB;

   // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
   // Keep clauses which seem to be usefull (their lbd was reduce during this sequence)

   int limit = learnts.size() / 2;

   for (i = j = 0; i < learnts.size(); i++)
   {
      Clause& c = ca[learnts[i]];
      if (c.lbd() > 2 && c.size() > 2 && c.removable() && !ig.locked(c) && (i < limit))
      {
         remove(learnts[i]);
      } else
      {
         if (!c.removable())
            limit++;  //we keep c, so we can delete an other clause
         c.removable(true);       // At the next step, c can be delete
         learnts[j++] = learnts[i];
      }
   }
   learnts.shrink(i - j);
}

template <typename Database>
void GlucoseReduce<Database>::addClause(CRef const ref)
{

   Clause & c = ca[ref];
   learnts.push(ref);
   claBumpActivity(c);
}

template <typename Database>
void GlucoseReduce<Database>::notifyConflictResolved()
{
   claDecayActivity();
}

template <typename Database>
void GlucoseReduce<Database>::adjustOnConflict()
{
}

template <typename Database>
inline void GlucoseReduce<Database>::clauseUsedInConflict(CRef const ref)
{
   Clause & c = ca[ref];
   // DYNAMIC NBLEVEL trick (see competition'09 companion paper)
   if (c.learnt() && c.lbd() > 2)
   {
      int const nblevels = ig.computeLBD(c);
      if (nblevels + 1 < c.lbd())
      {  // improve the LBD
         if (c.lbd() <= maxProtectableLbd)
            c.removable(false);

         // seems to be interesting : keep it for the next round
         c.set_lbd(nblevels);  // Update it
      }
   }
}

template <typename Database>
inline void GlucoseReduce<Database>::claDecayActivity()
{
   cla_inc *= (1 / clause_decay);
}

template <typename Database>
inline void GlucoseReduce<Database>::clauseImproved(const CRef ref)
{
}

template <typename Database>
inline void GlucoseReduce<Database>::claBumpActivity(Clause& c)
{
   if ((c.activity() += cla_inc) > 1e20)
   {
      for (int i = 0; i < learnts.size(); i++)
         ca[learnts[i]].activity() *= 1e-20;
      cla_inc *= 1e-20;
   }
}

template <typename Database>
inline int GlucoseReduce<Database>::nClauses() const
{
   return learnts.size();
}

template <typename Database>
inline void GlucoseReduce<Database>::relocAll(Database& to)
{

   // All learnt:
   //
   for (int i = 0; i < learnts.size(); i++)
      ca.reloc(learnts[i], to);
}

}

#endif /* SOURCES_REDUCE_GlucoseRestart_H_ */
