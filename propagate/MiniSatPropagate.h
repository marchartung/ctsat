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
#ifndef SOURCES_PROPAGATE_MINISATPROPAGATE_H_
#define SOURCES_PROPAGATE_MINISATPROPAGATE_H_

#include "core/Statistic.h"
#include "core/ImplicationGraph.h"
#include "mtl/Vec.h"
#include "mtl/OccLists.h"

namespace ctsat
{
template <typename DatabaseType>
class MinisatPropagate
{
   typedef typename DatabaseType::Lit Lit;
   typedef typename DatabaseType::Var Var;
   typedef typename DatabaseType::Clause Clause;
   typedef typename DatabaseType::CRef CRef;
   typedef typename DatabaseType::lbool lbool;
 public:
   typedef DatabaseType Database;

   MinisatPropagate(Statistic & stat, DatabaseType & db, ImplicationGraph<DatabaseType> & ig);

   Var newVar();
   void removeVar(Var const v);
   void relocAll(DatabaseType& to);

   void clear();

   void attachClause(CRef const cr);   // Attach a clause to watcher lists.
   int safeAttachClause(CRef const cr);  // returns either ig.nVars() when there is no problem, or a level lower than this, where a variable conflicts or propagates
   int attachLevel(CRef const cr) const;

   void detachClause(CRef const cr, bool const strict = false);  // Detach a clause to watcher lists.
   void swapWatched(CRef const cr, int const from, int const to);

   void simpleUncheckEnqueue(Lit p, CRef from = DatabaseType::npos());

   template <typename BranchType>
   void uncheckedEnqueue(BranchType & branch, Lit p, int level = 0, CRef from = DatabaseType::npos());
   // Test if fact 'p' contradicts current state, enqueue otherwise.
   template <typename BranchType>
   bool enqueue(BranchType & branch, Lit const p, CRef const from = DatabaseType::npos());

   template <typename BranchType>
   CRef propagate(BranchType & branch);
   CRef simplePropagate();

   void cancelUntilTrailRecord();

   template <typename BranchType>
   void cancelUntil(BranchType & branch, int const bLevel);

   bool isAttached(CRef const & ref) const;
   bool isBadAttached(CRef const & ref) const;

   bool binResMinimize(vec<Lit>& out_learnt);  // Further learnt clause minimization by binary resolution.
   bool extendedBinResMinimize(vec<Lit>& out_learnt);

   int qhead;  // Head of queue (as index into the trail -- no more explicit propagation queue in MiniSat).
   int trailRecord;
 private:
   Statistic & stat;
   DatabaseType & ca;
   ImplicationGraph<DatabaseType> & ig;
   vec<Lit> lowLevelLits;
   struct Watcher
   {
      CRef cref;
      Lit blocker;
      Watcher(CRef cr, Lit p)
            : cref(cr),
              blocker(p)
      {
      }
      bool operator==(const Watcher& w) const
      {
         return cref == w.cref;
      }
      bool operator!=(const Watcher& w) const
      {
         return cref != w.cref;
      }
   };

   struct WatcherDeleted
   {
      DatabaseType const & ca;
      WatcherDeleted(DatabaseType const& _ca)
            : ca(_ca)
      {
      }
      bool operator()(const Watcher& w) const
      {
         return ca[w.cref].mark() == 1;
      }
   };
   OccLists<Lit, vec<Watcher>, WatcherDeleted> watches_bin;   // Watches for binary clauses only.
   OccLists<Lit, vec<Watcher>, WatcherDeleted> watches;  // 'watches[lit]' is a list of constraints watching 'lit' (will go there if literal becomes true).

   bool unsee2BinaryImplied(Lit const & l);
};

template <typename DatabaseType>

template <typename BranchType>
inline typename MinisatPropagate<DatabaseType>::CRef MinisatPropagate<DatabaseType>::propagate(BranchType & branch)
{
   CRef confl = DatabaseType::npos();
   int num_props = 0;
   watches.cleanAll();
   watches_bin.cleanAll();

   while (qhead < ig.nAssigns())
   {
      Lit p = ig.getTrailLit(qhead++);   // 'p' is enqueued fact to propagate.
      int currLevel = ig.level(p.var());
      vec<Watcher> &ws = watches[p];
      Watcher *i, *j, *end;
      num_props++;

      vec<Watcher> &ws_bin = watches_bin[p];  // Propagate binary clauses first.
      for (int k = 0; k < ws_bin.size(); k++)
      {
         Lit the_other = ws_bin[k].blocker;
         if (ig.value(the_other).isFalse())
         {
            confl = ws_bin[k].cref;
            goto ExitProp;
         } else if (ig.value(the_other).isUndef())
         {
            uncheckedEnqueue(branch, the_other, currLevel, ws_bin[k].cref);
         }
      }

      for (i = j = (Watcher*) ws, end = i + ws.size(); i != end;)
      {
         // Try to avoid inspecting the clause:
         Lit blocker = i->blocker;
         if (ig.value(blocker).isTrue())
         {
            *j++ = *i++;
            continue;
         }

         // Make sure the false literal is data[1]:
         CRef const cr = i->cref;

         Clause& c = ca[cr];
         Lit false_lit = ~p;
         if (c[0] == false_lit)
            c[0] = c[1], c[1] = false_lit;
         assert(c[1] == false_lit);
         i++;

         // If 0th watch is true, then clause is already satisfied.
         Lit first = c[0];
         Watcher w = Watcher(cr, first);
         if (first != blocker && ig.value(first).isTrue())
         {
            *j++ = w;
            continue;
         }

         // Look for new watch:
         for (int k = 2; k < c.size(); k++)
            if (!ig.value(c[k]).isFalse())
            {
               c[1] = c[k];
               c[k] = false_lit;
               watches[~c[1]].push(w);
               goto NextClause;
            }

         // Did not find watch -- clause is unit under assignment:
         *j++ = w;
         if (ig.value(first).isFalse())
         {
            confl = cr;
            qhead = ig.nAssigns();
            // Copy the remaining watches:
            while (i < end)
               *j++ = *i++;
         } else
         {
            if (currLevel == ig.decisionLevel())
            {
               uncheckedEnqueue(branch, first, currLevel, cr);
            } else
            {
               int nMaxLevel = currLevel;
               int nMaxInd = 1;
               // pass over all the literals in the clause and find the one with the biggest level
               for (int nInd = 2; nInd < c.size(); ++nInd)
               {
                  int nLevel = ig.level(c[nInd].var());
                  if (nLevel > nMaxLevel)
                  {
                     nMaxLevel = nLevel;
                     nMaxInd = nInd;
                  }
               }

               if (nMaxInd != 1)
               {
                  std::swap(c[1], c[nMaxInd]);
                  --j;  // undo last watch
                  watches[~c[1]].push(w);
               }

               uncheckedEnqueue(branch, first, nMaxLevel, cr);
            }
         }

         NextClause: ;
      }
      ws.shrink(i - j);
   }

   ExitProp: ;
   stat.propagations += num_props;
   stat.simpDB_props -= num_props;

   return confl;
}

template <typename DatabaseType>
template <typename BranchType>
void MinisatPropagate<DatabaseType>::cancelUntil(BranchType & branch, int const bLevel)
{
   if (ig.decisionLevel() > bLevel)
   {
      lowLevelLits.clear();
      int const tEnd = ig.levelEnd(bLevel);
      for (int c = ig.nAssigns() - 1; c >= tEnd; c--)
      {
         Lit const l = ig.getTrailLit(c);
         Var const x = (l.var());

         if (ig.level(x) <= bLevel)
            lowLevelLits.push(l);
         else
         {
            branch.notifyVarUnassigned(x);
            ig.unassign(x);
         }
      }
      qhead = tEnd;
      ig.backtrack(bLevel);
      for (int nLitId = lowLevelLits.size() - 1; nLitId >= 0; --nLitId)
         ig.assign(lowLevelLits[nLitId]);

      lowLevelLits.clear();
   }
}

template <typename DatabaseType>
MinisatPropagate<DatabaseType>::MinisatPropagate(
                                             Statistic & stat,
                                             DatabaseType & db,
                                             ImplicationGraph<DatabaseType> & ig)
      : qhead(0),
        trailRecord(0),
        stat(stat),
        ca(db),
        ig(ig),
        watches_bin(WatcherDeleted(ca)),
        watches(WatcherDeleted(ca))
{

}

template <typename DatabaseType>
inline typename MinisatPropagate<DatabaseType>::CRef MinisatPropagate<DatabaseType>::simplePropagate()
{
   CRef confl = DatabaseType::npos();
   int num_props = 0;
   watches.cleanAll();
   watches_bin.cleanAll();
   while (qhead < ig.nAssigns())
   {
      Lit p = ig.getTrailLit(qhead++);   // 'p' is enqueued fact to propagate.
      vec<Watcher> &ws = watches[p];
      Watcher *i, *j, *end;
      num_props++;

      // First, Propagate binary clauses
      vec<Watcher> &wbin = watches_bin[p];

      for (int k = 0; k < wbin.size(); k++)
      {

         Lit imp = wbin[k].blocker;

         if (ig.value(imp).isFalse())
         {
            return wbin[k].cref;
         }

         if (ig.value(imp).isUndef())
         {
            simpleUncheckEnqueue(imp, wbin[k].cref);
         }
      }
      for (i = j = (Watcher*) ws, end = i + ws.size(); i != end;)
      {
         // Try to avoid inspecting the clause:
         Lit blocker = i->blocker;
         if (ig.value(blocker).isTrue())
         {
            *j++ = *i++;
            continue;
         }

         // Make sure the false literal is data[1]:
         CRef cr = i->cref;
         Clause& c = ca[cr];
         Lit false_lit = ~p;
         if (c[0] == false_lit)
            c[0] = c[1], c[1] = false_lit;
         assert(c[1] == false_lit);
         //  i++;

         // If 0th watch is true, then clause is already satisfied.
         // However, 0th watch is not the blocker, make it blocker using a new watcher w
         // why not simply do i->blocker=first in this case?
         Lit first = c[0];
         //  Watcher w     = Watcher(cr, first);
         if (first != blocker && ig.value(first).isTrue())
         {
            i->blocker = first;
            *j++ = *i++;
            continue;
         }

         else
         {
            for (int k = 2; k < c.size(); k++)
            {

               if (!ig.value(c[k]).isFalse())
               {
                  // watcher i is abandonned using i++, because cr watches now ~c[k] instead of p
                  // the blocker is first in the watcher. However,
                  // the blocker in the corresponding watcher in ~first is not c[1]
                  Watcher w = Watcher(cr, first);
                  i++;
                  c[1] = c[k];
                  c[k] = false_lit;
                  watches[~c[1]].push(w);
                  goto NextClause;
               }
            }
         }

         // Did not find watch -- clause is unit under assignment:
         i->blocker = first;
         *j++ = *i++;
         if (ig.value(first).isFalse())
         {
            confl = cr;
            qhead = ig.nAssigns();
            // Copy the remaining watches:
            while (i < end)
               *j++ = *i++;
         } else
         {
            simpleUncheckEnqueue(first, cr);
         }
         NextClause: ;
      }
      ws.shrink(i - j);
   }

   stat.s_propagations += num_props;

   return confl;
}

template <typename DatabaseType>
inline void MinisatPropagate<DatabaseType>::simpleUncheckEnqueue(Lit p, CRef from)
{
   assert(ig.value(p).isUndef());
   ig.assign(p, from);
}


template <typename DatabaseType>
template <typename BranchType>
inline void MinisatPropagate<DatabaseType>::uncheckedEnqueue(BranchType & branch, Lit p, int level, CRef from)
{

   assert(ig.value(p).isUndef());
   branch.notifyVarAssigned(p.var());
   ig.assign(p, from, level);
   assert(ig.value(p).isTrue());
}

// Test if fact 'p' contradicts current state, enqueue otherwise.
template <typename DatabaseType>
template <typename BranchType>
inline bool MinisatPropagate<DatabaseType>::enqueue(BranchType & branch, Lit const p, CRef const from)
{
   if (ig.value(p).isUndef())
   {
      uncheckedEnqueue(branch, p, ig.decisionLevel(), from);
      return true;
   } else
      return ig.value(p).isTrue();
}

template <typename DatabaseType>
inline void MinisatPropagate<DatabaseType>::cancelUntilTrailRecord()
{
   for (int c = ig.nAssigns() - 1; c >= trailRecord; c--)
   {
      Var x = (ig.getTrailLit(c).var());
      ig.unassign(x);

   }
   qhead = trailRecord;
   ig.shrink(ig.nAssigns() - trailRecord);
}

// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
template <typename DatabaseType>
typename MinisatPropagate<DatabaseType>::Var MinisatPropagate<DatabaseType>::newVar()
{
   int v = ig.nVars();
   watches_bin.init(Lit(v, false));
   watches_bin.init(Lit(v, true));
   watches.init(Lit(v, false));
   watches.init(Lit(v, true));

   return v;
}

template <typename DatabaseType>
void MinisatPropagate<DatabaseType>::removeVar(Var const v)
{
   watches_bin[Lit(v, false)].clear(true);
   watches_bin[Lit(v, true)].clear(true);
   watches[Lit(v, false)].clear(true);
   watches[Lit(v, true)].clear(true);
}
template <typename DatabaseType>
inline void MinisatPropagate<DatabaseType>::attachClause(CRef const cr)
{
   const Clause& c = ca[cr];
   assert(c.size() > 1);
   OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = c.size() == 2 ? watches_bin : watches;
   ws[~c[0]].push(Watcher(cr, c[1]));
   ws[~c[1]].push(Watcher(cr, c[0]));

   if (c.getLearnt() != 0)
   {
      stat.learnts_literals += c.size();
      ++stat.nWatchedLearnts;
   } else
   {
      stat.clauses_literals += c.size();
      ++stat.nWatchedClauses;
   }

}

template <typename DatabaseType>
int MinisatPropagate<DatabaseType>::attachLevel(CRef const cr) const
{
   Clause const& c = ca[cr];
   int undefCount = 0;
   for (int i = 0; i < c.size(); ++i)
   {
      lbool const val = ig.value(c[i]);
      undefCount += val.isUndef() + 2 * val.isTrue();
      if (undefCount > 1)
         return ig.nVars();
   }

   int highestLevel = 0, secondHighestLevel = 0;
   for (int i = 0; i < c.size(); ++i)
   {
      if (undefCount == 1 && ig.value(c[i]).isUndef())
         continue;
      int const level = ig.level(c[i].var());
      bool const isHighest = highestLevel < level;
      secondHighestLevel =
            (isHighest) ?
                  highestLevel : ((secondHighestLevel < level) ? level : secondHighestLevel);
      highestLevel = (isHighest) ? level : highestLevel;
   }
   return std::max(((undefCount == 0) ? secondHighestLevel : highestLevel) - 1, 0);

//   int undefCount = 0, highestLevel = 0, secondHighestLevel = 0;
//   for (int i = 0; i < c.size(); ++i)
//   {
//
//      lbool const val = ig.value(c[i]);
//      int const level = ig.level(c[i].var());
//      undefCount += val.isUndef() + 2*val.isTrue();
//      bool const isHighest = highestLevel < level;
//      secondHighestLevel =
//            (isHighest) ?
//                  highestLevel : ((secondHighestLevel < level) ? level : secondHighestLevel);
//      highestLevel = (isHighest) ? level : highestLevel;
//      if (undefCount > 1)
//         return ig.nVars();
//   }
//   return secondHighestLevel;
}

template <typename DatabaseType>
int MinisatPropagate<DatabaseType>::safeAttachClause(CRef const cr)
{
   Clause& c = ca[cr];
   int i = 1;
   while (i < c.size() && (ig.value(c[0]).isFalse() || ig.value(c[1]).isFalse()))
   {
      if (!ig.value(c[i]).isFalse())
      {
         if (ig.value(c[0]).isFalse())
            std::swap(c[0], c[i]);
         else if (ig.value(c[1]).isFalse())
            std::swap(c[1], c[i]);
      }
      ++i;
   }
   int res = ig.nVars();
   if (ig.value(c[0]).isFalse() || ig.value(c[1]).isFalse())
   {
      for (int j = 0; j < 2; ++j)
      {
         if (ig.value(c[j]).isFalse())
            for (int k = j + 1; k < c.size(); ++k)
               if (ig.level(c[j].var()) < ig.level(c[k].var()))
               {
                  std::swap(c[j], c[k]);
                  if (!ig.value(c[j]).isFalse())
                     break;
               }
      }
      res = ig.level(c[1].var());
   }

   attachClause(cr);
   return res;
}

template <typename DatabaseType>
bool MinisatPropagate<DatabaseType>::isAttached(CRef const & ref) const
{

   Clause const & c = ca[ref];
   OccLists<Lit, vec<Watcher>, WatcherDeleted> const & ws = c.size() == 2 ? watches_bin : watches;

   for (int k = 0; k < 2; ++k)
   {
      vec<Watcher> const & wats = ws[~c[k]];
      Watcher const lookedFor(ref, c[(k == 0) ? 1 : 0]);
      int i = 0;
      for (; i < wats.size(); ++i)
         if (wats[i] == lookedFor)
            break;
      if (i == wats.size())
         return false;
   }
   return true;
}

template <typename DatabaseType>
bool MinisatPropagate<DatabaseType>::isBadAttached(CRef const & ref) const
{
   Clause const & c = ca[ref];
   Watcher const lookedFor(ref,c[0]);
   assert(lookedFor == Watcher(ref,c[1]));
   OccLists<Lit, vec<Watcher>, WatcherDeleted> const & ws = c.size() == 2 ? watches_bin : watches;
   unsigned count = 0;
   for (int k = 0; k < c.size(); ++k)
   {
      vec<Watcher> const & wats = ws[~c[k]];
      int i = 0;
      for (; i < wats.size(); ++i)
         if (wats[i] == lookedFor)
         {
            ++count;
            if(count > 2 || i > 1)
               return true;
         }
   }
   return false;
}

template <typename DatabaseType>
inline void MinisatPropagate<DatabaseType>::detachClause(CRef const cr, bool const strict)
{
   const Clause& c = ca[cr];
   assert(c.size() > 1);
   OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = c.size() == 2 ? watches_bin : watches;

   if (strict)
   {
      for (int k = 0; k < 2; ++k)
      {
         vec<Watcher> & wats = ws[~c[k]];
         Watcher const lookFor(cr, c[(k == 0) ? 1 : 0]);
         remove(wats, lookFor);
      }
   } else
   {
      // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
      ws.smudge(~c[0]);
      ws.smudge(~c[1]);
   }

   if (c.getLearnt() != 0)
   {
      stat.learnts_literals -= c.size();
      --stat.nWatchedLearnts;
   } else
   {
      stat.clauses_literals -= c.size();
      --stat.nWatchedClauses;
   }
}

template <typename DatabaseType>
inline void MinisatPropagate<DatabaseType>::swapWatched(const CRef cr, const int from, const int to)
{
   assert(from < 2);
   Clause & conflCls = ca[cr];
   std::swap(conflCls[from], conflCls[to]);
   if (to > 1)
   {
      OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws =
            conflCls.size() == 2 ? watches_bin : watches;

      remove(ws[~conflCls[to]], Watcher(cr, conflCls[(from == 0) ? 1 : 0]));
      ws[~conflCls[from]].push(Watcher(cr, conflCls[(from == 0) ? 1 : 0]));
   }
}

template <typename DatabaseType>
inline bool MinisatPropagate<DatabaseType>::extendedBinResMinimize(vec<Lit>& c)
{
   ig.markSeen2(c, 1);
   bool removedSome = false;
   for (int i = 0; i < c.size() - 1; ++i)
      removedSome |= unsee2BinaryImplied(c[i]);
   if (removedSome)
   {
      if (!ig.removeNotSeen2(c, 1))
         assert(false);  // should have removed something
      return true;
   } else
      return false;
}

template <typename DatabaseType>
bool MinisatPropagate<DatabaseType>::unsee2BinaryImplied(Lit const & l)
{
   bool res = false;

   const vec<Watcher>& ws = watches_bin[~l];
   for (int i = 0; i < ws.size(); i++)
   {
      Lit const the_other = ws[i].blocker;
      // Does 'the_other' appear negatively in 'out_learnt'?
      if (ig.isSeen2(the_other.var()) && ig.value(the_other).isTrue())
      {
         res = true;
         ig.unsetSeen2(the_other.var());
      }
   }
   return res;
}

// Try further learnt clause minimization by means of binary clause resolution.
template <typename DatabaseType>
inline bool MinisatPropagate<DatabaseType>::binResMinimize(vec<Lit>& out_learnt)
{
   // Preparation: remember which false variables we have in 'out_learnt'.
   ig.markSeen2(out_learnt, 1);
   if (unsee2BinaryImplied(out_learnt[0]))
   {
      if (!ig.removeNotSeen2(out_learnt, 1))
         assert(false);  // should have removed something
      return true;
   } else
      return false;
}

template <typename DatabaseType>
void MinisatPropagate<DatabaseType>::clear()
{
   watches.clear(true);
   watches_bin.clear(true);
   lowLevelLits.clear(true);
}

template <typename DatabaseType>
void MinisatPropagate<DatabaseType>::relocAll(DatabaseType& to)
{
   // All watchers:
   //
   // for (int i = 0; i < watches.size(); i++)
   watches.cleanAll();
   watches_bin.cleanAll();
   for (int v = 0; v < ig.nVars(); v++)
      for (int s = 0; s < 2; s++)
      {
         Lit p(v, s);
         // printf(" >>> RELOCING: %s%d\n", sign(p)?"-":"", var(p)+1);
         vec<Watcher> &ws = watches[p];
         for (int j = 0; j < ws.size(); j++)
            ca.reloc(ws[j].cref, to);
         vec<Watcher> &ws_bin = watches_bin[p];
         for (int j = 0; j < ws_bin.size(); j++)
            ca.reloc(ws_bin[j].cref, to);
      }
}

}

#endif /* SOURCES_PROPAGATE_MINISATPROPAGATE_H_ */
