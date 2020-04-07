/*
 * DistLrbVsidsBranch.h
 *
 *  Created on: 03.03.2020
 *      Author: hartung
 */

#ifndef SOURCES_BRANCH_DISTLRBVSIDSBRANCH_H_
#define SOURCES_BRANCH_DISTLRBVSIDSBRANCH_H_

#include "core/ImplicationGraph.h"
#include "initial/SolverConfig.h"
#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "utils/Random.h"
#include "core/Statistic.h"
#include "core/SolveMode.h"
#include "utils/Timer.h"
#include "branch/Branch.h"

#define ANTI_EXPLORATION
namespace CTSat
{

template <typename Database>
class DistLrbVsidsBranch
{
 public:

   typedef typename Database::Var Var;
   typedef typename Database::Lit Lit;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
   typedef typename Database::Clause Clause;
   typedef typename Branch<Database>::BranchInputArgs BranchInputArgs;

   DistLrbVsidsBranch(typename Branch<Database>::BranchInputArgs args);

   Var newVar(bool const dvar = false);

   void insertVarOrder(Var const x);      // Insert a variable in the decision order priority queue.
   void rebuildOrderHeap();

   void setDecisionVar(Var const v, bool const b);
   bool isDecisionVar(Var const v) const;
   uint64_t nDecVars() const;

   void notifyConflictFound1(CRef const confl);
   void notifyConflictFound2(CRef const confl);
   void notifyVarSeenInConflict(Var const v);
   void notifyCreatedLearntClause(vec<Lit> const & out_learnt);
   void notifyConflictResolved();
   void notifyRestart();

   void notifyVarUnassigned(Var const v);
   void notifyVarAssigned(Var const v);

   Lit pickBranchLit();

   void varBumpActivity(Var const v, double const mult);
   void varDecayActivity();

 private:

   bool const initVarPolZero;
   bool const initRndActivity;
   bool const initRndPolarity;
   bool switched;
   bool DISTANCE;
   bool VSIDS;
   int vsids_var_decay_timer;
   double step_size;
   double step_size_dec;
   double min_step_size;

   double vsids_var_decay;
   double vsids_var_inc;          // Amount to bump next variable with.

   double dist_var_decay;
   double var_iLevel_inc;

   uint64_t dec_vars;

   vec<char> polarity;         // The preferred polarity of each variable.
   vec<char> decision;  // Declares if a variable is eligible for selection in the decision heuristic

   vec<double> activity_CHB;     // A heuristic measurement of the activity of a variable.
   vec<double> activity_VSIDS;
   vec<double> activity_distance;

   struct VarOrderLt
   {
      const vec<double>& activity;
      bool operator ()(Var x, Var y) const
      {
         return activity[x] > activity[y];
      }
      VarOrderLt(const vec<double>& act)
            : activity(act)
      {
      }
   };

   Heap<VarOrderLt> order_heap_CHB;  // A priority queue of variables ordered with respect to the variable activity.
   Heap<VarOrderLt> order_heap_VSIDS;
   Heap<VarOrderLt> order_heap_distance;

   vec<uint32_t> picked;
   vec<uint32_t> conflicted;
   vec<uint32_t> almost_conflicted;
#ifdef ANTI_EXPLORATION
   vec<uint32_t> canceled;
#endif

   vec<Lit> involved_lits;
   vec<double> var_iLevel;
   vec<double> var_iLevel_tmp;
   vec<int> pathCs;

   vec<Lit> add_tmp;
   vec<Var> bumpLater;
   Timer switchTimer;

   SolveMode & smode;
   Random & rand;
   Statistic & stat;
   Database & ca;
   ImplicationGraph<Database> & ig;

   bool collectFirstUIP(CRef confl);

   bool isVsids() const;
   void setVsids(bool const b);

   bool isDistance() const;
   void setDistance(bool const b);

};

template <typename Database>
inline void DistLrbVsidsBranch<Database>::notifyRestart()
{
   if (!switched && stat.conflicts > 10000)
   {
      switched = true;
      setVsids(false);
   } else if (!isVsids() && switchTimer.isOver())
   {
      setVsids(true);
      printf("c Switched to VSIDS.\n");
      fflush(stdout);
   }
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::setVsids(bool const b)
{
   VSIDS = b;
   smode.branchingLrb = !b;
}

template <typename Database>
inline bool DistLrbVsidsBranch<Database>::isDistance() const
{
   return DISTANCE;
}
template <typename Database>
inline void DistLrbVsidsBranch<Database>::setDistance(bool const b)
{
   DISTANCE = b;
}

template <typename Database>
inline bool DistLrbVsidsBranch<Database>::isVsids() const
{
   return VSIDS;
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::rebuildOrderHeap()
{
   vec<Var> vs;
   for (Var v = 0; v < ig.nVars(); v++)
      if (decision[v] && ig.value(v) == lbool::Undef())
         vs.push(v);

   order_heap_CHB.build(vs);
   order_heap_VSIDS.build(vs);
   order_heap_distance.build(vs);
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::varDecayActivity()
{
   vsids_var_inc *= (1 / vsids_var_decay);
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::varBumpActivity(Var const v, double const mult)
{
   if ((activity_VSIDS[v] += vsids_var_inc * mult) > 1e100)
   {
      // Rescale:
      for (int i = 0; i < ig.nVars(); i++)
         activity_VSIDS[i] *= 1e-100;
      vsids_var_inc *= 1e-100;
   }

   // Update order_heap with respect to new activity:
   if (order_heap_VSIDS.inHeap(v))
      order_heap_VSIDS.decrease(v);
}

template <typename Database>
inline bool DistLrbVsidsBranch<Database>::isDecisionVar(Var const v) const
{
   return decision[v];
}

template <typename Database>
inline uint64_t DistLrbVsidsBranch<Database>::nDecVars() const
{
   return dec_vars;
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::setDecisionVar(Var const v, bool const b)
{
   if (b && !decision[v])
      dec_vars++;
   else if (!b && decision[v])
      dec_vars--;

   decision[v] = b;
   if (b && !order_heap_CHB.inHeap(v))
   {
      order_heap_CHB.insert(v);
      order_heap_VSIDS.insert(v);
      order_heap_distance.insert(v);
   }
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::insertVarOrder(Var const x)
{
   //    Heap<VarOrderLt>& order_heap = VSIDS ? order_heap_VSIDS : order_heap_CHB;
   Heap<VarOrderLt>& order_heap =
         DISTANCE ? order_heap_distance : ((!VSIDS) ? order_heap_CHB : order_heap_VSIDS);
   if (!order_heap.inHeap(x) && decision[x])
      order_heap.insert(x);
}

template <typename Database>
inline typename DistLrbVsidsBranch<Database>::Var DistLrbVsidsBranch<Database>::newVar(bool const dvar)
{
   Var const v = activity_CHB.size();
   double const activityVal = (initRndActivity) ? rand.drand() * 0.0001 : 0.0;
   bool const polarityVal = (initRndPolarity) ? rand.irand(2) : !initVarPolZero;
   activity_CHB.push(activityVal);
   activity_VSIDS.push(activityVal);
   activity_distance.push(activityVal);
   var_iLevel.push(0);
   var_iLevel_tmp.push(0);
   pathCs.push(0);
   polarity.push(polarityVal);
   decision.push(false);
   picked.push(0);
   conflicted.push(0);
   almost_conflicted.push(0);
#ifdef ANTI_EXPLORATION
   canceled.push(0);
#endif
   if (dvar)
      setDecisionVar(v, true);
   return v;
}

// pathCs[k] is the number of variables assigned at level k,
// it is initialized to 0 at the begining and reset to 0 after the function execution
template <typename Database>
inline bool DistLrbVsidsBranch<Database>::collectFirstUIP(CRef confl)
{
   involved_lits.clear();
   int max_level = 1;
   Clause& c = ca[confl];
   int minLevel = ig.decisionLevel();
   for (int i = 0; i < c.size(); i++)
   {
      Var v = c[i].var();
      //        assert(!seen[v]);
      if (ig.level(v) > 0)
      {
         ig.setSeen(v);
         var_iLevel_tmp[v] = 1;
         pathCs[ig.level(v)]++;
         if (minLevel > ig.level(v))
         {
            minLevel = ig.level(v);
            assert(minLevel > 0);
         }
         //    varBumpActivity(v);
      }
   }

   int limit = ig.levelEnd(minLevel - 1);
   for (int i = ig.nAssigns() - 1; i >= limit; i--)
   {
      Lit const p = ig.getTrailLit(i);
      Var const v = p.var();
      if (ig.isSeen(v))
      {
         int const currentDecLevel = ig.level(v);
         //      if (currentDecLevel==decisionLevel())
         //       varBumpActivity(v);
         ig.unsetSeen(v);
         if (--pathCs[currentDecLevel] != 0)
         {
            Clause& rc = ca[ig.reason(v)];
            int reasonVarLevel = var_iLevel_tmp[v] + 1;
            if (reasonVarLevel > max_level)
               max_level = reasonVarLevel;
            if (rc.size() == 2 && ig.value(rc[0]) == lbool::False())
            {
               // Special case for binary clauses
               // The first one has to be SAT
               assert(ig.value(rc[1]) != lbool::False());
               Lit tmp = rc[0];
               rc[0] = rc[1], rc[1] = tmp;
            }
            for (int j = 1; j < rc.size(); j++)
            {
               Lit q = rc[j];
               Var v1 = q.var();
               if (ig.level(v1) > 0)
               {
                  if (minLevel > ig.level(v1))
                  {
                     minLevel = ig.level(v1);
                     limit = ig.levelEnd(minLevel - 1);
                     assert(minLevel > 0);
                  }
                  if (ig.isSeen(v1))
                  {
                     if (var_iLevel_tmp[v1] < reasonVarLevel)
                        var_iLevel_tmp[v1] = reasonVarLevel;
                  } else
                  {
                     var_iLevel_tmp[v1] = reasonVarLevel;
                     //   varBumpActivity(v1);
                     ig.setSeen(v1);
                     pathCs[ig.level(v1)]++;
                  }
               }
            }
         }
         involved_lits.push(p);
      }
   }
   double inc = var_iLevel_inc;
   vec<int> level_incs;
   level_incs.clear();
   for (int i = 0; i < max_level; i++)
   {
      level_incs.push(inc);
      inc = inc / dist_var_decay;
   }

   for (int i = 0; i < involved_lits.size(); i++)
   {
      Var v = involved_lits[i].var();
      //        double old_act=activity_distance[v];
      //        activity_distance[v] +=var_iLevel_inc * var_iLevel_tmp[v];
      activity_distance[v] += var_iLevel_tmp[v] * level_incs[var_iLevel_tmp[v] - 1];

      if (activity_distance[v] > 1e100)
      {
         for (int vv = 0; vv < ig.nVars(); vv++)
            activity_distance[vv] *= 1e-100;
         var_iLevel_inc *= 1e-100;
         for (int j = 0; j < max_level; j++)
            level_incs[j] *= 1e-100;
      }
      if (order_heap_distance.inHeap(v))
         order_heap_distance.decrease(v);

      //        var_iLevel_inc *= (1 / my_var_decay);
   }
   var_iLevel_inc = level_incs[level_incs.size() - 1];
   ig.clearSeen();
   return true;
}

template <typename Database>
inline DistLrbVsidsBranch<Database>::DistLrbVsidsBranch(
                                                        typename Branch<Database>::BranchInputArgs args)
      : initVarPolZero(args.config.initVarPolZero),
        initRndActivity(args.config.rnd_active),
        initRndPolarity(args.config.rnd_polarity),
        switched(false),
        DISTANCE(true),
        VSIDS(true),
        vsids_var_decay_timer(args.config.vsids_var_decay_timer),
        step_size(args.config.step_size),
        step_size_dec(args.config.step_size_dec),
        min_step_size(args.config.min_step_size),
        vsids_var_decay(args.config.vsids_var_decay),
        vsids_var_inc(1),
        dist_var_decay(args.config.dist_var_decay),
        var_iLevel_inc(1),
        dec_vars(0),
        order_heap_CHB(VarOrderLt(activity_CHB)),
        order_heap_VSIDS(VarOrderLt(activity_VSIDS)),
        order_heap_distance(VarOrderLt(activity_distance)),
        switchTimer(args.config.timeToBranchSwitch),
        smode(args.smode),
        rand(args.rand),
        stat(args.stat),
        ca(args.ca),
        ig(args.ig)
{

   setVsids(true);
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::notifyConflictFound1(CRef const confl)
{
   if (VSIDS)
   {
      if (--vsids_var_decay_timer == 0 && vsids_var_decay < 0.95)
         vsids_var_decay_timer = 5000, vsids_var_decay += 0.01;
   } else if (step_size > min_step_size)
      step_size -= step_size_dec;
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::notifyConflictFound2(CRef const confl)
{
   DISTANCE = stat.conflicts <= 50000;
   if (VSIDS && DISTANCE)
      collectFirstUIP(confl);
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::notifyVarSeenInConflict(const Var v)
{
   if (VSIDS)
   {
      varBumpActivity(v, .5);
      bumpLater.push(v);
   } else
      ++conflicted[v];
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::notifyCreatedLearntClause(const vec<Lit>& out_learnt)
{
   if (VSIDS)
   {
      int const btlevel = (out_learnt.size() > 1) ? ig.level(out_learnt[1].var()) - 1 : 0;
      for (int i = 0; i < add_tmp.size(); i++)
      {
         Var v = (add_tmp[i].var());
         if (ig.level(v) >= btlevel)
            varBumpActivity(v, 1);
      }
      add_tmp.clear();
   } else
   {
      ig.setSeen(out_learnt[0].var());
      ig.markSeenToClear(out_learnt[0]);
      for (int i = out_learnt.size() - 1; i >= 0; --i)
      {
         Var v = (out_learnt[i].var());
         CRef rea = ig.reason(v);
         if (rea != CRef_Undef)
         {
            const Clause& reaC = ca[rea];
            for (int i = 0; i < reaC.size(); ++i)
            {
               Lit l = reaC[i];
               if (!ig.isSeen(l.var()))
               {
                  ig.setSeen(l.var());
                  ++almost_conflicted[(l.var())];
                  ig.markSeenToClear(l);
               }
            }
         }
      }
   }
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::notifyConflictResolved()
{
   if (VSIDS)
      varDecayActivity();
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::notifyVarUnassigned(const Var x)
{
   if (!VSIDS)
   {
      uint32_t const age = stat.conflicts - picked[x];
      if (age > 0)
      {
         double const adjusted_reward = ((double) (conflicted[x] + almost_conflicted[x]))
            / ((double) age);
         double const old_activity = activity_CHB[x];
         activity_CHB[x] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
         if (order_heap_CHB.inHeap(x))
         {
            if (activity_CHB[x] > old_activity)
               order_heap_CHB.decrease(x);
            else
               order_heap_CHB.increase(x);
         }
      }
#ifdef ANTI_EXPLORATION
      canceled[x] = stat.conflicts;
#endif
   }

   polarity[x] = ig[x] == lbool::False();
   insertVarOrder(x);
}

template <typename Database>
inline void DistLrbVsidsBranch<Database>::notifyVarAssigned(const Var x)
{
   if (!VSIDS)
   {
      picked[x] = stat.conflicts;
      conflicted[x] = 0;
      almost_conflicted[x] = 0;
#ifdef ANTI_EXPLORATION
      uint32_t age = stat.conflicts - canceled[x];
      if (age > 0)
      {
         double const decay = pow(0.95, age);
         activity_CHB[x] *= decay;
         if (order_heap_CHB.inHeap(x))
            order_heap_CHB.increase(x);
      }
#endif
   }
}

template <typename Database>
inline typename DistLrbVsidsBranch<Database>::Lit DistLrbVsidsBranch<Database>::pickBranchLit()
{
   Var next = var_Undef;
   Heap<VarOrderLt>& order_heap =
         DISTANCE ? order_heap_distance : ((!VSIDS) ? order_heap_CHB : order_heap_VSIDS);

   // Activity based decision:
   while (next == var_Undef || ig.value(next) != lbool::Undef() || !decision[next])
      if (order_heap.empty())
         return Lit::Undef();
      else
      {
#ifdef ANTI_EXPLORATION
         if (!VSIDS)
         {
            Var v = order_heap_CHB[0];
            uint32_t age = stat.conflicts - canceled[v];
            while (age > 0)
            {
               double const decay = pow(0.95, age);
               activity_CHB[v] *= decay;
               if (order_heap_CHB.inHeap(v))
                  order_heap_CHB.increase(v);
               canceled[v] = stat.conflicts;
               v = order_heap_CHB[0];
               age = stat.conflicts - canceled[v];
            }
         }
#endif
         next = order_heap.removeMin();
      }

   return Lit(next, polarity[next]);
}
} /* namespace Minisat */

#endif /* SOURCES_BRANCH_DISTLRBVSIDSBRANCH_H_ */
