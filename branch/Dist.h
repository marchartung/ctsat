/*
 * DistBranch.h
 *
 *  Created on: 03.03.2020
 *      Author: hartung
 */

#ifndef SOURCES_BRANCH_DistBranch_H_
#define SOURCES_BRANCH_DistBranch_H_

#include "branch/Branch.h"
#include "core/ImplicationGraph.h"
#include "initial/Inputs.h"
#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "utils/Random.h"

namespace ctsat
{

template <typename Database>
class DistBranch : public Branch<Database>
{
   typedef Branch<Database> Super;
 public:

   typedef typename Database::Var Var;
   typedef typename Database::Lit Lit;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
   typedef typename Database::Clause Clause;
   typedef typename Branch<Database>::BranchInputArgs BranchInputArgs;

   DistBranch(typename Super::BranchInputArgs args);

   Var newVar(bool const dvar = true);

   void notifyConflictFound1(CRef const confl);
   void notifyConflictFound2(CRef const confl);
   void notifyVarSeenInConflict(Var const v);
   void notifyCreatedLearntClause(vec<Lit> const & out_learnt);
   void notifyConflictResolved();
   void notifyRestart();

   void notifyVarAssigned(Var const v);

   Lit pickBranchLit();

 private:

   double dist_var_decay;
   double var_iLevel_inc;

   vec<Lit> involved_lits;
   vec<double> var_iLevel;
   vec<double> var_iLevel_tmp;
   vec<int> pathCs;

   bool collectFirstUIP(CRef confl);

};
template <typename Database>
inline void DistBranch<Database>::notifyRestart()
{

}

template <typename Database>
inline typename DistBranch<Database>::Var DistBranch<Database>::newVar(bool const dvar)
{
   Var const v = Super::newVar(dvar);
   var_iLevel.push(0);
   var_iLevel_tmp.push(0);
   pathCs.push(0);
   return v;
}

// pathCs[k] is the number of variables assigned at level k,
// it is initialized to 0 at the begining and reset to 0 after the function execution
template <typename Database>
bool DistBranch<Database>::collectFirstUIP(CRef confl)
{
   involved_lits.clear();
   int max_level = 1;
   Clause& c = Super::ca[confl];
   int minLevel = Super::ig.decisionLevel();
   for (int i = 0; i < c.size(); i++)
   {
      Var v = c[i].var();
      //        assert(!seen[v]);
      if (Super::ig.level(v) > 0)
      {
         Super::ig.setSeen(v);
         var_iLevel_tmp[v] = 1;
         pathCs[Super::ig.level(v)]++;
         if (minLevel > Super::ig.level(v))
         {
            minLevel = Super::ig.level(v);
            assert(minLevel > 0);
         }
         //    varBumpActivity(v);
      }
   }

   int limit = Super::ig.levelEnd(minLevel - 1);
   for (int i = Super::ig.nAssigns() - 1; i >= limit; i--)
   {
      Lit const p = Super::ig.getTrailLit(i);
      Var const v = p.var();
      if (Super::ig.isSeen(v))
      {
         int const currentDecLevel = Super::ig.level(v);
         //      if (currentDecLevel==decisionLevel())
         //       varBumpActivity(v);
         Super::ig.unsetSeen(v);
         if (--pathCs[currentDecLevel] != 0)
         {
            Clause& rc = Super::ca[Super::ig.reason(v)];
            int reasonVarLevel = var_iLevel_tmp[v] + 1;
            if (reasonVarLevel > max_level)
               max_level = reasonVarLevel;
            if (rc.size() == 2 && Super::ig.value(rc[0]) == lbool::False())
            {
               // Special case for binary clauses
               // The first one has to be SAT
               assert(Super::ig.value(rc[1]) != lbool::False());
               Lit tmp = rc[0];
               rc[0] = rc[1], rc[1] = tmp;
            }
            for (int j = 1; j < rc.size(); j++)
            {
               Lit q = rc[j];
               Var v1 = q.var();
               if (Super::ig.level(v1) > 0)
               {
                  if (minLevel > Super::ig.level(v1))
                  {
                     minLevel = Super::ig.level(v1);
                     limit = Super::ig.levelEnd(minLevel - 1);
                     assert(minLevel > 0);
                  }
                  if (Super::ig.isSeen(v1))
                  {
                     if (var_iLevel_tmp[v1] < reasonVarLevel)
                        var_iLevel_tmp[v1] = reasonVarLevel;
                  } else
                  {
                     var_iLevel_tmp[v1] = reasonVarLevel;
                     //   varBumpActivity(v1);
                     Super::ig.setSeen(v1);
                     pathCs[Super::ig.level(v1)]++;
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
      //        double old_act=Super::activity[v];
      //        Super::activity[v] +=var_iLevel_inc * var_iLevel_tmp[v];
      Super::activity[v] += var_iLevel_tmp[v] * level_incs[var_iLevel_tmp[v] - 1];

      if (Super::activity[v] > 1e100)
      {
         for (int vv = 0; vv < Super::ig.nVars(); vv++)
            Super::activity[vv] *= 1e-100;
         var_iLevel_inc *= 1e-100;
         for (int j = 0; j < max_level; j++)
            level_incs[j] *= 1e-100;
      }
      if (Super::order_heap.inHeap(v))
         Super::order_heap.decrease(v);

      //        var_iLevel_inc *= (1 / my_var_decay);
   }
   var_iLevel_inc = level_incs[level_incs.size() - 1];
   Super::ig.clearSeen();
   return true;
}

template <typename Database>
inline DistBranch<Database>::DistBranch(typename Super::BranchInputArgs args)
      : Super(args),
        dist_var_decay(args.config.dist_var_decay),
        var_iLevel_inc(1)
{
}

template <typename Database>
inline void DistBranch<Database>::notifyConflictFound1(CRef const confl)
{
}

template <typename Database>
inline void DistBranch<Database>::notifyConflictFound2(CRef const confl)
{
   collectFirstUIP(confl);
}

template <typename Database>
inline void DistBranch<Database>::notifyVarSeenInConflict(const Var v)
{
}

template <typename Database>
inline void DistBranch<Database>::notifyCreatedLearntClause(const vec<Lit>& out_learnt)
{
}

template <typename Database>
inline void DistBranch<Database>::notifyConflictResolved()
{
}

template <typename Database>
inline void DistBranch<Database>::notifyVarAssigned(const Var x)
{

}

template <typename Database>
inline typename DistBranch<Database>::Lit DistBranch<Database>::pickBranchLit()
{
   Var next = var_Undef;
   // Activity based decision:
   while (next == var_Undef || Super::ig.value(next) != lbool::Undef() || !Super::decision[next])
      if (Super::order_heap.empty())
         return Lit::Undef();
      else
         next = Super::order_heap.removeMin();

   return Lit(next, Super::polarity[next]);
}
} /* namespace Minisat */

#endif /* SOURCES_BRANCH_DistBranch_H_ */
