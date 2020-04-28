/*
 * VsidsBranch.h
 *
 *  Created on: 03.03.2020
 *      Author: hartung
 */

#ifndef SOURCES_BRANCH_VsidsBranch_H_
#define SOURCES_BRANCH_VsidsBranch_H_

#include "branch/Branch.h"
#include "core/ImplicationGraph.h"
#include "initial/Inputs.h"
#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "utils/Random.h"

namespace ctsat
{

template <typename Database>
class VsidsBranch : public Branch<Database>
{
   typedef Branch<Database> Super;
 public:

   typedef typename Database::Var Var;
   typedef typename Database::Lit Lit;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
   typedef typename Database::Clause Clause;
   typedef typename Branch<Database>::BranchInputArgs BranchInputArgs;

   VsidsBranch(typename Super::BranchInputArgs args);

   void notifyConflictFound1(CRef const confl);
   void notifyConflictFound2(CRef const confl);
   void notifyVarSeenInConflict(Var const v);
   void notifyCreatedLearntClause(vec<Lit> const & out_learnt);
   void notifyConflictResolved();
   void notifyRestart();

   void notifyVarAssigned(Var const v);

   Lit pickBranchLit();

 private:
   int const init_var_decay_timer;
   int var_decay_timer;
   double const max_var_decay;
   double var_decay;
   double var_inc;          // Amount to bump next variable with.

//   vec<Lit> add_tmp;
   vec<Var> bumpLater;

   void varDecayActivity();
   void varBumpActivity(Var v, double mult);

};

template <typename Database>
inline void VsidsBranch<Database>::notifyRestart()
{

}

template <typename Database>
inline void VsidsBranch<Database>::varDecayActivity()
{
   var_inc *= (1 / var_decay);
}

template <typename Database>
inline void VsidsBranch<Database>::varBumpActivity(Var v, double mult)
{
   if ((Super::activity[v] += var_inc * mult) > 1e100)
   {
      // Rescale:
      for (int i = 0; i < Super::ig.nVars(); i++)
         Super::activity[i] *= 1e-100;
      var_inc *= 1e-100;
   }

   // Update order_heap with respect to new activity:
   if (Super::order_heap.inHeap(v))
      Super::order_heap.decrease(v);
}

template <typename Database>
inline VsidsBranch<Database>::VsidsBranch(typename Super::BranchInputArgs args)
        : Super(args),
          init_var_decay_timer(args.config.vsids_var_decay_timer),
        var_decay_timer(args.config.vsids_var_decay_timer),
        max_var_decay(args.config.vsids_max_var_decay),
        var_decay(args.config.vsids_var_decay),
        var_inc(1)
{
}

template <typename Database>
inline void VsidsBranch<Database>::notifyConflictFound1(CRef const confl)
{
   if (--var_decay_timer == 0 && var_decay < max_var_decay)
   {
      var_decay_timer = init_var_decay_timer;
      var_decay = std::min(var_decay + 0.01, max_var_decay);
   }
}

template <typename Database>
inline void VsidsBranch<Database>::notifyConflictFound2(CRef const confl)
{
}

template <typename Database>
inline void VsidsBranch<Database>::notifyVarSeenInConflict(const Var v)
{
   varBumpActivity(v, .5);
   bumpLater.push(v);
}

template <typename Database>
inline void VsidsBranch<Database>::notifyCreatedLearntClause(const vec<Lit>& out_learnt)
{
   int const btlevel = (out_learnt.size() > 1) ? Super::ig.level(out_learnt[1].var()) - 1 : 0;
   for (int i = 0; i < bumpLater.size(); i++)
   {
      Var v = bumpLater[i];
      if (Super::ig.level(v) >= btlevel)
         varBumpActivity(v, 1);
   }
   bumpLater.clear();
}

template <typename Database>
inline void VsidsBranch<Database>::notifyConflictResolved()
{
   varDecayActivity();
}

template <typename Database>
inline void VsidsBranch<Database>::notifyVarAssigned(const Var x)
{

}

template <typename Database>
inline typename VsidsBranch<Database>::Lit VsidsBranch<Database>::pickBranchLit()
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

#endif /* SOURCES_BRANCH_VsidsBranch_H_ */
