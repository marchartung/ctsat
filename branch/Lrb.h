/*
 * LrbBranch.h
 *
 *  Created on: 03.03.2020
 *      Author: hartung
 */

#ifndef SOURCES_BRANCH_LrbBranch_H_
#define SOURCES_BRANCH_LrbBranch_H_

#include "branch/Branch.h"
#include "core/ImplicationGraph.h"
#include "initial/Inputs.h"
#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "utils/Random.h"

#ifndef ANTI_EXPLORATION
#define ANTI_EXPLORATION
#endif

namespace ctsat
{

template <typename Database>
class LrbBranch : public Branch<Database>
{
   typedef Branch<Database> Super;
 public:

   typedef typename Database::Var Var;
   typedef typename Database::Lit Lit;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
   typedef typename Database::Clause Clause;
   typedef typename Branch<Database>::BranchInputArgs BranchInputArgs;

   LrbBranch(typename Super::BranchInputArgs args);

   Var newVar(bool const dvar = true);

   void notifyConflictFound1(CRef const confl);
   void notifyConflictFound2(CRef const confl);
   void notifyVarSeenInConflict(Var const v);
   void notifyCreatedLearntClause(vec<Lit> const & out_learnt);
   void notifyConflictResolved();
   void notifyRestart();

   void notifyVarUnassigned(Var const v);
   void notifyVarAssigned(Var const v);

   Lit pickBranchLit();

 private:
   double step_size;
   double step_size_dec;
   double min_step_size;

   vec<uint32_t> picked;
   vec<uint32_t> conflicted;
   vec<uint32_t> almost_conflicted;
#ifdef ANTI_EXPLORATION
   vec<uint32_t> canceled;
#endif

};
template <typename Database>
inline void LrbBranch<Database>::notifyRestart()
{

}

template <typename Database>
inline typename LrbBranch<Database>::Var LrbBranch<Database>::newVar(bool const dvar)
{
   Var const v = Super::newVar(dvar);

   picked.push(0);
   conflicted.push(0);
   almost_conflicted.push(0);
#ifdef ANTI_EXPLORATION
   canceled.push(0);
#endif
   return v;
}

template <typename Database>
inline LrbBranch<Database>::LrbBranch(typename Super::BranchInputArgs args)
    : Super(args),
        step_size(args.config.step_size),
        step_size_dec(args.config.step_size_dec),
        min_step_size(args.config.min_step_size)
{
}

template <typename Database>
inline void LrbBranch<Database>::notifyConflictFound1(CRef const confl)
{
   if (step_size > min_step_size)
      step_size -= step_size_dec;
}

template <typename Database>
inline void LrbBranch<Database>::notifyConflictFound2(CRef const confl)
{
}

template <typename Database>
inline void LrbBranch<Database>::notifyVarSeenInConflict(const Var v)
{
   ++conflicted[v];
}

template <typename Database>
inline void LrbBranch<Database>::notifyCreatedLearntClause(const vec<Lit>& out_learnt)
{
   Super::ig.setSeen(out_learnt[0].var());
   Super::ig.markSeenToClear(out_learnt[0]);
   for (int i = out_learnt.size() - 1; i >= 0; --i)
   {
      Var v = (out_learnt[i].var());
      CRef rea = Super::ig.reason(v);
      if (rea != CRef_Undef)
      {
         const Clause& reaC = Super::ca[rea];
         for (int i = 0; i < reaC.size(); ++i)
         {
            Lit l = reaC[i];
            if (!Super::ig.isSeen(l.var()))
            {
               Super::ig.setSeen(l.var());
               ++almost_conflicted[(l.var())];
               Super::ig.markSeenToClear(l);
            }
         }
      }
   }
}

template <typename Database>
inline void LrbBranch<Database>::notifyConflictResolved()
{
}

template <typename Database>
inline void LrbBranch<Database>::notifyVarUnassigned(const Var x)
{
   uint32_t const age = Super::stat.conflicts - picked[x];
   if (age > 0)
   {
      double adjusted_reward = ((double) (conflicted[x] + almost_conflicted[x])) / ((double) age);
      double old_activity = Super::activity[x];
      Super::activity[x] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
      if (Super::order_heap.inHeap(x))
      {
         if (Super::activity[x] > old_activity)
            Super::order_heap.decrease(x);
         else
            Super::order_heap.increase(x);
      }
   }
#ifdef ANTI_EXPLORATION
   canceled[x] = Super::stat.conflicts;
#endif
   Super::notifyVarUnassigned(x);
}

template <typename Database>
inline void LrbBranch<Database>::notifyVarAssigned(const Var x)
{
   picked[x] = Super::stat.conflicts;
   conflicted[x] = 0;
   almost_conflicted[x] = 0;
#ifdef ANTI_EXPLORATION
   uint32_t age = Super::stat.conflicts - canceled[x];
   if (age > 0)
   {
      double decay = pow(0.95, age);
      Super::activity[x] *= decay;
      if (Super::order_heap.inHeap(x))
      Super::order_heap.increase(x);
   }
#endif

}

template <typename Database>
inline typename LrbBranch<Database>::Lit LrbBranch<Database>::pickBranchLit()
{
   Var next = var_Undef;

   // Activity based decision:
   while (next == var_Undef || Super::ig.value(next) != lbool::Undef() || !Super::decision[next])
      if (Super::order_heap.empty())
         return Lit::Undef();
      else
      {
#ifdef ANTI_EXPLORATION
         Var v = Super::order_heap[0];
         uint32_t age = Super::stat.conflicts - canceled[v];
         while (age > 0)
         {
            double decay = pow(0.95, age);
            Super::activity[v] *= decay;
            if (Super::order_heap.inHeap(v))
            Super::order_heap.increase(v);
            canceled[v] = Super::stat.conflicts;
            v = Super::order_heap[0];
            age = Super::stat.conflicts - canceled[v];
         }
#endif
         next = Super::order_heap.removeMin();
      }

   return Lit(next, Super::polarity[next]);
}
} /* namespace Minisat */

#endif /* SOURCES_BRANCH_LrbBranch_H_ */
