/*
 * Branch.h
 *
 *  Created on: 06.03.2020
 *      Author: hartung
 */

#ifndef SOURCES_BRANCH_BRANCH_H_
#define SOURCES_BRANCH_BRANCH_H_

#include "core/ImplicationGraph.h"
#include "initial/SolverConfig.h"
#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "utils/Random.h"
#include "core/Statistic.h"
#include "core/SolveMode.h"

namespace ctsat
{

template <typename Database>
class Branch
{
 public:

   typedef typename Database::Var Var;
   typedef typename Database::Lit Lit;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
   typedef typename Database::Clause Clause;

   struct BranchInputArgs
   {
      BranchInputArgs(
                      SolverConfig const & config,
                      SolveMode & smode,
                      Random & rand,
                      Statistic & stat,
                      Database & ca,
                      ImplicationGraph<Database> & ig);
      SolverConfig const & config;
      SolveMode & smode;
      Random & rand;
      Statistic & stat;
      Database & ca;
      ImplicationGraph<Database> & ig;
   };

   Branch(BranchInputArgs args);

   ~Branch();

   Var newVar(bool const dvar = true);

   void rebuildOrderHeap();

   void setDecisionVar(Var const v, bool const b);
   bool isDecisionVar(Var const v) const;
   uint64_t nDecVars() const;

   void notifyVarUnassigned(Var const v);

   void notifyConflictFound1(CRef const confl)
   {
   }
   void notifyConflictFound2(CRef const confl)
   {
   }
   void notifyVarSeenInConflict(Var const v)
   {
   }
   void notifyCreatedLearntClause(vec<Lit> const & out_learnt)
   {
   }
   void notifyConflictResolved()
   {
   }
   void notifyVarAssigned(Var const v)
   {
   }

   void notifyRestart()
   {
   }

   Lit pickBranchLit()
   {
      return Lit::Undef();
   }

 protected:
   bool const initVarPolZero;
   bool const initRndActivity;
   bool const initRndPolarity;
   uint64_t dec_vars;

   vec<char> polarity;         // The preferred polarity of each variable.
   vec<char> decision;  // Declares if a variable is eligible for selection in the decision heuristic
   vec<double> activity;     // A heuristic measurement of the activity of a variable.

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
   Heap<VarOrderLt> order_heap;  // A priority queue of variables ordered with respect to the variable activity.

   SolveMode & smode;
   Random & rand;
   Statistic & stat;
   Database & ca;
   ImplicationGraph<Database> & ig;

   void insertVarOrder(Var x);            // Insert a variable in the decision order priority queue.

};

template <typename Database>
inline Branch<Database>::Branch(BranchInputArgs args)
      : initVarPolZero(args.config.initVarPolZero),
        initRndActivity(args.config.rnd_active),
        initRndPolarity(args.config.rnd_polarity),
        dec_vars(0),
        polarity(),
        decision(),
        activity(),
        order_heap(VarOrderLt(activity)),
        smode(args.smode),
        rand(args.rand),
        stat(args.stat),
        ca(args.ca),
        ig(args.ig)
{
}

template <typename Database>
inline Branch<Database>::~Branch()
{
}
template <typename Database>
inline void Branch<Database>::rebuildOrderHeap()
{
   vec<Var> vs;
   for (Var v = 0; v < ig.nVars(); v++)
      if (decision[v] && ig.value(v) == lbool::Undef())
         vs.push(v);

   order_heap.build(vs);
}

template <typename Database>
inline bool Branch<Database>::isDecisionVar(Var const v) const
{
   return decision[v];
}

template <typename Database>
inline uint64_t Branch<Database>::nDecVars() const
{
   return dec_vars;
}

template <typename Database>
inline void Branch<Database>::setDecisionVar(Var const v, bool const b)
{
   if (b != decision[v])
      dec_vars += (b) ? 1 : -1;

   decision[v] = b;
   if (b && !order_heap.inHeap(v))
      order_heap.insert(v);
}

template <typename Database>
inline void Branch<Database>::insertVarOrder(Var const x)
{
   if (!order_heap.inHeap(x) && decision[x])
      order_heap.insert(x);
}

template <typename Database>
inline typename Branch<Database>::Var Branch<Database>::newVar(bool const dvar)
{
   Var const v = activity.size();
   double const activityVal = (initRndActivity) ? rand.drand() * 0.0001 : 0.0;
   bool const polarityVal = (initRndPolarity) ? rand.irand(2) : !initVarPolZero;
   activity.push(activityVal);
   polarity.push(polarityVal);
   decision.push(false);  // need to be false or dec_vars is wrong
   if (dvar)
      setDecisionVar(v, true);
   return v;
}

template <typename Database>
inline Branch<Database>::BranchInputArgs::BranchInputArgs(
                                                          SolverConfig const & config,
                                                          SolveMode & smode,
                                                          Random & rand,
                                                          Statistic & stat,
                                                          Database & ca,
                                                          ImplicationGraph<Database> & ig)
      : config(config),
        smode(smode),
        rand(rand),
        stat(stat),
        ca(ca),
        ig(ig)
{
}

template <typename Database>
inline void Branch<Database>::notifyVarUnassigned(const Var x)
{
   polarity[x] = ig[x].isFalse();
   insertVarOrder(x);
}

} /* namespace Minisat */

#endif /* SOURCES_BRANCH_BRANCH_H_ */
