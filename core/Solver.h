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

#ifndef Minisat_Solver_h
#define Minisat_Solver_h

#define BIN_DRUP

#define GLUCOSE23
//#define INT_QUEUE_AVG
//#define LOOSE_PROP_STAT

#ifdef GLUCOSE23
#define INT_QUEUE_AVG
#endif

#include "mtl/Vec.h"
#include "mtl/Alg.h"
#include "mtl/OccLists.h"
#include "mtl/AvgQueue.h"
#include "utils/Options.h"
#include "core/Statistic.h"
#include "core/SolveMode.h"
#include "initial/SolverConfig.h"
#include "core/ImplicationGraph.h"
#include "utils/DratPrint.h"
#include "utils/Random.h"
#include "minimize/Vivification.h"
#include "analyze/FirstUipAnalyze.h"

#include <vector>

namespace ctsat
{

template <typename TConnector, typename TDatabase, typename TBranch, typename TRestart,
      typename TReduce, typename TPropEngine, typename TAnalyze, typename TExchanger>
struct TemplateConfiguration
{
   typedef TConnector Connector;
   typedef TDatabase Database;
   typedef TBranch Branch;
   typedef TRestart Restart;
   typedef TReduce Reduce;
   typedef TPropEngine PropEngine;
   typedef TAnalyze Anaylze;
   typedef TExchanger Exchanger;
};

//=================================================================================================
// Solver -- the main class:

template <typename TemplateConfig>
class Solver
{

   Solver() = delete;
   Solver(Solver<TemplateConfig> const &) = delete;
   Solver<TemplateConfig> & operator=(Solver<TemplateConfig> const&) = delete;

 public:
   typedef typename TemplateConfig::Database Database;
   typedef typename Database::Lit Lit;
   typedef typename Database::Var Var;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;

   // Constructor/Destructor:
   //
   Solver(
          SolverConfig const & config,
          typename TemplateConfig::Connector & connector,
          SatInstance const & db);
   Solver(
          SolverConfig const & config,
          typename TemplateConfig::Connector & connector,
          SatInstance && db);
   virtual ~Solver();

   bool hasInterrupt() const;
   bool setOk(bool const b);
   bool isOk() const;

   Var newVar(bool const dvar);
   void removeVar(Var const v);

   void removeClause(CRef const cr);   // Detach and free a clause.

   void relocAll(typename TemplateConfig::Database& to);

   bool simplify();   // Removes already satisfied clauses.
   lbool solve();   // Search without assumptions.
   bool okay() const;   // FALSE means solver is in a conflicting state

   void setModel(vec<lbool> & model, vec<bool> const & isDecisionVar);
   int nAssigns() const;   // The current number of assigned literals.
   int nClauses() const;   // The current number of original clauses.
   int nLearnts() const;   // The current number of learnt clauses.
   int nVars() const;   // The current number of variables.
   int nFreeVars() const;
   typename TemplateConfig::Database & getClauseAllocator();

   void setRemoveSatisfied(bool const b);
   void allocExtra(bool const b);

   Statistic & getStatistic()
   {
      return stat;
   }

   // Resource contraints:
   //
   void setConfBudget(int64_t x);
   void setPropBudget(int64_t x);
   void interrupt();   // Trigger a (potentially asynchronous) interruption of the solver.
   void printFinalStats() const;

   // Memory managment:
   //
   virtual void garbageCollect();
   void checkGarbage(double gf);
   void checkGarbage();

   bool useVivification;
   int verbosity;
   double garbage_frac;  // The fraction of wasted memory allowed before a garbage collection is triggered.

   double learntsize_factor;  // The intitial limit for learnt clauses is a factor of the original clauses.                (default 1 / 3)
   double learntsize_inc;  // The limit for learnt clauses is multiplied with this factor each restart.                 (default 1.1)

   int learntsize_adjust_start_confl;
   double learntsize_adjust_inc;

   SolveMode smode;
   Random randEngine;
   Statistic stat;
   DratPrint<Lit> drat;
   typename TemplateConfig::Database ca;
   ImplicationGraph<decltype(ca)> ig;
   typename TemplateConfig::Branch branch;
   typename TemplateConfig::Restart restart;
   typename TemplateConfig::Reduce reduce;
   typename TemplateConfig::PropEngine propEngine;
   typename TemplateConfig::Anaylze analyze;
   Vivification<typename TemplateConfig::PropEngine> vivification;
   typename TemplateConfig::Exchanger exchange;

   bool ok;  // If FALSE, the constraints are already unsatisfiable. No part of the solver state may be used!
   bool asynch_interrupt;
   bool remove_satisfied;  // Indicates whether possibly inefficient linear scan for satisfied clauses should be performed in 'simplify'.
   bool steady_simplify;

   vec<CRef> clauses;   // List of problem clauses.

 protected:
   Solver(SolverConfig const & config, typename TemplateConfig::Connector & connector);

   // move from SatInstance:
   Solver(
          SolverConfig const & config,
          typename TemplateConfig::Connector & connector,
          vec<bool> const & decisionVars,
          Database && db,
          vec<CRef> && clauses,
          DratPrint<Lit> && drat);

   // copy from SatInstance
   Solver(
          SolverConfig const & config,
          typename TemplateConfig::Connector & connector,
          vec<bool> const & decisionVars,
          decltype(SatInstance::ca) const& db,
          vec<decltype(SatInstance::ca)::CRef> const & clauses);

   struct ConflictData
   {
      ConflictData()
            : nHighestLevel(-1),
              bOnlyOneLitFromHighest(false)
      {
      }

      int nHighestLevel;
      bool bOnlyOneLitFromHighest;
   };

   struct LitAssertion
   {
      int const level;
      Lit const l;
      CRef const cr;
      LitAssertion(int const level, Lit const l, CRef const cr)
            : level(level),
              l(l),
              cr(cr)
      {
      }
   };

   int confl_to_chrono;
   int chrono;

   std::vector<LitAssertion> analyzeAssertions;

   // Temporaries (to reduce allocation overhead). Each variable is prefixed by the method in which it is
   // used, exept 'seen' wich is used in several places.
   //
   vec<Lit> add_tmp;

   // Resource contraints:
   int64_t conflict_budget;   // -1 means no budget.
   int64_t propagation_budget;   // -1 means no budget.
   uint64_t nbSimplifyAll;

   uint64_t curSimplify;
   uint64_t nbconfbeforesimplify;
   int incSimplify;

   // returns minimum backtrack level
   int addLearntClauses();

   void uncheckedEnqueue(Lit const p, int const lvl = 0, CRef const from = Database::npos());
   bool enqueue(Lit const p, CRef const from = Database::npos());
   CRef propagate();
   void cancelUntil(int const lvl);
   bool resolveConflict(CRef const confl);

   // callbacks for analyze
   void varUsedInConflict(Var const v);
   void clauseUsedInConflict(CRef const ref);
   void clauseCreatedInConflict(vec<Lit> const & c);

   CRef addLearnt(LearntClause<Database> const & lc, bool const additional = false);

   lbool search();   // Search for a given number of conflicts.
   void removeSatisfied(vec<CRef>& cs);   // Shrink 'cs' to contain only non-satisfied clauses.

   ConflictData FindConflictLevel(CRef cind);
   bool importClauses();
   bool withinBudget() const;

   bool simplifyAll();
   bool simplifyClause(CRef const & ref);

   bool removed(CRef cr);

   void attachClauses();

   bool repropagateCurrentSolution();

   void newVars(vec<bool> const & decisionVars)
   {
      assert(ig.nVars() == 0);
      for (int i = 0; i < decisionVars.size(); ++i)
         newVar(decisionVars[i]);
   }

   int getBacktrackLevel(int const highestLevel, int const assertingLevel);

};

//=================================================================================================
// Implementation of inline methods:

template <typename TemplateConfig>
inline void Solver<TemplateConfig>::uncheckedEnqueue(Lit const p, int const lvl, CRef const from)
{
   propEngine.uncheckedEnqueue(branch, p, lvl, from);
}

template <typename TemplateConfig>
inline bool Solver<TemplateConfig>::enqueue(Lit const p, CRef const from)
{
   return propEngine.enqueue(branch, p, from);
}

template <typename TemplateConfig>
inline typename Solver<TemplateConfig>::CRef Solver<TemplateConfig>::propagate()
{
   return propEngine.propagate(branch);
}

template <typename TemplateConfig>
inline void Solver<TemplateConfig>::cancelUntil(int const lvl)
{
   propEngine.cancelUntil(branch, lvl);
}

template <typename TemplateConfig>
inline void Solver<TemplateConfig>::allocExtra(bool const b)
{
   ca.extra_clause_field = b;
}

template <typename TemplateConfig> inline typename TemplateConfig::Database & Solver<TemplateConfig>::getClauseAllocator()
{
   return ca;
}

template <typename TemplateConfig> inline void Solver<TemplateConfig>::setRemoveSatisfied(
                                                                                          bool const b)
{
   remove_satisfied = b;
}

template <typename TemplateConfig> inline bool Solver<TemplateConfig>::isOk() const
{
   return ok;
}

template <typename TemplateConfig> inline bool Solver<TemplateConfig>::setOk(bool const b)
{
   return ok = b;
}

template <typename TemplateConfig>
inline bool Solver<TemplateConfig>::hasInterrupt() const
{
   return asynch_interrupt;
}

template <typename TemplateConfig> inline void Solver<TemplateConfig>::checkGarbage(void)
{
   return checkGarbage(garbage_frac);
}

template <typename TemplateConfig> inline void Solver<TemplateConfig>::checkGarbage(double gf)
{
   if (ca.wasted() > ca.size() * gf)
      garbageCollect();
}

template <typename TemplateConfig> inline int Solver<TemplateConfig>::nAssigns() const
{
   return ig.nAssigns();
}

template <typename TemplateConfig> inline int Solver<TemplateConfig>::nClauses() const
{
   return clauses.size();
}

template <typename TemplateConfig> inline int Solver<TemplateConfig>::nLearnts() const
{
   return reduce.nClauses();
}

template <typename TemplateConfig> inline int Solver<TemplateConfig>::nVars() const
{
   return ig.nVars();
}

template <typename TemplateConfig> inline int Solver<TemplateConfig>::nFreeVars() const
{
   return (int) branch.nDecVars() - (ig.decisionLevel() == 0 ? ig.nAssigns() : ig.levelEnd(0));
}

template <typename TemplateConfig> inline void Solver<TemplateConfig>::setConfBudget(int64_t x)
{
   conflict_budget = stat.conflicts + x;
}

template <typename TemplateConfig> inline void Solver<TemplateConfig>::setPropBudget(int64_t x)
{
   propagation_budget = stat.propagations + x;
}

template <typename TemplateConfig> inline void Solver<TemplateConfig>::interrupt()
{
   asynch_interrupt = true;
}

template <typename TemplateConfig> inline bool Solver<TemplateConfig>::withinBudget() const
{
   return !asynch_interrupt
      && (conflict_budget < 0 || stat.conflicts < (uint64_t) conflict_budget)
      && (propagation_budget < 0 || stat.propagations < (uint64_t) propagation_budget);
}

}

#endif
