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

#include <math.h>
#include <algorithm>
#include <signal.h>
#include <unistd.h>

#include "mtl/Sort.h"
#include "core/Solver.h"
#include "utils/Logger.h"

//#define PRINT_OUT

namespace ctsat
{

//=================================================================================================
// Options:

//=================================================================================================
// Constructor/Destructor:

template <typename TemplateConfig>
Solver<TemplateConfig>::Solver(
                               SolverConfig const & config,
                               typename TemplateConfig::Connector & connector)
      : useVivification(config.useVivification),
        verbosity(config.verbosity),
        garbage_frac(config.garbage_frac),

        learntsize_factor(config.learntsize_factor),
        learntsize_inc(config.learntsize_inc),
        learntsize_adjust_start_confl(),
        learntsize_adjust_inc(),

        smode(),
        randEngine(config.rnd_seed),
        stat(),
        drat(config.drup, config.drup_file),
        ca(),
        ig(ca),
        branch(
              typename TemplateConfig::Branch::BranchInputArgs(config, smode, randEngine, stat, ca,
                                                               ig)),
        restart(config, smode, stat),
        reduce(config, stat, ca, ig),
        propEngine(stat, ca, ig),
        analyze(config, ca, ig, propEngine),
        vivification(config, ca, ig, propEngine),
        exchange(config, stat, ca, ig, connector, propEngine),
        ok(true),
        asynch_interrupt(false),
        remove_satisfied(config.remove_satisfied),
        steady_simplify(config.onlyExportWhenMin),
        clauses(),

        confl_to_chrono(config.confl_to_chrono),
        chrono(config.chrono),
        add_tmp(),

        conflict_budget(config.conflict_budget),
        propagation_budget(config.propagation_budget),

        nbSimplifyAll(0),

        curSimplify(1),
        nbconfbeforesimplify(config.nbconfbeforesimplify),
        incSimplify(config.incSimplify)
{
}

template <typename TemplateConfig>
Solver<TemplateConfig>::Solver(
                               SolverConfig const & config,
                               typename TemplateConfig::Connector & connector,
                               SatInstance const & inst)
      : Solver(config, connector, inst.model, inst.ca, inst.clauses)
{
}

template <typename TemplateConfig>
Solver<TemplateConfig>::Solver(
                               SolverConfig const & config,
                               typename TemplateConfig::Connector & connector,
                               SatInstance && inst)
      : Solver(config, connector, inst.model, std::move(inst.ca), std::move(inst.clauses),
               std::move(inst.drat))
{
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::attachClauses()
{
   for (int i = 0; i < clauses.size(); ++i)
      propEngine.attachClause(clauses[i]);
}

template <typename TemplateConfig>
Solver<TemplateConfig>::Solver(
                               SolverConfig const & config,
                               typename TemplateConfig::Connector & connector,
                               vec<lbool> const & model,
                               Database && db,
                               vec<CRef> && clauses,
                               DratPrint<Lit> && drat)
      : Solver(config, connector)
{
   this->drat = std::move(drat);
   newVars(model);
   this->clauses = std::move(clauses);
   ca = std::move(db);
   for (int i = 0; i < clauses.size(); ++i)
   {
      Clause const & c = ca[clauses[i]];
      for (int j = 0; j < c.size(); ++j)
         assert(ig.value(c[j]).isUndef());
   }
   attachClauses();
}
// copy from SatInstance
template <typename TemplateConfig>
Solver<TemplateConfig>::Solver(
                               SolverConfig const & config,
                               typename TemplateConfig::Connector & connector,
                               vec<decltype(SatInstance::ca)::lbool> const & model,
                               decltype(SatInstance::ca) const & db,
                               vec<decltype(SatInstance::ca)::CRef> const & clauses)
      : Solver(config, connector)
{
   assert(this->clauses.size() == 0);

   newVars(model);
   for (int i = 0; i < clauses.size(); ++i)
   {
      CRef const ref = ca.alloc(db[clauses[i]], false);
      this->clauses.push(ref);
      Clause const & c = ca[ref];
      for (int j = 0; j < c.size(); ++j)
         assert(ig.value(c[j]).isUndef());
      propEngine.attachClause(ref);
   }
}

template <typename TemplateConfig>
Solver<TemplateConfig>::~Solver()
{
}

template <typename TemplateConfig>
bool Solver<TemplateConfig>::removed(CRef cr)
{
   return ca[cr].mark() == 1;
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::setModel(vec<lbool> & model)
{
LOG("setting model")
                              assert(ig.nVars() <= model.size());
   for (Var i = 0; i < ig.nVars(); ++i)
   {
      lbool const val = ig.value(i);
      if (!val.isUndef())
         model[i] = val;
   }
}

// return true, if the clause cr should be removed
template <typename TemplateConfig>
bool Solver<TemplateConfig>::simplifyClause(CRef const & cr)
{
   Clause & c = ca[cr];
   if (removed(cr) || ig.satisfied(c))
   {
      ca.remove(cr);
      return true;
   } else if (c.simplified())
      return false;
   else
   {
      if(!propEngine.isAttached(cr))
         std::cout << "is bad attached: " << propEngine.isBadAttached(cr) << "\n";
      propEngine.detachClause(cr, true);
      if (vivification.run(c))
         drat.addClause(c);

      assert(c.size() > 0);
      if (c.size() == 1)
      {
         exchange.unitLearnt(c[0]);
         uncheckedEnqueue(c[0]);
         if (propagate() != CRef_Undef)
         {
            return setOk(false);
         }
         ca.remove(cr);
         return true;
      } else
      {
         propEngine.attachClause(cr);
         c.setSimplified(true);
         reduce.clauseImproved(cr);

      }
   }
   return false;
}

template <typename TemplateConfig>
inline bool Solver<TemplateConfig>::simplifyAll()
{
   assert(isOk());
   // simplify
   //
   if (useVivification && stat.conflicts >= curSimplify * nbconfbeforesimplify)
   {
      LOG("Simplifies clauses");
      if (propagate() != CRef_Undef)
         setOk(false);
      else
      {
         ++nbSimplifyAll;

         reduce.applyRemoveGoodClauses([&](CRef const & cr)
         {
            if (withinBudget() && isOk())
            {
               // ensure polling import clauses also when load on simplification is high
               if (exchange.shouldFetch())
               exchange.fetchClauses();
               return simplifyClause(cr);
            }
            else
            return false;
         });
         checkGarbage();
         if (steady_simplify)
         {
            nbconfbeforesimplify = stat.conflicts * 10ul * nbconfbeforesimplify;
            assert(curSimplify == 1);
         } else
         {
            curSimplify = (stat.conflicts / nbconfbeforesimplify) + 1;
            nbconfbeforesimplify += incSimplify;
         }
      }

      LOG("Simplify finished");
   }

   return isOk();
}

//=================================================================================================
// Minor methods:

// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//

template <typename TemplateConfig>
inline typename Solver<TemplateConfig>::Var Solver<TemplateConfig>::newVar(bool const dvar)
{
   Var const v = branch.newVar(dvar);
   ig.newVar();
   branch.setDecisionVar(v, dvar);
   propEngine.newVar();
   return v;
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::removeVar(Var const v)
{
   propEngine.removeVar(v);
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::removeClause(CRef const cr)
{
   Clause& c = ca[cr];
   if (c.mark() != 1)
      drat.removeClause(c);  // FIXME drat must be in database

   propEngine.detachClause(cr);
// Don't leave pointers to free'd memory!
   if (ig.locked(c))
   {
      Lit implied = c.size() != 2 ? c[0] : (ig.value(c[0]).isTrue() ? c[0] : c[1]);
      ig.reason(implied) = CRef_Undef;
   }
   ca.remove(cr);
}

template <typename TemplateConfig>
inline typename Solver<TemplateConfig>::ConflictData Solver<TemplateConfig>::FindConflictLevel(
                                                                                               CRef cind)
{
   ConflictData data;
   Clause& conflCls = ca[cind];
   data.nHighestLevel = ig.level(conflCls[0]);
   if (data.nHighestLevel == ig.decisionLevel() && ig.level(conflCls[1]) == ig.decisionLevel())
   {
      return data;
   }

   int highestId = 0;
   data.bOnlyOneLitFromHighest = true;
// find the largest decision level in the clause
   for (int nLitId = 1; nLitId < conflCls.size(); ++nLitId)
   {
      int const nLevel = ig.level(conflCls[nLitId]);
      if (nLevel > data.nHighestLevel)
      {
         highestId = nLitId;
         data.nHighestLevel = nLevel;
         data.bOnlyOneLitFromHighest = true;
      } else if (nLevel == data.nHighestLevel && data.bOnlyOneLitFromHighest == true)
         data.bOnlyOneLitFromHighest = false;

   }

   if (highestId != 0)
      propEngine.swapWatched(cind, 0, highestId);
   if (data.bOnlyOneLitFromHighest)
   {
      int secondHighest = ig.level(conflCls[1]);
      highestId = 1;
      for (int nLitId = 1; nLitId < conflCls.size(); ++nLitId)
      {
         int const nLevel = ig.level(conflCls[nLitId]);
         if (nLevel > secondHighest)
         {
            highestId = nLitId;
            secondHighest = nLevel;
         }
      }
      if (highestId != 1)
         propEngine.swapWatched(cind, 1, highestId);
      assert(ig.level(conflCls[1]) < ig.level(conflCls[0]));
      for (int i = 2; i < conflCls.size(); ++i)
         assert(ig.level(conflCls[i]) <= ig.level(conflCls[1]));
   }

   return data;
}

template <typename TemplateConfig>
inline void Solver<TemplateConfig>::clauseUsedInConflict(CRef const ref)
{
   reduce.clauseUsedInConflict(ref);
   exchange.clauseUsedInConflict(ref);
}

template <typename TemplateConfig>
inline void Solver<TemplateConfig>::varUsedInConflict(Var const v)
{
   branch.notifyVarSeenInConflict(v);
}
template <typename TemplateConfig>
inline void Solver<TemplateConfig>::clauseCreatedInConflict(vec<Lit> const & c)
{
   branch.notifyCreatedLearntClause(c);
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::removeSatisfied(vec<CRef>& cs)
{
   int i, j;
   for (i = j = 0; i < cs.size(); i++)
   {
      Clause& c = ca[cs[i]];
      if (ig.satisfied(c))
         removeClause(cs[i]);
      else
         cs[j++] = cs[i];
   }
   cs.shrink(i - j);
}

/*_________________________________________________________________________________________________
 |
 |  simplify : [void]  ->  [bool]
 |
 |  Description:
 |    Simplify the clause database according to the current top-level assigment. Currently, the only
 |    thing done here is the removal of satisfied clauses, but more things can be put here.
 |________________________________________________________________________________________________@*/
template <typename TemplateConfig>
bool Solver<TemplateConfig>::simplify()
{
   assert(ig.decisionLevel() == 0);

   if (!isOk() || propagate() != CRef_Undef)
   {
      return setOk(false);
   }

   if (nAssigns() == stat.simpDB_assigns || (stat.simpDB_props > 0))
      return true;

// Remove satisfied clauses:
   reduce.removeSatisfied([&](CRef const & ref)
   {  removeClause(ref);});
   if (remove_satisfied)        // Can be turned off.
      removeSatisfied(clauses);
   checkGarbage();
   branch.rebuildOrderHeap();

   stat.simpDB_assigns = nAssigns();
   stat.simpDB_props = stat.clauses_literals + stat.learnts_literals;  // (shouldn't depend on stats really, but it will do for now)

   return true;
}

template <typename TemplateConfig>
bool Solver<TemplateConfig>::importClauses()
{
LOG("Imports clauses")
            assert(ig.decisionLevel() == 0);
   Lit u;
   bool attached;
   CRef ref;
   exchange.fetchClauses();
   while ((u = exchange.getImportUnit()) != Lit::Undef())
      if (!enqueue(u))
         return setOk(false);
   while (std::get<1>((std::tie(attached, ref) = exchange.getImportClause())) != Database::npos())
   {
      reduce.addClause(ref);
      if (!attached
         && propEngine.safeAttachClause(ref) < ig.nVars()
         && ig.value(ca[ref][0]).isUndef())
      {
         uncheckedEnqueue(ca[ref][0], ref, 0);
         if (propagate() != Database::npos())
            return setOk(false);
      }
   }
   LOG("Imports clauses end")
   return exchange.isOk();
}
template <typename TemplateConfig>
inline int Solver<TemplateConfig>::getBacktrackLevel(
                                                     int const highestLevel,
                                                     int const assertingLevel)
{
   bool const shouldChrono = (confl_to_chrono < 0
      || static_cast<unsigned>(confl_to_chrono) <= stat.conflicts)
      && chrono > -1
      && (ig.decisionLevel() - assertingLevel) >= chrono;
   int backtrackLevel = assertingLevel;

   if (shouldChrono)
   {
      backtrackLevel = highestLevel - 1;
      ++stat.chrono_backtrack;
   }
   stat.non_chrono_backtrack += !shouldChrono;
//   LOG("Backtrack to " + std::to_string(backtrackLevel))
   return backtrackLevel;
}

template <typename VarVisit, typename ClauseVisit, typename ClauseCreated>
struct AnalyseCallbacks
{
   VarVisit const onVarVisit;
   ClauseVisit const onClauseVisit;
   ClauseCreated const onClauseCreated;

   AnalyseCallbacks(
                    VarVisit const & onVarVisit,
                    ClauseVisit const & onClauseVisit,
                    ClauseCreated const & onClauseCreated);
};
template <typename VarVisit, typename ClauseVisit, typename ClauseCreated>
inline AnalyseCallbacks<VarVisit, ClauseVisit, ClauseCreated>::AnalyseCallbacks(
                                                                                VarVisit const & onVarVisit,
                                                                                ClauseVisit const & onClauseVisit,
                                                                                ClauseCreated const & onClauseCreated)
      : onVarVisit(onVarVisit),
        onClauseVisit(onClauseVisit),
        onClauseCreated(onClauseCreated)
{
}

// C++11 workarround:
template <typename VarVisit, typename ClauseVisit, typename ClauseCreated>
inline AnalyseCallbacks<VarVisit, ClauseVisit, ClauseCreated> getCallbacks(
                                                                           VarVisit const & onVarVisit,
                                                                           ClauseVisit const & onClauseVisit,
                                                                           ClauseCreated const & onClauseCreated)
{
   return AnalyseCallbacks<VarVisit, ClauseVisit, ClauseCreated>(onVarVisit, onClauseVisit,
                                                                 onClauseCreated);
}

template <typename TemplateConfig>
inline bool Solver<TemplateConfig>::resolveConflict(CRef const confl)
{
   auto const callbacks = getCallbacks([&](Var const v)
   {  varUsedInConflict(v);},
                                       [&](CRef const cref)
                                       {  clauseUsedInConflict(cref);},
                                       [&](vec<Lit> const & c)
                                       {  clauseCreatedInConflict(c);});
   LOG("Conflict found")
   int assertingLevel;
   analyzeAssertions.clear();
   reduce.adjustOnConflict();
   branch.notifyConflictFound1(confl);
   restart.notifyConflictFound();
   exchange.conflictFound();

   ConflictData const conflInfo = FindConflictLevel(confl);
//   LOG("Conflict level " + std::to_string(conflInfo.nHighestLevel))
   if (conflInfo.nHighestLevel == 0)
   {
      LOG("Level 0 conflict")
      return setOk(false);
   }

   if (conflInfo.bOnlyOneLitFromHighest)
   {
      Clause const & c = ca[confl];
      analyzeAssertions.emplace_back(ig.level(c[1]), c[0], confl);
      assertingLevel = analyzeAssertions.back().level;
   } else
   {

      branch.notifyConflictFound2(confl);

      analyze.run(confl, callbacks);
      assertingLevel = addLearntClauses();

      branch.notifyConflictResolved();
      reduce.notifyConflictResolved();
   }
   cancelUntil(getBacktrackLevel(conflInfo.nHighestLevel, assertingLevel));
   stat.level_backtracked += conflInfo.nHighestLevel - ig.decisionLevel();

   assert(analyzeAssertions.size() > 0);
   for (auto const & prop : analyzeAssertions)
      if (prop.level <= ig.decisionLevel() && ig.value(prop.l).isUndef())
      {
         if (prop.cr != Database::npos())
         {
            Clause const & c = ca[prop.cr];
            assert(prop.l == c[0]);
            for (int i = 1; i < c.size(); ++i)
               assert(ig.value(c[i]).isFalse());
         }
         uncheckedEnqueue(prop.l, prop.level, prop.cr);
      }
   LOG("Conflict resolved")
   return true;
}

template <typename TemplateConfig>
inline int Solver<TemplateConfig>::addLearntClauses()
{
   int assertingLevel = ig.decisionLevel(), lbd = std::numeric_limits<int>::max(), nAdded = 0;
   assert(analyze.hasLearntClause());

   if (analyze.learntsMultipleClauses())
      while (analyze.hasLearntClause())
      {
         LearntClause<Database> const & lc = analyze.getLearntClause();
         CRef const cr = addLearnt(lc);
         ++nAdded;
         if (lc.isAsserting)
         {
            analyzeAssertions.emplace_back((lc.c.size() > 1) ? ig.level(lc.c[1]) : 0, lc.c[0], cr);
            assertingLevel = std::min(assertingLevel, analyzeAssertions.back().level);
            assert(
                  analyzeAssertions.back().level > 0
                     || analyzeAssertions.back().cr == Database::npos());
            lbd = std::min(lbd, lc.lbd);
         }
      }
   else
   {
      nAdded = 1;
      LearntClause<Database> const & lc = analyze.getLearntClause();
      analyzeAssertions.emplace_back((lc.c.size() > 1) ? ig.level(lc.c[1]) : 0, lc.c[0],
                                     addLearnt(lc));
      assertingLevel = analyzeAssertions.back().level;
      lbd = lc.lbd;
   }
//   LOG("Added " + std::to_string(nAdded) + " clauses")
   stat.nAdditionalLearnt += nAdded - 1;
   restart.clauseLearnt(lbd);  // only use heuristic on one clause
   return assertingLevel;
}

template <typename TemplateConfig>
typename Solver<TemplateConfig>::CRef Solver<TemplateConfig>::addLearnt(
                                                                        LearntClause<Database> const & lc)
{
   CRef cr = Database::npos();
   drat.addClause(lc.c);
   if (lc.c.size() == 1)
      exchange.unitLearnt(lc.c[0]);
   else
   {
      cr = ca.alloc(lc.c, true);
      Clause & c = ca[cr];
      c.set_lbd(lc.lbd);
      propEngine.attachClause(cr);
      reduce.addClause(cr);
      exchange.clauseLearnt(cr);
   }
   return cr;
}

/*_________________________________________________________________________________________________
 |
 |  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
 |
 |  Description:
 |    Search for a model the specified number of conflicts.
 |
 |  Output:
 |    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
 |    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
 |    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
 |________________________________________________________________________________________________@*/
template <typename TemplateConfig>
typename Solver<TemplateConfig>::lbool Solver<TemplateConfig>::search()
{
   assert(isOk());
   bool firstAfterConflict = false;
   while (true)
   {
      CRef const confl = propagate();

      if (confl != CRef_Undef)
      {
         ++stat.conflicts;
         if (verbosity > 0 && stat.conflicts % 20000 == 0)
            stat.print();
         if (!resolveConflict(confl))
            return lbool::False();
         firstAfterConflict = true;
         continue;
      }

      if (firstAfterConflict)
      {
         firstAfterConflict = false;
         if (restart.shouldRestart() || exchange.isFinished() || !withinBudget())
         {
            cancelUntil(0);
            return lbool::Undef();
         }

         if (ig.decisionLevel() == 0)
         {
            if (!importClauses() || !simplify())
               return lbool::False();
         }

         if (reduce.run([&](CRef const ref)
         {  removeClause(ref);}))
            checkGarbage();
      }

      Lit const next = branch.pickBranchLit();
      if (next == Lit::Undef())
      {
         // Model found:
         return lbool::True();
      }

// Increase decision level and enqueue 'next'
      ig.newDecisionLevel();
      uncheckedEnqueue(next, ig.decisionLevel());
      ++stat.decisions;
   }
}

// NOTE: assumptions passed in member-variable 'assumptions'.
template <typename TemplateConfig>
typename Solver<TemplateConfig>::lbool Solver<TemplateConfig>::solve()
{
   LOG("Starts solving")
   if (!isOk())
      return lbool::False();
   lbool status = lbool::Undef();

   add_tmp.clear();

// Search:
   while (status == lbool::Undef() && withinBudget() && !exchange.isFinished())
   {
LOG("Start")
                                       assert(ig.decisionLevel() == 0);
      branch.notifyRestart();
      restart.notifyRestart();
      analyze.notifyRestart();

      importClauses();
      if (!exchange.isOk() || !isOk() || !simplifyAll())
         status = lbool::False();
      else
      {
         ++stat.restarts;
         status = search();
      }
   }

   if (verbosity >= 1)
      printf("c ===============================================================================\n");

   if (status != lbool::Undef() && exchange.setFinished(status))
   {
      LOG("Solver won")
      if (status == lbool::False())
      {
         if (!exchange.isOk() && verbosity > -1)
            std::cout << "c solved through import\n";

         drat.addEmptyClause();
         drat.flush();
      }
   } else
      status = lbool::Undef();
   LOG("Ends solving")
   return status;
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::printFinalStats() const
{
   stat.printFinal();
}

//=================================================================================================
// Garbage Collection methods:

template <typename TemplateConfig>
void Solver<TemplateConfig>::relocAll(typename TemplateConfig::Database& to)
{
// All watchers:
   propEngine.relocAll(to);

// All reasons:
//
   for (int i = 0; i < ig.nAssigns(); i++)
   {
      Var const v = ig.getTrailLit(i).var();
      CRef const cr = ig.reason(v);
      if (cr != CRef_Undef && (ca[cr].reloced() || ig.locked(cr)))
         ca.reloc(ig.reason(v), to);
   }

// All learnt:
//
   reduce.relocAll(to);

// All original:
//
   int i, j;
   for (i = j = 0; i < clauses.size(); i++)
      if (ca[clauses[i]].mark() != 1)
      {
         ca.reloc(clauses[i], to);
         clauses[j++] = clauses[i];
      }

   exchange.relocAll(to);
   clauses.shrink(i - j);
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::garbageCollect()
{
// Initialize the next region to a size corresponding to the estimated utilization degree. This
// is not precise but should avoid some unnecessary reallocations for the new region:
   ClauseAllocator to(ca.size() - ca.wasted());

   relocAll(to);
   if (verbosity >= 2)
      printf("c |  Garbage collection:   %12d bytes => %12d bytes             |\n",
             ca.size() * ClauseAllocator::Unit_Size, to.size() * ClauseAllocator::Unit_Size);
   to.moveTo(ca);
}
}
