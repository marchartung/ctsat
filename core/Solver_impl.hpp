/*****************************************************************************************[Main.cc]
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

//#define PRINT_OUT

namespace CTSat
{

//=================================================================================================
// Options:

//=================================================================================================
// Constructor/Destructor:

template <typename TemplateConfig>
Solver<TemplateConfig>::Solver(
                               SolverConfig const & config,
                               typename TemplateConfig::Connector & connector)
      : verbosity(config.verbosity),
        ccmin_mode(config.ccmin_mode),
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
        exchange(config, stat, ca, ig, connector, propEngine),
        ok(true),
        asynch_interrupt(false),
        remove_satisfied(config.remove_satisfied),
        clauses(),

        confl_to_chrono(config.confl_to_chrono),
        chrono(config.chrono),
        analyze_stack(),
        add_tmp(),

        conflict_budget(config.conflict_budget),
        propagation_budget(config.propagation_budget),

        nbSimplifyAll(0),

        simp_learnt_clause(),

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
   assert(ig.nVars() <= model.size());
   for (Var i = 0; i < ig.nVars(); ++i)
   {
      lbool const val = ig.value(i);
      if (!val.isUndef())
         model[i] = val;
   }
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::simpleAnalyze(CRef confl, vec<Lit>& out_learnt, bool True_confl)
{
   int pathC = 0;
   Lit p = Lit::Undef();
   int index = ig.nAssigns() - 1;

   do
   {
      if (confl != CRef_Undef)
      {
         Clause& c = ca[confl];
         // Special case for binary clauses
         // The first one has to be SAT
         if (p != Lit::Undef() && c.size() == 2 && ig.value(c[0]).isFalse())
         {

            assert(ig.value(c[1]).isTrue());
            Lit tmp = c[0];
            c[0] = c[1], c[1] = tmp;
         }
         // if True_confl==true, then choose p begin with the 1th index of c;
         for (int j = (p == Lit::Undef() && True_confl == false) ? 0 : 1; j < c.size(); j++)
         {
            Var const v = c[j].var();
            if (!ig.isSeen(v))
            {
               ig.setSeen(v);
               pathC++;
            }
         }
      } else if (confl == CRef_Undef)
      {
         out_learnt.push(~p);
      }
      // if not break, while() will come to the index of trail blow 0, and fatal error occur;
      if (pathC == 0)
         break;
      // Select next clause to look at:
      while (!ig.isSeen(ig.getTrailLit(index--).var()))
         ;
      // if the reason cr from the 0-level assigned var, we must break avoid move forth further;
      // but attention that maybe seen[x]=1 and never be clear. However makes no matter;
      if (propEngine.trailRecord > index + 1)
         break;
      p = ig.getTrailLit(index + 1);
      Var const v = p.var();
      confl = ig.reason(v);
      ig.unsetSeen(v);
      pathC--;

   } while (pathC >= 0);
}

// return true, if the clause cr should be removed
template <typename TemplateConfig>
bool Solver<TemplateConfig>::simplifyClause(CRef const cr)
{
   if(exchange.shouldFetch())
      exchange.fetchClauses();
   bool sat, false_lit;
   Clause & c = ca[cr];

   if (removed(cr))
      return true;
   else if (c.simplified() || !isOk())
      return false;
   else
   {
      int saved_size = c.size();
      sat = false_lit = false;
      for (int i = 0; i < c.size(); i++)
      {
         lbool const val = ig.value(c[i]);
         if (val.isTrue())
         {
            sat = true;
            break;
         } else if (val.isFalse())
            false_lit = true;
      }
      if (sat)
      {
         removeClause(cr);
         return true;
      } else
      {
         propEngine.detachClause(cr, true);

         if (false_lit)
         {
            int li = 0, lj = 0;
            for (; li < c.size(); li++)
            {
               if (!ig.value(c[li]).isFalse())
               {
                  c[lj++] = c[li];
               }
            }
            c.shrink(li - lj);
         }

         assert(c.size() > 1);
         // simplify a learnt clause c
         simplifyLearnt(c);
         assert(c.size() > 0);

         if (c.size() == 1)
         {
            // when unit clause occur, enqueue and propagate
            exchange.unitLearnt(c[0]);
            drat.addClause(c);
            uncheckedEnqueue(c[0]);
            if (propagate() != CRef_Undef)
            {
               return setOk(false);
            }
            // delete the clause memory in logic
            c.mark(1);
            ca.free(cr);
            return true;
         } else
         {
            propEngine.attachClause(cr);

            int nblevels = ig.computeLBD(c);
            if (nblevels < c.lbd())
               c.set_lbd(nblevels);

            c.setSimplified(true);
            reduce.clauseImproved(cr);
            if (saved_size != c.size())
            {
               exchange.clauseImproved(cr);
               drat.addClause(c);
            }

         }
      }
   }
   return false;
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::simplifyLearnt(Clause& c)
{

   propEngine.trailRecord = ig.nAssigns();      // record the start pointer

   vec<Lit> falseLit;
   falseLit.clear();

   //sort(&c[0], c.size(), VarOrderLevelLt(vardata));

   bool True_confl = false;
   int i, j;
   CRef confl;

   for (i = 0, j = 0; i < c.size(); i++)
   {
      if (ig.value(c[i]).isUndef())
      {
         propEngine.simpleUncheckEnqueue(~c[i]);
         c[j++] = c[i];
         confl = propEngine.simplePropagate();
         if (confl != CRef_Undef)
         {
            break;
         }
      } else
      {
         if (ig.value(c[i]).isTrue())
         {
            c[j++] = c[i];
            True_confl = true;
            confl = ig.reason((c[i].var()));
            break;
         } else
         {
            falseLit.push(c[i]);
         }
      }
   }
   c.shrink(c.size() - j);

   if (confl != CRef_Undef || True_confl == true)
   {
      simp_learnt_clause.clear();
      if (True_confl == true)
      {
         simp_learnt_clause.push(c.last());
      }
      simpleAnalyze(confl, simp_learnt_clause, True_confl);

      if (simp_learnt_clause.size() < c.size())
      {
         for (i = 0; i < simp_learnt_clause.size(); i++)
         {
            c[i] = simp_learnt_clause[i];
         }
         c.shrink(c.size() - i);
      }
   }
   propEngine.cancelUntilTrailRecord();
}

template <typename TemplateConfig>
bool Solver<TemplateConfig>::simplifyAll()
{

   if (!isOk() || propagate() != CRef_Undef)
      return setOk(false);

   reduce.applyRemoveGoodClauses([&](CRef const cr)
   {  return simplifyClause(cr);});
   checkGarbage();

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
void Solver<TemplateConfig>::removeClause(CRef cr)
{
   Clause& c = ca[cr];

   if (c.mark() != 1)
      drat.removeClause(c);
   else
      printf("c Bug. I don't expect this to happen.\n");

   propEngine.detachClause(cr);
   // Don't leave pointers to free'd memory!
   if (ig.locked(c))
   {
      Lit implied = c.size() != 2 ? c[0] : (ig.value(c[0]).isTrue() ? c[0] : c[1]);
      ig.reason((implied.var())) = CRef_Undef;
   }
   c.mark(1);
   ca.free(cr);
}

template <typename TemplateConfig>
inline typename Solver<TemplateConfig>::ConflictData Solver<TemplateConfig>::FindConflictLevel(
                                                                                               CRef cind)
{
   ConflictData data;
   Clause& conflCls = ca[cind];
   data.nHighestLevel = ig.level((conflCls[0].var()));
   if (data.nHighestLevel == ig.decisionLevel()
      && ig.level((conflCls[1].var())) == ig.decisionLevel())
   {
      return data;
   }

   int highestId = 0;
   data.bOnlyOneLitFromHighest = true;
   // find the largest decision level in the clause
   for (int nLitId = 1; nLitId < conflCls.size(); ++nLitId)
   {
      int const nLevel = ig.level((conflCls[nLitId].var()));
      if (nLevel > data.nHighestLevel)
      {
         highestId = nLitId;
         data.nHighestLevel = nLevel;
         data.bOnlyOneLitFromHighest = true;
      } else if (nLevel == data.nHighestLevel && data.bOnlyOneLitFromHighest == true)
      {
         data.bOnlyOneLitFromHighest = false;
      }
   }

   if (highestId != 0)
      propEngine.swapWatched(cind, 0, highestId);
   if (data.bOnlyOneLitFromHighest)
   {
      int secondHighest = ig.level((conflCls[1].var()));
      highestId = 1;
      for (int nLitId = 1; nLitId < conflCls.size(); ++nLitId)
      {
         int const nLevel = ig.level((conflCls[nLitId].var()));
         if (nLevel > secondHighest)
         {
            highestId = nLitId;
            secondHighest = nLevel;
         }
      }
      if (highestId != 1)
         propEngine.swapWatched(cind, 1, highestId);
      assert(ig.level(conflCls[1].var()) < ig.level(conflCls[0].var()));
      for (int i = 2; i < conflCls.size(); ++i)
         assert(ig.level(conflCls[i].var()) <= ig.level(conflCls[1].var()));
   }

   return data;
}
template <typename TemplateConfig>
inline void Solver<TemplateConfig>::clauseUsedInConflict(CRef const ref)
{
   reduce.clauseUsedInConflict(ref);
   exchange.clauseUsedInConflict(ref);
}

/*_________________________________________________________________________________________________
 |
 |  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
 |
 |  Description:
 |    Analyze conflict and produce a reason clause.
 |
 |    Pre-conditions:
 |      * 'out_learnt' is assumed to be cleared.
 |      * Current decision level must be greater than root level.
 |
 |    Post-conditions:
 |      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
 |      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the
 |        rest of literals. There may be others from the same level though.
 |
 |________________________________________________________________________________________________@*/
template <typename TemplateConfig>
void Solver<TemplateConfig>::analyze(
                                     CRef confl,
                                     vec<Lit>& out_learnt,
                                     int& out_btlevel,
                                     int& out_lbd)
{
   int pathC = 0;
   Lit p = Lit::Undef();

   // Generate conflict clause:
   //
   out_learnt.push();      // (leave room for the asserting literal)
   int index = ig.nAssigns() - 1;
   int nDecisionLevel = ig.level((ca[confl][0].var()));
   assert(nDecisionLevel == ig.level((ca[confl][0].var())));

   do
   {
      assert(confl != CRef_Undef);  // (otherwise should be UIP)
      Clause & c = ca[confl];

      // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
      if (p != Lit::Undef() && c.size() == 2 && ig.value(c[0]).isFalse())
      {
         assert(ig.value(c[1]).isTrue());
         Lit const tmp = c[0];
         c[0] = c[1], c[1] = tmp;
      }

      clauseUsedInConflict(confl);

      for (int j = (p == Lit::Undef()) ? 0 : 1; j < c.size(); j++)
      {
         Lit const q = c[j];
         Var const v = c[j].var();

         if (!ig.isSeen(v) && ig.level(v) > 0)
         {
            branch.notifyVarSeenInConflict(v);
            ig.setSeen(v);
            if (ig.level(v) >= nDecisionLevel)
            {
               pathC++;
            } else
            {
               out_learnt.push(q);
               ig.markSeenToClear(q);
            }
         }
      }

      // Select next clause to look at:
      do
      {
         while (!ig.isSeen(ig.getTrailLit(index--).var()))
            ;
         p = ig.getTrailLit(index + 1);
      } while (ig.level((p.var())) < nDecisionLevel);

      Var const pv = p.var();
      confl = ig.reason(pv);
      ig.unsetSeen(pv);
      pathC--;

   } while (pathC > 0);
   out_learnt[0] = ~p;

   // Simplify conflict clause:
   //
   int i, j;

   if (ccmin_mode == 2)
   {
      uint32_t abstract_level = 0;
      for (i = 1; i < out_learnt.size(); i++)
         abstract_level |= ig.abstractLevel((out_learnt[i].var()));  // (maintain an abstraction of levels involved in conflict)

      for (i = j = 1; i < out_learnt.size(); i++)
         if (ig.reason((out_learnt[i].var())) == CRef_Undef
            || !litRedundant(out_learnt[i], abstract_level))
            out_learnt[j++] = out_learnt[i];

   } else if (ccmin_mode == 1)
   {
      for (i = j = 1; i < out_learnt.size(); i++)
      {
         Var const x = (out_learnt[i].var());

         if (ig.reason(x) == CRef_Undef)
            out_learnt[j++] = out_learnt[i];
         else
         {
            Clause& c = ca[ig.reason((out_learnt[i].var()))];
            for (int k = c.size() == 2 ? 0 : 1; k < c.size(); k++)
               if (!ig.isSeen(c[k].var()) && ig.level((c[k].var())) > 0)
               {
                  out_learnt[j++] = out_learnt[i];
                  break;
               }
         }
      }
   } else
      i = j = out_learnt.size();

   stat.max_literals += out_learnt.size();
   out_learnt.shrink(i - j);
   stat.tot_literals += out_learnt.size();

   out_lbd = ig.computeLBD(out_learnt);
   if (out_lbd <= 6 && out_learnt.size() <= 30)  // Try further minimization?
      if (propEngine.binResMinimize(out_learnt))
         out_lbd = ig.computeLBD(out_learnt);  // Recompute LBD if minimized.

   // Find correct backtrack level:
   //
   if (out_learnt.size() == 1)
      out_btlevel = 0;
   else
   {
      int max_i = 1;
      // Find the first literal assigned at the next-highest level:
      for (int i = 2; i < out_learnt.size(); i++)
         if (ig.level((out_learnt[i].var())) > ig.level((out_learnt[max_i].var())))
            max_i = i;
      // Swap-in this literal at index 1:
      Lit const p = out_learnt[max_i];
      out_learnt[max_i] = out_learnt[1];
      out_learnt[1] = p;
      out_btlevel = ig.level((p.var()));
   }

   branch.notifyCreatedLearntClause(out_learnt);

   ig.clearSeen();
}

// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
template <typename TemplateConfig>
bool Solver<TemplateConfig>::litRedundant(Lit p, uint32_t abstract_levels)
{
   analyze_stack.clear();
   analyze_stack.push(p);
   int const top = ig.nToClear();
   while (analyze_stack.size() > 0)
   {
      assert(ig.reason((analyze_stack.last().var())) != CRef_Undef);
      Clause& c = ca[ig.reason((analyze_stack.last().var()))];
      analyze_stack.pop();

      // Special handling for binary clauses like in 'analyze()'.
      if (c.size() == 2 && ig.value(c[0]).isFalse())
      {
         assert(ig.value(c[1]).isTrue());
         Lit const tmp = c[0];
         c[0] = c[1], c[1] = tmp;
      }

      for (int i = 1; i < c.size(); i++)
      {
         Lit const p = c[i];
         Var const v = c[i].var();
         if (!ig.isSeen(v) && ig.level(v) > 0)
         {
            if (ig.reason(v) != CRef_Undef && (ig.abstractLevel(v) & abstract_levels) != 0)
            {
               ig.setSeen(v);
               ig.markSeenToClear(p);
               analyze_stack.push(p);

            } else
            {
               ig.clearSeen(top);
               return false;
            }
         }
      }
   }

   return true;
}

/*_________________________________________________________________________________________________
 |
 |  analyzeFinal : (p : Lit)  ->  [void]
 |
 |  Description:
 |    Specialized analysis procedure to express the final conflict in terms of assumptions.
 |    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
 |    stores the result in 'out_conflict'.
 |________________________________________________________________________________________________@*/
template <typename TemplateConfig>
void Solver<TemplateConfig>::analyzeFinal(Lit p, vec<Lit>& out_conflict)
{
   out_conflict.clear();
   out_conflict.push(p);

   if (ig.decisionLevel() == 0)
      return;

   ig.setSeen(p.var());

   for (int i = ig.nAssigns() - 1; i >= ig.levelEnd(0); i--)
   {
      Lit const l = ig.getTrailLit(i);
      Var x = l.var();
      if (ig.isSeen(x))
      {
         if (ig.reason(x) == CRef_Undef)
         {
            assert(ig.level(x) > 0);
            out_conflict.push(~l);
         } else
         {
            Clause& c = ca[ig.reason(x)];
            for (int j = c.size() == 2 ? 0 : 1; j < c.size(); j++)
               if (ig.level(c[j].var()) > 0)
                  ig.setSeen(c[j].var());
         }
         ig.unsetSeen(x);
      }
   }
   ig.unsetSeen(p.var());
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

template <typename SolverType>
struct UIPOrderByILevel_Lt
{
   SolverType& solver;
   const vec<double>& var_iLevel;
   bool operator ()(Lit x, Lit y) const
   {
      return var_iLevel[x.var()] < var_iLevel[y.var()]
         || (var_iLevel[x.var()] == var_iLevel[y.var()]
            && solver.level(x.var()) > solver.level(y.var()));
   }
   UIPOrderByILevel_Lt(const vec<double>& iLevel, SolverType& para_solver)
         : solver(para_solver),
           var_iLevel(iLevel)
   {
   }
};

template <typename TemplateConfig>
bool Solver<TemplateConfig>::importClauses()
{
   assert(ig.decisionLevel() == 0);
   Lit u;
   bool attached;
   CRef ref;
   exchange.fetchClauses();
   while ((u = exchange.getImportUnit()) != Lit::Undef())
      if (!enqueue(u))
      {
         std::cout << "c unit import conflict not caught\n";
      }
   while (std::get<1>((std::tie(attached, ref) = exchange.getImportClause())) != Database::npos())
   {
      reduce.addClause(ref);
      if (!attached
         && propEngine.safeAttachClause(ref) < ig.nVars()
         && ig.value(ca[ref][0]).isUndef())
      {
         uncheckedEnqueue(ca[ref][0], ref, 0);
      }
   }
   return exchange.isOk();
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
   int backtrack_level;
   int lbd;
   vec<Lit> learnt_clause;

   stat.starts++;
   importClauses();
   if (!exchange.isOk())
   {
      return lbool::False();
   }
   // simplify
   //
   if (stat.conflicts >= curSimplify * nbconfbeforesimplify)
   {
      nbSimplifyAll++;
      if (!simplifyAll())
      {
         return lbool::False();
      }
      curSimplify = (stat.conflicts / nbconfbeforesimplify) + 1;
      nbconfbeforesimplify += incSimplify;
   }

   for (;;)
   {
      CRef confl = propagate();

      if (confl != CRef_Undef)
      {
         // CONFLICT
         stat.conflicts++;
         reduce.adjustOnConflict();
         branch.notifyConflictFound1(confl);
         restart.notifyConflictFound();
         ConflictData data = FindConflictLevel(confl);
         if (data.nHighestLevel == 0)
         {
            setOk(false);
            return lbool::False();
         }
         if (exchange.isFinished() || !withinBudget())
            return lbool::Undef();

         if (data.bOnlyOneLitFromHighest)
         {
            Clause const & c = ca[confl];
            int const blevel = (chrono > -1) ? data.nHighestLevel - 1 : ig.level(c[1].var());
            cancelUntil(blevel);
            uncheckedEnqueue(c[0], ig.level(c[1].var()), confl);
            continue;
         }
         branch.notifyConflictFound2(confl);

         learnt_clause.clear();
         analyze(confl, learnt_clause, backtrack_level, lbd);

         // check chrono backtrack condition
         if ((confl_to_chrono < 0 || static_cast<unsigned>(confl_to_chrono) <= stat.conflicts)
            && chrono > -1
            && (ig.decisionLevel() - backtrack_level) >= chrono)
         {
            ++stat.chrono_backtrack;
            cancelUntil(data.nHighestLevel - 1);
         } else  // default behavior
         {
            ++stat.non_chrono_backtrack;
            cancelUntil(backtrack_level);
         }

         lbd--;

         if (learnt_clause.size() == 1)
         {
            exchange.unitLearnt(learnt_clause[0]);
            uncheckedEnqueue(learnt_clause[0]);
         } else
         {
            CRef cr = ca.alloc(learnt_clause, true);
            Clause & c = ca[cr];
            c.set_lbd(lbd);
            reduce.addClause(cr);
            exchange.clauseLearnt(cr);
            propEngine.attachClause(cr);

            uncheckedEnqueue(learnt_clause[0], backtrack_level, cr);

         }
         restart.notifyConflictResolved(lbd);

         drat.addClause(learnt_clause);
         branch.notifyConflictResolved();
         reduce.notifyConflictResolved();

      } else
      {
         // NO CONFLICT
         if (restart.shouldRestart())
         {
            cancelUntil(0);
            return lbool::Undef();
         }

         // Simplify the set of problem clauses:
         if (ig.decisionLevel() == 0 && !simplify())
            return lbool::False();

         if (reduce.run([&](CRef const ref)
         {  removeClause(ref);}))
            checkGarbage();

         Lit next = Lit::Undef();
         // New variable decision:
         stat.decisions++;
         next = branch.pickBranchLit();

         if (next == Lit::Undef())
            // Model found:
            return lbool::True();

         // Increase decision level and enqueue 'next'
         ig.newDecisionLevel();
         uncheckedEnqueue(next, ig.decisionLevel());
      }
   }
}

// NOTE: assumptions passed in member-variable 'assumptions'.
template <typename TemplateConfig>
typename Solver<TemplateConfig>::lbool Solver<TemplateConfig>::solve()
{
   if (!isOk())
      return lbool::False();
   lbool status = lbool::Undef();

   if (verbosity >= 1)
   {
      printf("c ============================[ Search Statistics ]==============================\n");
      printf("c | Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
      printf("c |           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
      printf("c ===============================================================================\n");
   }

   add_tmp.clear();

// Search:
   while (status == lbool::Undef() && withinBudget() && !exchange.isFinished())
   {

      assert(ig.decisionLevel() == 0);
      branch.notifyRestart();
      restart.notifyRestart();
      status = search();
   }

   if (verbosity >= 1)
      printf("c ===============================================================================\n");

   if (status != lbool::Undef() && exchange.setFinished(status))
   {
      if (status == lbool::False())
      {
         if (!exchange.isOk() && verbosity > -1)
            std::cout << "c solved through import\n";

         drat.addEmptyClause();
         drat.flush();
      }
   } else
      status = lbool::Undef();

   return status;
}

template <typename TemplateConfig>
void Solver<TemplateConfig>::printStats() const
{
   stat.print();
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
