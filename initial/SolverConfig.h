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
#ifndef SOURCES_CORE_SOLVERCONFIG_H_
#define SOURCES_CORE_SOLVERCONFIG_H_

#include "initial/Inputs.h"
#include <cstdint>
#include <string>

namespace CTSat
{
enum class BranchHeuristic
{
   DIST,
   LRB,
   VSIDS,
   DISTLRBVSIDS
};

inline BranchHeuristic getBranch()
{
   std::string str(Inputs::branch);
   if (str == "dist_mixed")
      return BranchHeuristic::DISTLRBVSIDS;
   else if (str == "dist")
      return BranchHeuristic::DIST;
   else if (str == "lrb")
      return BranchHeuristic::LRB;
   else if (str == "vsids")
      return BranchHeuristic::VSIDS;
   assert(false);  // FIXME throw proper exception
   return BranchHeuristic::DISTLRBVSIDS;
}

enum class ReduceHeuristic
{
   CHANSEOKOH,
   GLUCOSE
};
inline ReduceHeuristic getReduce()
{
   std::string str(Inputs::reduce);
   if (str == "chanseok")
      return ReduceHeuristic::CHANSEOKOH;
   else if (str == "glucose")
      return ReduceHeuristic::GLUCOSE;
   assert(false);  // FIXME throw proper exception
   return ReduceHeuristic::CHANSEOKOH;
}

enum class RestartHeuristic
{
   LUBY,
   GLUCOSE,
   MIXED
};
inline RestartHeuristic getRestart()
{
   std::string str(Inputs::restart);
   if (str == "mixed")
      return RestartHeuristic::MIXED;
   else if (str == "glucose")
      return RestartHeuristic::GLUCOSE;
   else if (str == "luby")
      return RestartHeuristic::LUBY;
   assert(false);  // FIXME throw proper exception
   return RestartHeuristic::MIXED;
}

enum class ExchangeHeuristic
{
   NOEXCHANGE,
   SIMPLE,
   IMPORTBUFFER
};

inline ExchangeHeuristic getExchange()
{
   std::string str(Inputs::exchange);
   if (str == "none")
      return ExchangeHeuristic::NOEXCHANGE;
   else if (str == "simple")
      return ExchangeHeuristic::SIMPLE;
   else if (str == "onewatched")
   {
      assert(false);
      return ExchangeHeuristic::SIMPLE;
   } else if (str == "importbuff")
      return ExchangeHeuristic::IMPORTBUFFER;
   assert(false);  // FIXME throw proper exception
   return ExchangeHeuristic::NOEXCHANGE;
}

enum class DatabaseImplementation
{
   MINISAT,
   STICKY
};
inline DatabaseImplementation getDatabase()
{
   std::string str(Inputs::database);
   if (str == "minisat")
      return DatabaseImplementation::MINISAT;
   else if (str == "sticky")
      return DatabaseImplementation::STICKY;
   assert(false);  // FIXME throw proper exception
   return DatabaseImplementation::MINISAT;
}

enum class PropagateStyle
{
   MINISAT,
   STICKY
};
inline PropagateStyle getPropagate()
{
   std::string str(Inputs::database);  // still only depended from database, future implement cadical i.e. lingeling watcher
   if (str == "minisat")
      return PropagateStyle::MINISAT;
   else if (str == "sticky")
      return PropagateStyle::STICKY;
   assert(false);  // FIXME throw proper exception
   return PropagateStyle::MINISAT;
}

class SolverConfig
{
 public:
   bool initialSolver;
   BranchHeuristic branch;
   RestartHeuristic restart;
   ReduceHeuristic reduce;
   ExchangeHeuristic exchange;
   DatabaseImplementation database;
   PropagateStyle propagateStyle;

   // Reduce
   // glucose
   int firstReduceDB;
   int incReduceDB;
   int specialIncReduceDB;
   unsigned int maxProtectableLbd;
   //chanseok oh
   int core_lbd_cut;
   double clause_decay;
   uint64_t next_T2_reduce;
   uint64_t next_L_reduce;

   // branch
   bool rnd_active;
   bool rnd_polarity;
   bool initVarPolZero;
   int vsids_var_decay_timer;
   double step_size;
   double step_size_dec;
   double min_step_size;
   double vsids_var_decay;
   double dist_var_decay;
   double timeToBranchSwitch;

   // Restart
   int lbd_queue;
   int restart_first;
   double restart_inc;

   // Backtrack
   int chrono;
   int confl_to_chrono;

   // Minimization
   bool eliminate;
   bool remove_satisfied;
   int ccmin_mode;
   double learntsize_factor;
   double learntsize_inc;
   int incSimplify;
   uint64_t nbconfbeforesimplify;

   // Preprocessor
   bool elim;
   int grow;
   int clause_lim;
   int subsumption_lim;
   double simp_garbage_frac;

   // Divers configs
   double garbage_frac;
   double rnd_seed;

   bool drup;
   std::string drup_file;
   int verbosity;

   uint64_t conflict_budget;
   uint64_t propagation_budget;

   // Parallel
   bool minimize_import_cl;
   int nThreads;
   int max_export_lbd;
   int max_export_sz;
   int numConflictsToDelete;

   static SolverConfig getInputConfig()
   {
      return SolverConfig();
   }

   static SolverConfig getThreadConfig(int const threadId, int const rank = 0)
   {
      SolverConfig res;

      // 1. reduce: chanseok, glucose
      // 2. restart: mixed, luby glucose
      // 3. branch: mixed, vsids, lrb, dist, chb

      //FIXME meaningful combis:

      // balanced:
      // 1. chanseok 2. mixed 3. mixed false init
      // 1. chanseok 2. mixed 3. mixed true init
      // 1. chanseok 2. luby 3. dist-vsids
      // 1- glucose 2. glucose 3. lrb

      //unsat:
      // 1. glucose 2.glucose 3. dist-vsids
      // 1. chanseok 2.glucose 3. dist-vsids

      // sat:
      // 1. chanseok 2. luby 3. lrb
      // 1. chanseok 2. luby 3 chb


      // first: just modify branch init
      res.initVarPolZero = threadId % 2;
      if (rank == 0 && threadId < 2)
      {
         res.rnd_active = false;
         res.rnd_polarity = false;
      } else
      {
         int const globId = res.nThreads*rank+threadId;
         res.rnd_active = true;
         res.rnd_polarity = globId%4 > 1;
         res.rnd_seed = res.rnd_seed+globId;
      }
      res.initialSolver = false;
      res.verbosity = -1;
      return res;
   }

   SolverConfig()
         : initialSolver(true),
           branch(getBranch()),
           restart(getRestart()),
           reduce(getReduce()),
           exchange(getExchange()),
           database(getDatabase()),
           propagateStyle(getPropagate()),

           firstReduceDB(Inputs::first_reduce_db),
           incReduceDB(Inputs::inc_reduce_db),
           specialIncReduceDB(Inputs::spec_inc_reduce_db),
           maxProtectableLbd(Inputs::maxProtectableLbd),

           core_lbd_cut(3),
           clause_decay(Inputs::clause_decay),
           next_T2_reduce(10000),
           next_L_reduce(15000),

           rnd_active(Inputs::rnd_init_act),
           rnd_polarity(Inputs::rnd_polarity),
           initVarPolZero(Inputs::initVarPolZero),
           vsids_var_decay_timer(5000),
           step_size(Inputs::step_size),
           step_size_dec(Inputs::step_size_dec),
           min_step_size(Inputs::min_step_size),
           vsids_var_decay(Inputs::var_decay),
           dist_var_decay(0.6),
           timeToBranchSwitch(2500),

           lbd_queue(50),
           restart_first(Inputs::restart_first),
           restart_inc(Inputs::restart_inc),
           chrono(Inputs::chrono),
           confl_to_chrono(Inputs::conf_to_chrono),
           eliminate(true),
           remove_satisfied(true),
           ccmin_mode(Inputs::ccmin_mode),
           learntsize_factor(1.0 / 3.0),
           learntsize_inc(1.1),
           incSimplify(1000),
           nbconfbeforesimplify(1000),

           elim(Inputs::use_elim),
           grow(Inputs::grow),
           clause_lim(Inputs::clause_lim),
           subsumption_lim(Inputs::subsumption_lim),
           simp_garbage_frac(Inputs::simp_garbage_frac),

           garbage_frac(Inputs::garbage_frac),
           rnd_seed(Inputs::random_seed),
           drup(static_cast<bool>(Inputs::drup)),
           drup_file(static_cast<std::string>(Inputs::drup_file)),
           verbosity(Inputs::verb),
           conflict_budget(Inputs::conflict_budget),
           propagation_budget(Inputs::propagation_budget),
           minimize_import_cl(Inputs::minimize_import_cl),
           nThreads(Inputs::nThreads),
           max_export_lbd(Inputs::max_export_lbd),
           max_export_sz(Inputs::max_export_sz),
           numConflictsToDelete(Inputs::numConflictsToDelete)
   {
   }
};

}

#endif /* SOURCES_CORE_SOLVERCONFIG_H_ */
