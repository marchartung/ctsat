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
#ifndef SOURCES_CORE_SOLVERCONFIG_H_
#define SOURCES_CORE_SOLVERCONFIG_H_

#include "initial/Inputs.h"
#include "utils/CPUBind.h"

#include <cstdint>
#include <string>
#include <cassert>
#include <array>

namespace ctsat
{

enum class AnalyzeHeuristic
{
   FUIP,
   MUIP,
   LEVELAWARE
};

inline AnalyzeHeuristic getAnalyze()
{
   std::string str(Inputs::analyze);
   if (str == "muip")
      return AnalyzeHeuristic::MUIP;
   else if (str == "laa")
      return AnalyzeHeuristic::LEVELAWARE;
   return AnalyzeHeuristic::FUIP;
}

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
   IMPORTBUFFER
};

inline ExchangeHeuristic getExchange()
{
   std::string str(Inputs::exchange);
   if (str == "none")
      return ExchangeHeuristic::NOEXCHANGE;
   else if (str == "onewatched")
   {
      assert(false && "not implemented yet");
      return ExchangeHeuristic::IMPORTBUFFER;
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

enum class MpiNodeType
{
   SINGLE_THREAD,
   HALF_THREADS,
   FULL_THREADS
};

class SolverConfig
{
 public:

   static MpiNodeType getNodeType(int const rank)
   {
      int idx = rank%10;
      if(idx == 10)
         return MpiNodeType::SINGLE_THREAD ;
      else if(idx%3 == 0)
         return MpiNodeType::FULL_THREADS;
      else
         return MpiNodeType::HALF_THREADS;
   }

   static int countExitingNodesTypes(int const rank)
   {
      MpiNodeType const type = getNodeType(rank);
      switch(type)
      {
         case MpiNodeType::SINGLE_THREAD:
               return rank/10;
         default:
            assert(false);
            // FIXME only supported for single threaded nodes
            return 0;
      }
   }

   static int getNumRecommandedThreads(int const rank)
   {
      NumaAwareSet const & numa = NumaAwareSet::instance;
      int res = Inputs::nThreads;
      if (Inputs::mpiAutoThreads)
      {
         switch (getNodeType(rank))
         {
            case MpiNodeType::SINGLE_THREAD:
               res = 1;
               break;
            case MpiNodeType::FULL_THREADS:
               res = numa.getNumCores()-1; // communication thread
               break;
            default:
               res = numa.getNumCores() / 2;
         }
      }
      return std::max(1,res);
   }

   static int getMaxNumThreads()
   {
      if (Inputs::mpiAutoThreads)
      {
         NumaAwareSet const & numa = NumaAwareSet::instance;
         return numa.getNumCores();
      } else
         return Inputs::nThreads;
   }

   bool initialSolver;
   BranchHeuristic branch;
   RestartHeuristic restart;
   ReduceHeuristic reduce;
   ExchangeHeuristic exchange;
   DatabaseImplementation database;
   PropagateStyle propagateStyle;
   AnalyzeHeuristic analyze;

   // Analyze
   // LAA:
   bool LAA_alwaySwap;
   int LAA_levelDiffEnforce;
   int LAA_numInitialConflicts;
   int LAA_levelQueueSz;

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
   double vsids_max_var_decay;
   double dist_var_decay;

   // Restart
   int lbd_queue;
   int luby_base_factor;
   double luby_inc_factor;

   // Backtrack
   int chrono;
   int confl_to_chrono;

   // Minimization
   bool eliminate;
   bool useVivification;
   bool remove_satisfied;
   int ccmin_mode;
   int maxEntendedBinaryResolutionSz;
   int maxFullImplicationMinLbd;
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
   int secToSwitchHeuristic;

   uint64_t conflict_budget;
   uint64_t propagation_budget;

   // Parallel
   bool pinSolver;
   bool minimize_import_cl;
   bool onlyExportWhenMin;
   int nThreads;
   int max_export_lbd;
   int max_import_lbd;
   int max_export_sz;
   int numConflictsToDelete;

   static SolverConfig getInputConfig()
   {
      return SolverConfig();
   }

   static void makeUnsatTuned(SolverConfig & res, int num)
   {
      res.max_import_lbd = 5;
      if (num > 2)
      {
         res.vsids_var_decay_timer = std::max((num / 2) * 2000 + 5000, 10000);
         num %= 2;
      }
      switch (num)
      {
         case 0:
            res.branch = BranchHeuristic::VSIDS;
            res.reduce = ReduceHeuristic::GLUCOSE;
            res.restart = RestartHeuristic::GLUCOSE;
            break;
         case 1:
            res.branch = BranchHeuristic::VSIDS;
            res.reduce = ReduceHeuristic::CHANSEOKOH;
            res.restart = RestartHeuristic::GLUCOSE;
            break;
         default:
            assert(false);
      }
   }

   static void makeSatTuned(SolverConfig & res, int num)
   {
      res.branch = BranchHeuristic::LRB;
      res.reduce = ReduceHeuristic::CHANSEOKOH;
      res.restart = RestartHeuristic::LUBY;
      res.max_import_lbd = 3;
      if (num == 1)
         res.luby_inc_factor = 80;
      if (num > 1)
         res.luby_inc_factor += 20 * (num - 1);
   }

   static void makeMixed(SolverConfig & res, int num)
   {

      res.max_import_lbd = 4;
      switch (num % 2)
      {
         case 0:
            res.branch = BranchHeuristic::VSIDS;
            res.reduce = ReduceHeuristic::CHANSEOKOH;
            res.restart = RestartHeuristic::LUBY;
            break;
         case 1:
            res.branch = BranchHeuristic::LRB;
            res.reduce = ReduceHeuristic::GLUCOSE;
            res.restart = RestartHeuristic::LUBY;
            break;
         default:
            assert(false);
      }
   }

   static SolverConfig getThreadConfig(int threadId, double const seedAdd = 0)
   {
      SolverConfig res;
      res.initialSolver = false;
      res.verbosity = -1;
      res.rnd_active = true;
      res.rnd_seed += 9919.4853094755497L * (threadId + seedAdd);

      if (threadId > 1)
      {
         res.initialSolver = false;
         res.verbosity = -1;
         threadId -= 2;
         switch (threadId % 3)
         {
            case 0:
               makeUnsatTuned(res, threadId / 3);
               break;
            case 1:
               makeSatTuned(res, threadId / 3);
               break;
            default:
               makeMixed(res, threadId / 3);
         }
      } else
      {
         // the first two threads use the input heuristic except:
         res.max_import_lbd = 4;
         res.initVarPolZero = threadId % 2;
      }
      return res;

   }

   static SolverConfig getMpiThreadConfig(int threadId, int const rank)
   {
      if (Inputs::mpiAutoThreads
         && getNodeType(rank) == MpiNodeType::SINGLE_THREAD)
         threadId = countExitingNodesTypes(rank);

      return getThreadConfig(threadId, NumaAwareSet::instance.getNumCores() * rank);
   }

   SolverConfig()
         : initialSolver(true),
           branch(getBranch()),
           restart(getRestart()),
           reduce(getReduce()),
           exchange(getExchange()),
           database(getDatabase()),
           propagateStyle(getPropagate()),
           analyze(getAnalyze()),

           LAA_alwaySwap(Inputs::LAA_alwaySwap),
           LAA_levelDiffEnforce(Inputs::LAA_levelDiffEnforce),
           LAA_numInitialConflicts(Inputs::LAA_numInitialConflicts),
           LAA_levelQueueSz(Inputs::LAA_levelQueueSz),

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
           vsids_var_decay_timer(Inputs::vsids_var_decay_timer),
           step_size(Inputs::step_size),
           step_size_dec(Inputs::step_size_dec),
           min_step_size(Inputs::min_step_size),
           vsids_var_decay(Inputs::vsids_var_decay),
           vsids_max_var_decay(Inputs::vsids_max_var_decay),
           dist_var_decay(0.6),

           lbd_queue(50),
           luby_base_factor(Inputs::luby_base_factor),
           luby_inc_factor(Inputs::luby_inc_factor),
           chrono(Inputs::chrono),
           confl_to_chrono(Inputs::conf_to_chrono),

           eliminate(true),
           useVivification(Inputs::useVivification),
           remove_satisfied(true),
           ccmin_mode(Inputs::ccmin_mode),
           maxEntendedBinaryResolutionSz(Inputs::maxEntendedBinaryResolutionSz),
           maxFullImplicationMinLbd(Inputs::maxFullImplicationMinLbd),
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
           secToSwitchHeuristic(Inputs::secToSwitchHeuristic),
           conflict_budget(Inputs::conflict_budget),
           propagation_budget(Inputs::propagation_budget),
           pinSolver(Inputs::pinSolver),
           minimize_import_cl(Inputs::minimize_import_cl),
           onlyExportWhenMin(Inputs::onlyExportWhenMin),
           nThreads(Inputs::nThreads),
           max_export_lbd(Inputs::max_export_lbd),
           max_import_lbd(Inputs::max_import_lbd),
           max_export_sz(Inputs::max_export_sz),
           numConflictsToDelete(Inputs::numConflictsToDelete)
   {
   }
};

}

#endif /* SOURCES_CORE_SOLVERCONFIG_H_ */
