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

#ifndef SOURCES_CORE_SOLVERRUNNER_H_
#define SOURCES_CORE_SOLVERRUNNER_H_

#include <errno.h>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/ResourceLimits.h"
#include "utils/Options.h"
#include "initial/Inputs.h"

#include "initial/SatInstance.h"
#include "initial/SolverConfig.h"

#include "initial/Preprocessor.h"
#include "database/BasicTypes.h"

#include "core/ModelChecker.h"

#include "core/Solver.h"

#include "database/MinisatAllocatorDb.h"

#include "branch/DistLrbVsidsBranch.h"
#include "branch/Dist.h"
#include "branch/Vsids.h"
#include "branch/Lrb.h"

#include "reduce/ChanseokOh.h"
#include "reduce/Glucose.h"

#include "restart/GlucoseLuby.h"
#include "restart/Glucose.h"
#include "restart/Luby.h"

#include "core/NoConnector.h"

#include "exchange/NoClauseExchanger.h"
#include "exchange/SimpleClauseExchanger.h"
#include "exchange/ConflictExchange.h"

#include "propagate/MiniSatPropagate.h"

#include "core/Solver_impl.hpp"

namespace CTSat
{

template <typename Connector>
struct ThreadData
{
   ThreadData() = delete;
   ThreadData(ThreadData const &) = delete;
   ThreadData& operator=(ThreadData const &) = delete;
   ThreadData& operator=(ThreadData &&) = delete;

   int const threadId;

   pthread_t pthreadId;
   SatInstance const & inst;
   Connector & conn;
   Statistic * stat;
   SolverConfig config;

   ThreadData(ThreadData && in)
         : threadId(in.threadId),
           pthreadId(in.pthreadId),
           inst(in.inst),
           conn(in.conn),
           stat(in.stat),
           config(in.config)
   {

   }

   ThreadData(SatInstance const& inst, Connector & conn, int const threadId, int const rank = 0)
         : threadId(threadId),
           pthreadId(),
           inst(inst),
           conn(conn),
           stat(nullptr),
           config(SolverConfig::getThreadConfig(threadId, rank))
   {
   }
};

class SolverRunner
{
 public:

   static void printSolverAnnouncement()
   {
      printf("c ################################  CTSat  #############################\n");
      printf("c CTSat is based on MapleLCMChronoBT. So thanks for all the great ideas\n");
      printf("c and implementations (Minisat, Glucose, Maple,...) enabling this solver\n");
   }

   static lbool run(SolverConfig const & config)
   {
      printSolverAnnouncement();
      SatInstance inst;
      return runDatabase<NoConnector>(config, nullptr);
   }

   template <typename Connector>
   static lbool runDatabase(SolverConfig const & config, void * threadData = nullptr)
   {
      switch (config.database)
      {
         case DatabaseImplementation::MINISAT:
            return runBranch<Connector, ClauseAllocator>(config, threadData);
         default:
            assert(false);
            return lbool::Undef();
      }
   }

   template <typename Connector, typename Database>
   static lbool runBranch(SolverConfig const & config, void * threadData = nullptr)
   {
      switch (config.branch)
      {
         case BranchHeuristic::DISTLRBVSIDS:
            return runRestart<Connector, Database, DistLrbVsidsBranch<Database>>(config, threadData);
         case BranchHeuristic::DIST:
            return runRestart<Connector, Database, DistBranch<Database>>(config, threadData);
         case BranchHeuristic::VSIDS:
            return runRestart<Connector, Database, VsidsBranch<Database>>(config, threadData);
         case BranchHeuristic::LRB:
            return runRestart<Connector, Database, LrbBranch<Database>>(config, threadData);
         default:
            assert(false);
            return lbool::Undef();
      }
   }

   template <typename Connector, typename Database, typename Branch>
   static lbool runRestart(SolverConfig const & config, void * threadData = nullptr)
   {
      switch (config.restart)
      {
         case RestartHeuristic::MIXED:
            return runReduce<Connector, Database, Branch, GlucoseLuby>(config, threadData);
         case RestartHeuristic::LUBY:
            return runReduce<Connector, Database, Branch, LubyRestart>(config, threadData);
         case RestartHeuristic::GLUCOSE:
            return runReduce<Connector, Database, Branch, GlucoseRestart>(config, threadData);
         default:
            assert(false);
            return lbool::Undef();

      }
   }

   template <typename Connector, typename Database, typename Branch, typename Restart>
   static lbool runReduce(SolverConfig const & config, void * threadData = nullptr)
   {
      switch (config.reduce)
      {
         case ReduceHeuristic::CHANSEOKOH:
            return runPropagate<Connector, Database, Branch, Restart, ChanseokOhReduce<Database>>(
                  config, threadData);
         case ReduceHeuristic::GLUCOSE:
            return runPropagate<Connector, Database, Branch, Restart, GlucoseReduce<Database>>(
                  config, threadData);
         default:
            assert(false);
            return lbool::Undef();

      }
   }

   template <typename Connector, typename Database, typename Branch, typename Restart,
         typename Reduce>
   static lbool runPropagate(SolverConfig const & config, void * threadData = nullptr)
   {
      switch (config.propagateStyle)
      {
         case PropagateStyle::MINISAT:

            return runExchange<Connector, Database, Branch, Restart, Reduce,
                  MinisatPropagate<Database>>(config, threadData);
         default:
            assert(false);
            return lbool::Undef();

      }

   }

   template <typename Connector, typename Database, typename Branch, typename Restart,
         typename Reduce, typename Propagate>
   static lbool runExchange(SolverConfig const & config, void * threadData = nullptr)
   {
      switch (config.exchange)
      {
         case ExchangeHeuristic::NOEXCHANGE:
            return run_main<Connector, Database, Branch, Restart, Reduce, Propagate,
                  NoClauseExchanger<Database, Connector,Propagate>>(config, threadData);
         case ExchangeHeuristic::SIMPLE:
            return run_main<Connector, Database, Branch, Restart, Reduce, Propagate,
                  SimpleClauseExchanger<Database, Connector,Propagate>>(config, threadData);
         case ExchangeHeuristic::IMPORTBUFFER:
            return run_main<Connector, Database, Branch, Restart, Reduce, Propagate,
                  ConflictExchange<Database, Connector,Propagate>>(config, threadData);
         default:
            assert(false);
            return lbool::Undef();

      }

   }

   template <typename Connector, typename Database, typename Branch, typename Restart,
         typename Reduce, typename Propagate, typename Exchanger>
   static lbool run_main(SolverConfig const & config, void * threadData = nullptr)
   {

      if (threadData == nullptr)
      {
         assert(config.initialSolver);
         return run_sequential<Database, Branch, Restart, Reduce, Propagate>(config);
      } else
      {

         assert(!config.initialSolver);
         return run_parallel<Connector, Database, Branch, Restart, Reduce, Propagate, Exchanger>(
               config, threadData);
      }
   }

   static SatInstance getInstance(std::string const & filename)
   {

      Preprocessor prepro;
      SatInstance res = prepro.getInstance(filename);
//      ModelChecker::printSatisfiedClauses(res.model, filename);
      return res;
   }

   template <typename Database, typename Branch, typename Restart, typename Reduce,
         typename Propagate>
   static lbool run_sequential(SolverConfig const & config)
   {
      lbool ret = lbool::Undef();
      SatInstance inst = getInstance((*Inputs::argv)[1]);

      if (!inst.isOk())
         ret = lbool::False();
      else
      {
         NoConnector connector;
         Solver<
               TemplateConfiguration<NoConnector, Database, Branch, Restart, Reduce, Propagate,
                     NoClauseExchanger<Database, NoConnector,Propagate>>> solver(config, connector,
                                                                       std::move(inst));

         if (solver.isOk())
            ret = solver.solve();
         if (config.verbosity > 0)
         {
            solver.printStats();
            printf("\n");
         }
         if (ret.isTrue())
         {

            EliminatedClauseDatabase & elimDb = inst.elimDb;
            vec<lbool> & model = inst.model;
            solver.setModel(model);
            elimDb.extendModel(model);

            if (Inputs::model)
               elimDb.printModel(model);
            if (Inputs::verifySat)
            {
               if (ModelChecker::checkSat(model, (*Inputs::argv)[1]))
                  std::cout << "c SAT solution is correct\n";
               else
                  std::cout << "c Solution is WRONG!!!!\n";
            }
         }
      }
      printf(
            ret.isTrue() ? "s SATISFIABLE\n" : ret.isFalse() ? "s UNSATISFIABLE\n" : "s UNKNOWN\n");

      return ret;
   }

   template <typename Connector, typename Database, typename Branch, typename Restart,
         typename Reduce, typename Propagate, typename Exchanger>
   static lbool run_parallel(SolverConfig const & config, void * threadData)
   {
      assert(!config.drup);
      lbool ret = lbool::Undef();

      ThreadData<Connector> & data = *reinterpret_cast<ThreadData<Connector>*>(threadData);

      data.conn.notifyThreadStart();
      Solver<
            TemplateConfiguration<Connector, Database, Branch, Restart, Reduce, Propagate, Exchanger>> solver(
            config, data.conn, data.inst);
      data.stat = &solver.getStatistic();
      data.conn.notifyThreadInitialized();

      if (solver.isOk())
         ret = solver.solve();

      if (!ret.isUndef())
      {
         if (ret.isTrue())
         {
            typedef typename Database::lbool bool_type;
            vec<bool_type> model(solver.nVars(), bool_type::Undef());
            solver.setModel(model);
            data.conn.commitModel(std::move(model));
         }
      }

      data.conn.notifyThreadEnd();
      data.stat = nullptr;
      return ret;
   }

};
}

#endif /* SOURCES_CORE_SOLVERRUNNER_H_ */
