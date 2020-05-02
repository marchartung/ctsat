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

#ifndef SOURCES_PARALLEL_PARALLELSOLVERUNNER_H_
#define SOURCES_PARALLEL_PARALLELSOLVERUNNER_H_

#include "core/SolverRunner.h"
#include "parallel/ParallelStatistic.h"
#include "core/Statistic.h"
#include "parallel/AtomicConnector.h"
#include "utils/Timer.h"
#include "utils/Logger.h"

#include <pthread.h>
#include <vector>
#include <memory>

namespace ctsat
{

template <typename Connector>
class ParallelSolverRunner : public SolverRunner
{
 public:
   typedef Connector ConnType;
   typedef ThreadData<ConnType> TData;

   static void printSolverAnnouncement()
   {
      SolverRunner::printSolverAnnouncement();
      printf("c Parallel version of CTSat. In most cases, passed arguments will only\n");
      printf("c effect the first two threads.\n");
   }
   static lbool run()
   {
      Timer initTime;
      printSolverAnnouncement();
      std::shared_ptr<SolverMemory> mem = getSolveMemory();
      startThreads(Inputs::nThreads, mem->tdata, mem->conn, mem->inst);
      runLoop(mem->conn, mem->tdata, mem->inst);
      lbool res = finalizeResult(mem->conn, mem->inst);
      joinThreads(mem->conn, mem->tdata);
      LOG_DEINIT
      std::cout << "c complete time: " << initTime.getPassedTime() << std::endl;
      return res;
   }

 protected:

   struct SolverMemory
   {
      AtomicConnector conn;
      SatInstance inst;
      std::vector<TData> tdata;

      SolverMemory()
            : conn(Inputs::mbExchangeBufferPerThread * 1024ul * 1024ul, Inputs::nThreads)
      {
      }
   };

   static void printStats(std::vector<TData> & tdata)
   {
      PStatistic stat;
      for (size_t i = 0; i < tdata.size(); ++i)
         stat.add(*tdata[i].stat);
      stat.print();
   }

   static std::shared_ptr<SolverMemory> getSolveMemory()
   {
      LOG_INIT(0)
      std::shared_ptr<SolverMemory> res = std::make_shared<SolverMemory>();
      res->inst = getInstance((*Inputs::argv)[1], SolverConfig::getInputConfig());
      return res;
   }

   static void startThreads(
                            size_t const threadCount,
                            std::vector<TData> & tdata,
                            ConnType & conn,
                            SatInstance const & inst,
                            int const rank = 0)
   {
      LOG("Starting " + std::to_string(threadCount) + " threads")
      tdata.reserve(threadCount);

      for (size_t i = 0; i < threadCount; ++i)
      {
         tdata.emplace_back(TData(inst, conn, i, rank));
         if (!startThread(tdata[i]))
            throw InputException();  // FIXME stop started threads
      }
   }

   static lbool finalizeResult(ConnType & conn, SatInstance & inst)
   {
      lbool ret = conn.getResult();

      printf(
            ret.isTrue() ? "s SATISFIABLE\n" : ret.isFalse() ? "s UNSATISFIABLE\n" : "s UNKNOWN\n");
      if (ret.isTrue())
      {
         EliminatedClauseDatabase const & elimDb = inst.elimDb;
         vec<lbool> & model = conn.getModel();
         for (int i = 0; i < inst.model.size(); ++i)
            model[i] = (inst.model[i].isUndef()) ? model[i] : inst.model[i];
         elimDb.extendModel(model);
         if (Inputs::model)
            elimDb.printModel(model);
         if (Inputs::verifySat)
         {
            if (ModelChecker::checkSat(model, (*Inputs::argv)[1],inst.isDecisionVar))
               std::cout << "c SAT solution is correct\n";
            else
               std::cout << "c Solution is WRONG!!!!\n";
         }
      }
      return ret;
   }

   static void runLoop(ConnType & conn, std::vector<TData> & tdata, SatInstance & inst)
   {
      int const verb = Inputs::verb;
      conn.waitInitialize(tdata.size());
      inst.ca.clear(true);
      Timer printTimer(Inputs::print_interval);
      while (conn.nRunningThreads() > 0 && !conn.isFinished())
      {
         conn.sleep();
         if (verb > 0 && printTimer.isOver())
         {
            printStats(tdata);
            printTimer.reset();
         }
      }

      printStats(tdata);
   }

   static void joinThreads(ConnType & conn, std::vector<TData> & tdata)
   {
      conn.allowMemFree();
      for (size_t i = 0; i < tdata.size(); ++i)
      {
         if (pthread_join(tdata[i].pthreadId, NULL) != 0)
            std::cout
               << "c Warning: Thread join failed on thread "
               << tdata[i].threadId
               << "("
               << tdata[i].pthreadId
               << ")"
               << std::endl;
      }
   }

   static void * pthreadStart(void * v)
   {
      ThreadData<ConnType> & data = *reinterpret_cast<ThreadData<ConnType>*>(v);
      LOG_INIT(data.rank)
      LOG("Thread " + std::to_string(pthread_self()) + " started")
      SolverRunner::runDatabase<ConnType>(SolverConfig::getThreadConfig(data.threadId), v);
      LOG("Thread finished")
      LOG_DEINIT
      pthread_exit(NULL);
   }

   static bool startThread(TData & d)
   {
      void * data = &d;
      int res = pthread_create(&d.pthreadId, NULL, ParallelSolverRunner::pthreadStart, data);
      LOG("Started thread " + std::to_string(d.pthreadId));
      return res == 0;
   }
};

}

#endif /* SOURCES_PARALLEL_PARALLELSOLVERUNNER_H_ */
