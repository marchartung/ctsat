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

#ifndef SOURCES_MPI_MPISOLVERUNNER_H_
#define SOURCES_MPI_MPISOLVERUNNER_H_

#include "core/SolverRunner.h"
#include "mpi/MPIBroadcastConnector.h"
#include "mpi/MPIExportFilter.h"
#include "mpi/MPISatDistributor.h"
#include "mpi/MPIPartitioner.h"
#include "utils/Logger.h"

#include <memory>

namespace ctsat
{

class MPISolverRunner : public SolverRunner
{
 public:
   typedef MPIBroadcastConnector ConnType;
   typedef ThreadData<ConnType> TData;

   struct MPIMemory
   {
      MPIPartitioner partitioner;
      int nthreads;
      MPISatDistributor distributor;
      MPIBroadcastConnector conn;
      SatInstance inst;
      std::vector<TData> tdata;
      MPIMemory()
            : partitioner(),
              nthreads(SolverConfig::getNumRecommandedThreads(partitioner.getGlobalRank())),
              distributor(),
              conn((partitioner.isFilter()) ?
                    0 : Inputs::mbExchangeBufferPerThread * 1024ul * 1024ul,
                   Inputs::mpiMbBufferSize * 1024ul * 1024ul, nthreads + 1,  // +1 one because mpi thread reads too
                   SolverConfig::getMaxNumThreads(), partitioner.getPartitionComm())
      {
      }
   };

   static void printSolverAnnouncement()
   {
      SolverRunner::printSolverAnnouncement();
      printf("c MPI version of CTSat. In most cases, passed arguments will only\n");
      printf("c effect the first two threads of each MPI process.\n");
   }

   static lbool run()
   {
      Timer initTime;
      std::shared_ptr<MPIMemory> mem = init();
      if (mem->partitioner.isRoot())
         std::cout << "c init time: " << initTime.getPassedTime() << std::endl;
      // start solving
      lbool res = lbool::Undef();
      if (mem->tdata.empty())
      {
         LOG("Solved through preprocessor")
         res = lbool::False();
      } else
      {
         if (!mem->partitioner.isFilter())
            runLoop(*mem);
         else
            runLoopFilter(*mem);
      }
      res = finalizeResult(*mem);
      if (mem->partitioner.isRoot())
         std::cout << "c complete time: " << initTime.getPassedTime() << std::endl;
      deinit(mem);
      return res;
   }

 protected:

   static void deinit(std::shared_ptr<MPIMemory> & mem)
   {
      LOG("Joining threads")
      joinThreads(mem->conn, mem->tdata);
      LOG("Clear memory")
      mem.reset();
      LOG("Finalize mpi")
      MPI_Finalize();
      LOG_DEINIT;
   }

   static std::shared_ptr<MPIMemory> init()
   {
      MPI_Init(Inputs::argc, Inputs::argv);
      std::shared_ptr<MPIMemory> mem = std::make_shared<MPIMemory>();

      LOG_INIT(mem->partitioner.getGlobalRank());
      if (mem->partitioner.isRoot() && Inputs::verb > 0)
         printSolverAnnouncement();
      if (mem->partitioner.isRoot())
      {
         if (!sendInstance(*mem))
            return mem;
      } else if (!receiveInstance(*mem))
         return mem;
      if (!mem->partitioner.isFilter())
         startThreads(SolverConfig::getNumRecommandedThreads(mem->partitioner.getGlobalRank()),
                      mem->tdata, mem->conn, mem->inst, mem->partitioner.getGlobalRank());
      else
         startFilterThreads(SolverConfig::getMaxNumThreads(), mem->tdata, mem->conn, mem->inst,
                            mem->partitioner.getGlobalRank());
      mem->distributor.freeInstance();
      mem->conn.waitInitialize(mem->tdata.size());
      if (!Inputs::verifySat)
         mem->inst.ca.clear(true);
      return mem;
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

   static void joinThreads(ConnType & conn, std::vector<TData> & tdata)
   {
      conn.allowMemFree();
      for (size_t i = 0; i < tdata.size(); ++i)
      {
         if (pthread_join(tdata[i].pthreadId, NULL) != 0)
            std::cout << "c Warning: Thread join failed on thread " << tdata[i].threadId << "("
                      << tdata[i].pthreadId << ")" << std::endl;
      }
   }

   static void startFilterThreads(
                                  size_t const threadCount,
                                  std::vector<TData> & tdata,
                                  ConnType & conn,
                                  SatInstance const & inst,
                                  int const rank = 0)
   {
      LOG("Starting " + std::to_string(threadCount) + " filter threads")
      tdata.reserve(threadCount);

      for (size_t i = 0; i < threadCount; ++i)
      {
         tdata.emplace_back(TData(inst, conn, i, rank));
         if (!startFilterThread(tdata[i]))
            throw InputException();  // FIXME stop started threads
      }
   }

   static void * pthreadStartFilter(void * v)
   {
//      TData & data = *reinterpret_cast<TData*>(v);
      LOG_INIT(data.rank)
LOG("Filter thread " + std::to_string(pthread_self()) + " started")
            assert(false);
      LOG("Thread finished")
      LOG_DEINIT
      pthread_exit(NULL);
   }

   static bool startFilterThread(TData & d)
   {
      void * data = &d;
      int res = pthread_create(&d.pthreadId, NULL, MPISolverRunner::pthreadStartFilter, data);
      LOG("Started thread " + std::to_string(d.pthreadId));
      return res == 0;
   }

   static void * pthreadStart(void * v)
   {
      TData & data = *reinterpret_cast<TData*>(v);
      LOG_INIT(data.rank)
      LOG("Thread " + std::to_string(pthread_self()) + " started")
      SolverRunner::runDatabase<ConnType>(
            SolverConfig::getMpiThreadConfig(data.threadId, data.rank), v);
      LOG("Thread finished")
      LOG_DEINIT
      pthread_exit(NULL);
   }

   static bool startThread(TData & d)
   {
      void * data = &d;
      int res = pthread_create(&d.pthreadId, NULL, MPISolverRunner::pthreadStart, data);
      LOG("Started thread " + std::to_string(d.pthreadId));
      return res == 0;
   }

   static bool sendInstance(MPIMemory& mem)
   {
      LOG("Sending instance")
      uint64_t nBytes = 0;
      void const * data = nullptr;
      mem.inst = getInstance((*Inputs::argv)[1], SolverConfig::getInputConfig());
      if (mem.inst.isOk())
      {
         nBytes = mem.inst.ca.nBytes();
         data = mem.inst.ca.data();
         assert(mem.inst == SatInstance(data, nBytes));
      }
      mem.distributor.sendInstance(data, nBytes);
      return mem.inst.isOk();
   }

   static bool receiveInstance(MPIMemory& mem)
   {
      LOG("Receive instance")
      std::tuple<void const *, uint64_t> rData = mem.distributor.receiveInstance();
      if (std::get<1>(rData) == 0)
         return false;
      else
      {
         mem.inst = SatInstance(std::get<0>(rData), std::get<1>(rData));
         return true;
      }
   }

   static void updateStat(MPIStatistc & header, std::vector<TData> const & tdata)
   {
      PStatistic & stat = *dynamic_cast<PStatistic*>(&header);
      stat.reset();
      for (size_t i = 0; i < tdata.size(); ++i)
         stat.add(*tdata[i].stat);
   }
   static void runLoop(MPIMemory& mem)
   {
      LOG("Starts main loop")
      MPIBroadcastConnector & conn = mem.conn;
      int const verb = (mem.partitioner.isRoot()) ? Inputs::verb : -1;
      Timer printTimer(Inputs::print_interval), sendTimer(Inputs::mpi_send_interval);
      MPIStatistc h;
      MPIExportFilter<ClauseAllocator> filter(
            conn, Inputs::mpiHashClauseFilter && !mem.partitioner.hasFilterNodes());
      bool readClauses = true;

      while (conn.nRunningThreads() > 0 && !conn.isFinished())
      {
         if (verb > 0 && printTimer.isOver())
         {
            printTimer.reset();
            conn.getAccumulatedHeader().print();
            filter.printState();
         }
         filter.updateLocalClauses();
         conn.progress();
         if (conn.isAllowedToSend())
         {
            if (!readClauses)
            {
               LOG("Reads received clauses")
               filter.readRecvClauses();
               readClauses = true;
            } else if (sendTimer.isOver())
            {
               updateStat(h, mem.tdata);
               sendTimer.reset();
               filter.addLocalClausesToConn();
               LOG("Sends clauses")
               conn.send(h);
               readClauses = false;
            }
         } else
            conn.sleep();
      }
      LOG("Exit main loop")
      updateStat(h, mem.tdata);
      LOG("Sends abort")
      conn.abort(h);
   }

   static void runLoopFilter(MPIMemory& mem)
   {
      LOG("Starts filter loop")
      MPIBroadcastConnector & conn = mem.conn;
      assert(!mem.partitioner.isRoot());
      Timer printTimer(Inputs::print_interval), sendTimer(Inputs::mpi_send_interval);
      MPIStatistc h;

      bool readClauses = true;

      while (conn.nRunningThreads() > 0 && !conn.isFinished())
      {
         conn.progress();
         if (conn.isAllowedToSend())
         {
            if (!readClauses)
            {
               LOG("Reads received clauses")
               readClauses = true;
            } else if (sendTimer.isOver())
            {
               updateStat(h, mem.tdata);
               sendTimer.reset();
               LOG("Sends clauses")
               conn.send(h);
               readClauses = false;
            }
         } else
            conn.sleep();
      }
      LOG("Exit main loop")
      updateStat(h, mem.tdata);
      LOG("Sends abort")
      conn.abort(h);
   }
   static lbool finalizeResult(MPIMemory & mem)
   {
      LOG("Finalize result")
      bool isWinner = false;
      lbool ret = mem.conn.getResult();
      if (!ret.isUndef() && mem.distributor.trySetRootResult(ret))
      {
         isWinner = true;

         LOG("Won the race")
         if (!mem.partitioner.isRoot() && ret.isTrue())
         {
            LOG("Sending model")
            mem.distributor.sendModel(mem.conn.getModel());
         }
      }
      if (mem.partitioner.isRoot())
      {
         if (!isWinner)
         {
            LOG("Receives result")
            ret = mem.distributor.getRootResult();
         }

         if (ret.isTrue())
         {

            vec<lbool> model;
            if (!isWinner)
            {
               LOG("Receives model")
               model.growTo(mem.inst.model.size(), lbool::Undef());
               mem.distributor.receiveModel(model);
            } else
               mem.conn.getModel().swap(model);
            if (ret.isTrue())
            {
               LOG("Sets up model")
               std::cout << "c Winner rank: " << mem.distributor.getWinner() << std::endl;
               SatInstance & inst = mem.inst;
               EliminatedClauseDatabase const & elimDb = inst.elimDb;
               for (int i = 0; i < mem.inst.model.size(); ++i)
               {
                  assert(!inst.isDecisionVar[i] || !model[i].isUndef());
                  inst.model[i] = (inst.isDecisionVar[i]) ? model[i] : inst.model[i];
               }
               if (Inputs::verifySat
                  && !ModelChecker::checkSat(inst.model, (*Inputs::argv)[1], inst.isDecisionVar,
                                             true))
               {
                  std::cout << "c SAT solution is before extend wrong\n";
                  if(!inst.checkModel())
                     std::cout << "c created instance was already wrong\n";
                  else
                     std::cout << "Created instance was correct. Is there a bug in the preprocessor?\n";
                  assert(false);
               }

               elimDb.extendModel(inst.model);
               if (Inputs::model)
                  elimDb.printModel(inst.model);
               if (Inputs::verifySat)
               {
                  if (ModelChecker::checkSat(inst.model, (*Inputs::argv)[1], inst.isDecisionVar))
                     std::cout << "c SAT solution is correct\n";
                  else
                  {
                     std::cout << "c Solution is WRONG!!!!\n";
                     assert(false);
                  }
               }
            }
            printf(
                  ret.isTrue() ? "s SATISFIABLE\n" :
                  ret.isFalse() ? "s UNSATISFIABLE\n" : "s UNKNOWN\n");
         }
      }
      return ret;
   }
};

}

#endif /* SOURCES_MPI_MPISOLVERUNNER_H_ */
