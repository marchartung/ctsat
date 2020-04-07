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

#ifndef SOURCES_MPI_MPISOLVERUNNER_H_
#define SOURCES_MPI_MPISOLVERUNNER_H_

#include "parallel/ParallelSolveRunner.h"
#include "mpi/MPIBroadcastConnector.h"
#include "mpi/MPIExportFilter.h"

namespace CTSat
{

class MPISolverRunner : public ParallelSolverRunner<MPIBroadcastConnector>
{
 public:
   typedef ParallelSolverRunner<MPIBroadcastConnector> Super;
   typedef Super::ConnType ConnType;
   typedef Super::TData TData;
   typedef Super::SolverMemory SolverMemory;

   static void printSolverAnnouncement()
   {
      SolverRunner::printSolverAnnouncement();
      printf("c MPI version of CTSat. In most cases, passed arguments will only\n");
      printf("c effect the first two threads of each MPI process.\n");
   }

   static lbool run()
   {
      Timer initTime;
      std::shared_ptr<SolverMemory> mem = std::make_shared<SolverMemory>();
      if(mem->conn.isRoot() && Inputs::verb > 0)
         printSolverAnnouncement();
      if (mem->conn.isRoot())
      {
         if (!sendInstance(mem))
            return lbool::False();
      } else if (!receiveInstance(mem))
         return lbool::False();

      // start solving
      startThreads(mem, mem->conn.getRank());
      mem->conn.freeInstance();
      if (mem->conn.isRoot())
         std::cout << "c init time: " << initTime.getPassedTime() << std::endl;
      runLoop(mem->conn, mem->tdata, mem->inst);
      lbool res = finalizeResult(mem->conn, mem->inst);
      joinThreads(mem->conn, mem->tdata);
      if (mem->conn.isRoot())
         std::cout << "c complete time: " << initTime.getPassedTime() << std::endl;
      return res;
   }

 protected:

   static bool sendInstance(std::shared_ptr<SolverMemory> mem)
   {
      uint64_t nBytes = 0;
      void const * data = nullptr;
      mem->inst = getInstance((*Inputs::argv)[1]);
      if (mem->inst.isOk())
      {
         nBytes = mem->inst.ca.nBytes();
         data = mem->inst.ca.data();
         assert(mem->inst == SatInstance(data, nBytes));
      }
      mem->conn.sendInstance(data, nBytes);
      return mem->inst.isOk();
   }

   static bool receiveInstance(std::shared_ptr<SolverMemory> mem)
   {
      std::tuple<void const *, uint64_t> rData = mem->conn.receiveInstance();
      if (std::get<1>(rData) == 0)
         return false;
      else
      {
         mem->inst = SatInstance(std::get<0>(rData), std::get<1>(rData));
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
   static void runLoop(ConnType & conn, std::vector<TData> & tdata, SatInstance & inst)
   {
      int const verb = (conn.isRoot()) ? Inputs::verb : -1;
      conn.waitInitialize(tdata.size());
      Timer printTimer(Inputs::print_interval), sendTimer(Inputs::mpi_send_interval);
      MPIStatistc h(conn.getRank());
      MPIExportFilter<ClauseAllocator> filter(conn);
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
               filter.readRecvClauses();
               readClauses = true;
            } else if (sendTimer.isOver())
            {
               updateStat(h, tdata);
               sendTimer.reset();
               filter.addLocalClausesToConn();
               conn.send(h);
               readClauses = false;
            }
         }
         else
            conn.sleep();
      }
      if (!conn.getResult().isUndef() || conn.isAborted())
      {
         h.abort = conn.isAborted();
         h.res = conn.getResult();
         updateStat(h, tdata);
         conn.sendResult(h);
         h.res = conn.getResult();
         if (h.res.isTrue())
            conn.sendModel();
      }
      if (verb > 0)
      {
         conn.getAccumulatedHeader().print();
         filter.printState();
      }
   }

   static lbool finalizeResult(ConnType & conn, SatInstance & inst)
   {
      lbool ret = conn.getResult();
      if (conn.isRoot())
      {
         ret = conn.receiveResult(inst.model.size());
         printf(
               ret.isTrue() ? "s SATISFIABLE\n" :
               ret.isFalse() ? "s UNSATISFIABLE\n" : "s UNKNOWN\n");
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
               if (ModelChecker::checkSat(model, (*Inputs::argv)[1]))
                  std::cout << "c SAT solution is correct\n";
               else
                  std::cout << "c Solution is WRONG!!!!\n";
            }
         }
      }
      return ret;
   }
};

}

#endif /* SOURCES_MPI_MPISOLVERUNNER_H_ */
