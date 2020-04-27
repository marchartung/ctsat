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

#ifndef MPI_MPIPARTITIONER_H_
#define MPI_MPIPARTITIONER_H_

#include "initial/Inputs.h"
#include <mpi/mpi.h> // FIXME

namespace ctsat
{

class MPIPartitioner
{
 public:
   MPIPartitioner();

   bool isRoot() const
   {
      return gRank == 0;
   }
   bool isFilter() const
   {
      return filter;
   }

   int getGlobalRank() const
   {
      return gRank;
   }

   MPI_Comm getPartitionComm() const
   {
      return pComm;
   }
   MPI_Comm getFilterComm() const
   {
      assert(filter);
      return fComm;
   }
   bool hasFilterNodes() const
   {
      return hasFilters;
   }

 private:
   bool filter;
   bool hasFilters;
   int gRank;

   MPI_Comm pComm;
   MPI_Comm fComm;

   void setPartition(int const np, int const nranks);
};

MPIPartitioner::MPIPartitioner()
      : filter(false),
        hasFilters(false),
        gRank(-1),
        pComm(MPI_COMM_WORLD),
        fComm(MPI_Comm())
{
   int gRanks;
   int np = Inputs::nMpiPartitions;
   MPI_Comm_size(MPI_COMM_WORLD, &gRanks);
   MPI_Comm_rank(MPI_COMM_WORLD, &gRank);
   // calc partion:
   if (np * 2 >= gRanks)
   {
      np = gRanks / 2;
      std::cout
         << "c Warning: Reduced number of partitions to "
         << np
         << ". The minimum number of nodes per partition is 2"
         << std::endl;
   }
   if(np > 1)
      setPartition(np,gRanks);
}


void MPIPartitioner::setPartition(int const np, int const nranks)
{
   assert(np > 1);
   // create partition comm:
   int const nPerPart = nranks/np;
   assert(nPerPart >= 2);
   int color = gRank/nPerPart;
   if(color >= np)
      color = np-1; // special case caused by rounding, the last one will be bigger
   MPI_Comm_split(MPI_COMM_WORLD,color,gRank,&pComm);
   int pRank;
   MPI_Comm_rank(pComm,&pRank);

   // create filter comm
   filter = pRank == nPerPart-1; // ensure here that gRank == 0 is not filter
   assert(!filter || gRank != 0);
   MPI_Comm_split(MPI_COMM_WORLD,filter,gRank,&fComm);
   hasFilters = true;
}

}

#endif /* MPI_MPIPARTITIONER_H_ */
