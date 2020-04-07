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

#ifndef SOURCES_MPI_MPIEXPORTFILTER_H_
#define SOURCES_MPI_MPIEXPORTFILTER_H_

#include "mpi/MPIBroadcastConnector.h"
#include "exchange/SimpleClauseExchanger.h"
#include "parallel/AtomicConnector.h"
#include "mtl/Alloc.h"

namespace CTSat
{
template <typename Database>
class MPIExportFilter
{
   typedef ExportClause<Database> EClause;
   typedef typename Database::base_type base_type;
   typedef RegionAllocator<base_type> Allocator;
   typedef typename Database::Ref Ref;
 public:

   MPIExportFilter(MPIBroadcastConnector & conn);

   void readRecvClauses();
   void updateLocalClauses();
   void addLocalClausesToConn();

   void printState() const;
 private:
   MPICommunicationStat statistic;
   Allocator db;
   vec<vec<Ref>> lbdRefs;
   AtomicConnector::size_type readPos;
   MPIBroadcastConnector & conn;
   AtomicRingAllocator & importExchange;
   AtomicRingAllocator & exportExchange;

   void clear()
   {
      db.clear();
      for (int i = 0; i < lbdRefs.size(); ++i)
         lbdRefs[i].clear();
   }

   EClause const & get(Ref const ref) const
   {
      return *reinterpret_cast<EClause const *>(db.lea(ref));
   }

   Ref alloc(EClause const & c)
   {
      assert(EClause::nbytes(c) % sizeof(base_type) == 0);
      Ref const res = db.alloc(EClause::nbytes(c) / sizeof(base_type));
      void * pos = db.lea(res);
      if (pos != new (pos) EClause(c))
         assert(false && "The current value_type is not the systems alignment");
      return res;
   }
};

template <typename Database>
inline MPIExportFilter<Database>::MPIExportFilter(MPIBroadcastConnector& conn)
      : statistic(conn.getNumRanks(), conn.getSendBufferSize()),
        db(1024),
        readPos(0),
        conn(conn),
        importExchange(conn.getMPIExchangeBuffer()),
        exportExchange(conn.getLocalExchangeBuffer())
{
}

template <typename Database>
inline void MPIExportFilter<Database>::readRecvClauses()
{
   for (int i = 0; i < conn.getNumRanks(); ++i)
   {
      if (i == conn.getRank())
         continue;
      uint64_t sumBytes = 0, nCl = 0, sumLbd = 0;
      uint8_t const * p = conn.getRecvBuffer(i);
      while (*reinterpret_cast<base_type const*>(p) != std::numeric_limits<base_type>::max())
      {
         EClause const & e = *reinterpret_cast<EClause const*>(p);
         sumLbd += e.lbd;
         assert(e.size() <= 30 || e.lbd < 3);
         uint64_t const nBytes = EClause::nbytes(e);
         assert(nBytes % sizeof(base_type) == 0);
         importExchange.allocConstruct<EClause, EClause const &>(nBytes, e);
         p += nBytes;
         sumBytes += nBytes;
         ++nCl;
      }
      statistic.addComm(false, sumBytes, nCl, sumLbd);
   }

}

template <typename Database>
inline void MPIExportFilter<Database>::updateLocalClauses()
{
   while (exportExchange.isValid(readPos))
   {
      EClause const & e = exportExchange.get<EClause>(readPos);
      if (lbdRefs.size() <= e.lbd)
         lbdRefs.growTo(e.lbd + 1);
      lbdRefs[e.lbd].push(alloc(e));
      readPos = exportExchange.getNextPos(readPos);
   }
}

template <typename Database>
inline void MPIExportFilter<Database>::printState() const
{
   statistic.print();
}

template <typename Database>
inline void MPIExportFilter<Database>::addLocalClausesToConn()
{
   bool bufferFull = false;
   uint64_t sumBytes = sizeof(base_type), nCl = 0, sumLbd = 0;  // for marker
   for (int i = 0; i < lbdRefs.size(); ++i)
   {
      for (int j = 0; j < lbdRefs[i].size(); ++j)
      {
         EClause const & e = get(lbdRefs[i][j]);
         assert(e.size() <= 30 || e.lbd < 3);
         sumLbd += e.lbd;
         uint64_t const nBytes = EClause::nbytes(e);
         if (conn.nFreeSendBytes() >= nBytes + sizeof(base_type))
         {
            if (!conn.addToSend(&e, nBytes))
               assert(false);
         } else
         {
            bufferFull = true;
            break;
         }
         sumBytes += nBytes;
         ++nCl;
      }
      if (bufferFull)
         break;
   }

   // mark end:
   base_type b = std::numeric_limits<base_type>::max();
   if (!conn.addToSend(&b, sizeof(base_type)))
      assert(false);
   clear();
   statistic.addComm(true, sumBytes, nCl, sumLbd);
}
}

#endif /* SOURCES_MPI_MPIEXPORTFILTER_H_ */
