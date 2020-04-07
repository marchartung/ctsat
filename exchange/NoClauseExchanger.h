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

#ifndef SOURCES_PARALLEL_NOCLAUSEEXCHANGER_H_
#define SOURCES_PARALLEL_NOCLAUSEEXCHANGER_H_

#include <tuple>
#include <atomic>
#include "core/ImplicationGraph.h"
#include "initial/SolverConfig.h"
#include "mtl/Vec.h"
#include "database/BasicTypes.h"

namespace CTSat
{
template <typename Database, typename Connector, typename PropEngine>
class NoClauseExchanger
{
   typedef typename Database::Clause Clause;
   typedef typename Database::Lit Lit;
   typedef typename Database::CRef CRef;
 public:

   NoClauseExchanger(
                     SolverConfig const & config,
                     Statistic & stat,
                     Database & db,
                     ImplicationGraph<Database> & ig,
                     Connector & conn,
                     PropEngine & propEngine);
   ~NoClauseExchanger();

   bool setFinished(lbool const res);
   bool isFinished() const;

   bool shouldFetch() const
   {
      return false;
   }

   void clauseUsedInConflict(CRef const ref);
   void clauseLearnt(CRef const ref);
   void clauseImproved(CRef const ref);
   void unitLearnt(Lit const l);

   bool hasImportClauses() const;
   std::tuple<bool,CRef> getImportClause();
   Lit getImportUnit();

   bool isOk() const
   {
      return ok;
   }

   void fetchClauses()
   {

   }

   void relocAll(Database& to)
   {

   }

 protected:
   bool ok;
   Statistic & stat;
   Connector & conn;
};

template <typename Database, typename Connector, typename PropEngine>
typename NoClauseExchanger<Database, Connector, PropEngine>::Lit NoClauseExchanger<Database,
      Connector, PropEngine>::getImportUnit()
{
   Lit res = Lit::Undef();
   return res;
}

template <typename Database, typename Connector, typename PropEngine>
inline bool NoClauseExchanger<Database, Connector, PropEngine>::hasImportClauses() const
{
   return false;
}

template <typename Database, typename Connector, typename PropEngine>
std::tuple<bool,typename NoClauseExchanger<Database, Connector, PropEngine>::CRef> NoClauseExchanger<Database,
      Connector, PropEngine>::getImportClause()
{
   return std::make_tuple(true,Database::npos());
}

template <typename Database, typename Connector, typename PropEngine>
inline NoClauseExchanger<Database, Connector, PropEngine>::NoClauseExchanger(
                                                                             SolverConfig const & config,
                                                                             Statistic & stat,
                                                                             Database & db,
                                                                             ImplicationGraph<
                                                                                   Database> & ig,
                                                                             Connector & conn,
                                                                             PropEngine & propEngine)
      : ok(true),
        stat(stat),
        conn(conn)
{

}
template <typename Database, typename Connector, typename PropEngine>
inline NoClauseExchanger<Database, Connector, PropEngine>::~NoClauseExchanger()
{

}

template <typename Database, typename Connector, typename PropEngine>
inline bool NoClauseExchanger<Database, Connector, PropEngine>::setFinished(lbool const res)
{
   assert(!res.isUndef());
   CTSat::lbool r = (res.isTrue()) ? CTSat::lbool::True() : CTSat::lbool::False();
   return conn.setFinished(r);
}
template <typename Database, typename Connector, typename PropEngine>
inline bool NoClauseExchanger<Database, Connector, PropEngine>::isFinished() const
{
   return conn.isFinished();
}

template <typename Database, typename Connector, typename PropEngine>
inline void NoClauseExchanger<Database, Connector, PropEngine>::clauseUsedInConflict(CRef const ref)
{
}

template <typename Database, typename Connector, typename PropEngine>
inline void NoClauseExchanger<Database, Connector, PropEngine>::clauseLearnt(CRef const ref)
{
}

template <typename Database, typename Connector, typename PropEngine>
inline void NoClauseExchanger<Database, Connector, PropEngine>::clauseImproved(CRef const ref)
{
}

template <typename Database, typename Connector, typename PropEngine>
inline void NoClauseExchanger<Database, Connector, PropEngine>::unitLearnt(const Lit l)
{
}
}

#endif /* SOURCES_PARALLEL_NOCLAUSEEXCHANGER_H_ */
