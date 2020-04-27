/*
 * ClauseFilter.h
 *
 *  Created on: 14.04.2020
 *      Author: hartung
 */

#ifndef MPI_CLAUSEFILTER_H_
#define MPI_CLAUSEFILTER_H_

#include "minimize/Vivification.h"
#include "initial/SolverConfig.h"
#include "initial/SatInstance.h"
#include "core/ImplicationGraph.h"
#include "core/Statistic.h"
#include "mtl/Vec.h"
#include "mpi/AtomicStealResultQueue.h"
#include "exchange/ExportClause.h"

#include <atomic>

namespace ctsat
{


template<typename Database>
struct TestedClause
{
   typedef typename Database::Lit Lit;
   uint16_t initLbd;
   uint16_t foundLbd;
   double propPropability;
   Lit data[1];
};

template<typename Database>
struct StealDistributor
{
   std::atomic<bool> abort;
   std::atomic<unsigned> maxNumViviTries;
   AtomicStealResultQueue<ExportClause<Database> const *, double> stealQueue;

   std::atomic<uint32_t> importRound;
   std::atomic<unsigned> numSolverImported;
   ExportClause<Database> const * importCl;

};

template <typename Database, typename Propagate>
class PropagationFilter
{
   typedef typename Database::Lit Lit;
   typedef typename Database::Var Var;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
   typedef StealDistributor<Database> SD;
   typedef typename StealDistributor<Database>::size_type size_type;

 public:

   PropagationFilter(
                     SolverConfig const & config,
                     SD & distributor,
                     vec<decltype(SatInstance::ca)::lbool> const & model,
                     decltype(SatInstance::ca) const& db,
                     vec<decltype(SatInstance::ca)::CRef> const & clauses);


   bool setOk(bool const b);
   bool isOk() const;

   void run();

 private:
   uint64_t importRound;
   uint64_t numRounds;
   uint64_t lastReduce;
   Statistic stat;
   Database ca;
   ImplicationGraph<Database> ig;
   Vivification<Propagate> vivi;
   SD & distributor;

   vec<CRef> clauses;
   vec<CRef> imports;

   void reduce();
   void import();
   void testClauses();

};
template <typename Database, typename Propagate>
inline void ctsat::PropagationFilter<Database, Propagate>::import()
{
}

template <typename Database, typename Propagate>
inline ctsat::PropagationFilter<Database, Propagate>::PropagationFilter(
                                                                        const SolverConfig& config,
                                                                        StealDistributor<Database>& distributor,
                                                                        const vec<
                                                                              decltype(SatInstance::ca)::lbool>& model,
                                                                        const decltype(SatInstance::ca)& db,
                                                                        const vec<
                                                                              decltype(SatInstance::ca)::CRef>& clauses)
{
}

template <typename Database, typename Propagate>
inline bool ctsat::PropagationFilter<Database, Propagate>::setOk(const bool b)
{
}

template <typename Database, typename Propagate>
inline bool ctsat::PropagationFilter<Database, Propagate>::isOk() const
{
}

template <typename Database, typename Propagate>
inline void ctsat::PropagationFilter<Database, Propagate>::run()
{
   size_type pos = SD::npos();
   while(!distributor.abort)
   {
      testClauses();
      import();
      reduce();

   }
}

template <typename Database, typename Propagate>
inline void ctsat::PropagationFilter<Database, Propagate>::reduce()
{
}
}

#endif /* MPI_CLAUSEFILTER_H_ */
