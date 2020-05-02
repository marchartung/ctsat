/*
 * BasicInprocess.h
 *
 *  Created on: 11.04.2020
 *      Author: hartung
 */

#ifndef MINIMIZE_BASICINPROCESS_H_
#define MINIMIZE_BASICINPROCESS_H_

#include "initial/SolverConfig.h"
#include "core/ImplicationGraph.h"

namespace ctsat
{
template <typename Propagate>
class BasicInprocessMinimize
{

   typedef typename Propagate::Database Database;
   typedef typename Database::Lit Lit;
   typedef typename Database::Var Var;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
 public:

   BasicInprocessMinimize(
                          SolverConfig const & config,
                          ImplicationGraph<Database> & ig,
                          Propagate & propEngine);

   // returns lbd of minimized clause
   template <typename LitVec>
   int run(LitVec & c, int const startIdx = 1, bool const containsTrue = false);

 private:
   ImplicationGraph<Database> & ig;
   Propagate & propEngine;

   int const ccmin_mode;   // Controls conflict clause minimization (0=none, 1=basic, 2=deep).
   int const maxEntendedBinaryResolutionSz;
   int const maxFullImplicationMinLbd;

};

template <typename Propagate>
inline BasicInprocessMinimize<Propagate>::BasicInprocessMinimize(
                                                                 const SolverConfig& config,
                                                                 ImplicationGraph<Database>& ig,
                                                                 Propagate& propEngine)
      : ig(ig),
        propEngine(propEngine),
        ccmin_mode(config.ccmin_mode),
        maxEntendedBinaryResolutionSz(config.maxEntendedBinaryResolutionSz),
        maxFullImplicationMinLbd(config.maxFullImplicationMinLbd)
{
}

template <typename Propagate>
template <typename LitVec>
inline int BasicInprocessMinimize<Propagate>::run(
                                                  LitVec & c,
                                                  int const startIdx,
                                                  bool const containsTrue)
{
   int lbd = ig.computeLBD(c);
   // minimize conflict clause:
   if (ccmin_mode == 2)
      ig.minimizeImpliedRecursive(c, lbd <= maxFullImplicationMinLbd, startIdx);
   else if (ccmin_mode == 1)
      ig.minimizeImplied(c);

   // TODO write version for true literals:
   if (!containsTrue)
   {
      if (c.size() <= maxEntendedBinaryResolutionSz)
      {
         propEngine.extendedBinResMinimize(c, startIdx);
         lbd = std::min(ig.computeLBD(c), lbd + 1);
      } else if (lbd <= 6 && c.size() <= 30)  // Try further minimization?
         if (propEngine.binResMinimize(c, startIdx))
            lbd = std::min(ig.computeLBD(c), lbd + 1);
   }
   return lbd;
}

}

#endif /* MINIMIZE_BASICINPROCESS_H_ */
