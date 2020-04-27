/*
 * Vivification.h
 *
 *  Created on: 11.04.2020
 *      Author: hartung
 */

#ifndef MINIMIZE_VIVIFICATION_H_
#define MINIMIZE_VIVIFICATION_H_

#include "initial/SolverConfig.h"
#include "core/ImplicationGraph.h"

#include <tuple>

namespace ctsat
{

template <typename Propagate>
class Vivification
{

   typedef typename Propagate::Database Database;
   typedef typename Database::Lit Lit;
   typedef typename Database::Var Var;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
 public:
   struct PropagateResult
   {
      bool isTruePropagation;
      CRef confl;
      int numProps;
   };

   Vivification(
                SolverConfig const & config,
                Database & ca,
                ImplicationGraph<Database> & ig,
                Propagate & propEngine);

   bool run(Clause& c);

   template <typename ClauseType>
   PropagateResult propagateClause(ClauseType const & c);

 private:
   Database & ca;
   ImplicationGraph<Database> & ig;
   Propagate & propEngine;

   vec<Lit> simp_learnt_clause;

   // returns lbd
   int simpleAnalyze(CRef confl, vec<Lit>& c, bool const True_confl);
};

template <typename Propagate>
inline Vivification<Propagate>::Vivification(
                                             const SolverConfig& config,
                                             Database & ca,
                                             ImplicationGraph<Database>& ig,
                                             Propagate& propEngine)
      : ca(ca),
        ig(ig),
        propEngine(propEngine)
{
}

template <typename Propagate>
template <typename ClauseType>
inline typename Vivification<Propagate>::PropagateResult Vivification<Propagate>::propagateClause(
                                                                                           ClauseType const & c)
{
   PropagateResult res;
   res.isTruePropagation = false;
   res.confl = Database::npos();
   res.numProps = 0;

   for (int i = 0; i < c.size(); ++i)
   {
      if (ig.value(c[i]).isUndef())
      {
         ++res.numProps;
         propEngine.simpleUncheckEnqueue(~c[i]);
         if ((res.confl = propEngine.simplePropagate()) != Database::npos())
            break;
      } else
      {
         if (ig.value(c[i]).isTrue())
         {
            res.isTruePropagation = true;
            res.confl = ig.reason(c[i]);
            Clause & conflC = ca[res.confl];
            if(!ig.value(conflC[0]).isTrue())
                  std::swap(conflC[0],conflC[1]);
            assert(ig.value(conflC[0]).isTrue());
            break;
         }
      }
   }
   return res;
}

template <typename Propagate>
inline bool Vivification<Propagate>::run(Clause& c)
{

   const int initSize = c.size();
   propEngine.trailRecord = ig.nAssigns();      // record the start pointer

   PropagateResult pr = propagateClause(c);

   if (pr.confl != Database::npos())
   {
      simp_learnt_clause.clear();
      if (pr.isTruePropagation)
         simp_learnt_clause.push(ca[pr.confl][0]);
      simpleAnalyze(pr.confl, simp_learnt_clause, pr.isTruePropagation);

      if (simp_learnt_clause.size() < c.size())
      {
         int i = 0;
         for (; i < simp_learnt_clause.size(); i++)
            c[i] = simp_learnt_clause[i];
         c.shrink(c.size() - i);
      }
   }
   else if(pr.numProps < c.size())
   {
      int i = 0, j = 0;
      for(;i<c.size();++i)
      {
         assert(ig.value(c[i]).isFalse());
         c[j] = c[i];
         j += ig.reason(c[i]) == Database::npos();
      }
      c.shrink(i-j);
   }
   propEngine.cancelUntilTrailRecord();

   if (pr.numProps < c.lbd())
      c.set_lbd(pr.numProps);
   return c.size() < initSize;
}

template <typename Propagate>
inline int Vivification<Propagate>::simpleAnalyze(CRef confl, vec<Lit>& c, const bool True_confl)
{
   int pathC = 0;
   Lit p = Lit::Undef();
   int index = ig.nAssigns() - 1;

   do
   {
      if (confl != Database::npos())
      {
         Clause& c = ca[confl];
         if (p != Lit::Undef() && c.size() == 2 && ig.value(c[0]).isFalse())
         {
            assert(ig.value(c[1]).isTrue());
            Lit tmp = c[0];
            c[0] = c[1], c[1] = tmp;
         }
         // if True_confl==true, then choose p begin with the 1th index of c;
         for (int j = (p == Lit::Undef() && True_confl == false) ? 0 : 1; j < c.size(); j++)
         {
            Var const v = c[j].var();
            if (!ig.isSeen(v))
            {
               ig.setSeen(v);
               pathC++;
            }
         }
      } else if (confl == Database::npos())
      {
         c.push(~p);
      }
      // if not break, while() will come to the index of trail blow 0, and fatal error occur;
      if (pathC == 0)
         break;
      // Select next clause to look at:
      while (!ig.isSeen(ig.getTrailLit(index--).var()))
         ;
      // if the reason cr from the 0-level assigned var, we must break avoid move forth further;
      // but attention that maybe seen[x]=1 and never be clear. However makes no matter;
      if (propEngine.trailRecord > index + 1)
         break;
      p = ig.getTrailLit(index + 1);
      Var const v = p.var();
      confl = ig.reason(v);
      ig.unsetSeen(v);
      pathC--;

   } while (pathC >= 0);
   return c.size() - 1;
}

}

#endif /* MINIMIZE_VIVIFICATION_H_ */
