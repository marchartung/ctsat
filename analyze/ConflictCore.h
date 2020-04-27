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

#ifndef ANALYZE_CONFLICTCORE_H_
#define ANALYZE_CONFLICTCORE_H_

#include "minimize/BasicInprocess.h"
#include "analyze/FirstUipAnalyze.h"

namespace ctsat
{
template <typename Propagate>
class ConflictCoreAnalyze;
template <typename Propagate>
class LevelAwareAnalyze;

// this is just a helper, is cannot analyze a conflict with asserting clauses
template <typename Database>
class ConflictCore
{
   typedef typename Database::Lit Lit;
   typedef typename Database::Var Var;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;
 public:
   template <typename Propagate>
   friend class ConflictCoreAnalyze;
   template <typename Propagate>
   friend class LevelAwareAnalyze;

   ConflictCore(SolverConfig const & config, Database & ca, ImplicationGraph<Database> & ig);

   template <typename InprocessMinimize>
   void run(CRef const confl, int const conflictLevel, InprocessMinimize & inConflictMinimize);

   bool hasLearntClause() const;
   LearntClause<Database> const & getLearntClause();

 protected:
   bool hasClause;
   Database & ca;
   ImplicationGraph<Database> & ig;

   LearntClause<Database> lc;
};

template <typename Database>
ConflictCore<Database>::ConflictCore(
                                     SolverConfig const & config,
                                     Database & ca,
                                     ImplicationGraph<Database> & ig)
      : hasClause(false),
        ca(ca),
        ig(ig)
{

}

template <typename Database>
template <typename InprocessMinimize>
inline void ConflictCore<Database>::run(
                                        const CRef confl,
                                        int const conflictLevel,
                                        InprocessMinimize & inConflictMinimize)
{
   int pathC = 0;
   unsigned numResolvents = 0, numBinResolvents = 0, numSkipped = 0;
   vec<Lit> & outC = lc.c;
   outC.clear();
   Lit p = Lit::Undef();
   CRef cref = confl;
   int index = ig.nAssigns();

   do
   {
      assert(cref != Database::npos());  // (otherwise should be UIP)
      Clause & c = ca[cref];

      // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
      if (p != Lit::Undef() && c.size() == 2 && ig.value(c[0]).isFalse())
      {
         assert(ig.value(c[1]).isTrue());
         Lit const tmp = c[0];
         c[0] = c[1], c[1] = tmp;
      }

      bool use = true;
      if (cref != confl)
         for (int i = 0; i < c.size(); ++i)
         {
            Var const v = c[i].var();
            int const level = ig.level(v);
            if (!ig.isSeen(v) && level < conflictLevel && level > 0)
            {
               outC.push(~p);  // collect not resolved lits
               assert(ig.isSeen(p.var()));
               ++numSkipped;
               use = false;
               break;
            }
         }

      if (use)
      {
         ++numResolvents;
         numBinResolvents += c.size() == 2;

         for (int j = (p == Lit::Undef()) ? 0 : 1; j < c.size(); j++)
         {
            Var const v = c[j].var();

            if (!ig.isSeen(v) && ig.level(v) > 0)
            {
               ig.setSeen(v);
               pathC += ig.level(v) >= conflictLevel;
               ig.markSeenToClear(c[j]);
            }
         }
      }
      if (pathC == 0)
         break;
      // Select next clause to look at:
      do
      {
         while (!ig.isSeen(ig.getTrailLit(--index).var()))
            ;
         p = ig.getTrailLit(index);
      } while (ig.level((p.var())) < conflictLevel);

      Var const pv = p.var();
      cref = ig.reason(pv);
      --pathC;
   } while (ig.reason(p) != Database::npos());
   // Generate conflict clause:
   //
   if (numSkipped > 0 && numResolvents > 1 && numResolvents > numBinResolvents)
   {
      if (ig.reason(p) == Database::npos())
         outC.push(~p);
      // add lits from conflict clause:
      Clause const & conflClause = ca[confl];
      for (int i = 0; i < conflClause.size(); ++i)
      {
         int const lvl = ig.level(conflClause[i]);
         if (lvl > 0 && lvl < conflictLevel)
            outC.push(conflClause[i]);
      }

      lc.lbd = inConflictMinimize.run(outC);
      lc.isAsserting = outC.size() == 1 || ig.level(outC[1]) < conflictLevel;  // all but one lit were eliminated during minimize
      // prepare for asserting
      if (lc.isAsserting && outC.size() > 1)
      {
         assert(ig.level(outC[1]) > 0);
         int max_i = 1;
         for (int i = 2; i < outC.size(); i++)
         {
            if (ig.level(outC[i]) > ig.level(outC[max_i]))
               max_i = i;
            assert(ig.level(outC[i]) > 0);
         }
         std::swap(outC[1], outC[max_i]);
      }
      hasClause = true;
   } else
   {
      hasClause = false;
      outC.clear();
   }
}

template <typename Database>
inline bool ConflictCore<Database>::hasLearntClause() const
{
   return hasClause;
}

template <typename Database>
inline const LearntClause<Database>& ConflictCore<Database>::getLearntClause()
{
   assert(hasClause);
   hasClause = false;
   return lc;
}
}

#endif /* ANALYZE_CONFLICTCORE_H_ */
