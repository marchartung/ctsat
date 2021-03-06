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

#ifndef ANALYZE_MULTIUIPANALYZE_H_
#define ANALYZE_MULTIUIPANALYZE_H_

#include "mtl/Vec.h"
#include "core/ImplicationGraph.h"
#include "initial/SolverConfig.h"

#include "minimize/BasicInprocess.h"
#include "analyze/FirstUipAnalyze.h"

namespace ctsat
{
template <typename Propagate>
class MultiUipAnalyze : public FirstUipAnalyze<Propagate>
{
   typedef FirstUipAnalyze<Propagate> FUA;
   typedef typename Propagate::Database Database;
   typedef typename Database::Lit Lit;
   typedef typename Database::Var Var;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;

 public:
   MultiUipAnalyze(
                   SolverConfig const & config,
                   Database & ca,
                   ImplicationGraph<Database> & ig,
                   Propagate & propEngine);

   template <typename Callbacks>
   void run(CRef const confl, Callbacks const & in);

   bool hasLearntClause() const;
   LearntClause<Database> const & getLearntClause();

   constexpr bool learntsMultipleClauses() const;

 protected:
   int lcReadPos;
   int lcPos;

   vec<LearntClause<Database>> lcs;

   void collectMuipLearnts(int index, int const conflictLevel);

   int analyzeMuip(int index, int const conflictLevel, LearntClause<Database> & ls);
};

template <typename Propagate>
inline MultiUipAnalyze<Propagate>::MultiUipAnalyze(
                                                   const SolverConfig& config,
                                                   Database& ca,
                                                   ImplicationGraph<Database>& ig,
                                                   Propagate& propEngine)
      : FirstUipAnalyze<Propagate>(config, ca, ig, propEngine),
        lcReadPos(0),
        lcPos(0)
{
}

template <typename Propagate>
template <typename Callbacks>
inline void MultiUipAnalyze<Propagate>::run(const CRef confl, const Callbacks& in)
{
   int const curIndex = FUA::analyzeFistUip(confl, in);
   int const conflictLevel = FUA::ig.level(FUA::lc.c[0]);
   FUA::ig.clearSeen();
   collectMuipLearnts(curIndex,conflictLevel);
}

template <typename Propagate>
inline void MultiUipAnalyze<Propagate>::collectMuipLearnts(int curIndex, const int conflictLevel)
{
   lcPos = 0;
   while (FUA::ig.reason(FUA::ig.getTrailLit(curIndex)) != Database::npos())
   {
      if (lcPos >= lcs.size())
         lcs.growTo(lcs.size() + 1);
      curIndex = analyzeMuip(curIndex, conflictLevel, lcs[lcPos]);
      if (lcs[lcPos].c.size() > 0)
      {
         LearntClause<Database> & mc = lcs[lcPos];
         if(mc.isAsserting && mc.c.size() < FUA::lc.c.size())
            mc.swap(FUA::lc);
         ++lcPos;
      }
   }
   lcReadPos = 0;
}

template <typename Propagate>
inline bool MultiUipAnalyze<Propagate>::hasLearntClause() const
{
   return FirstUipAnalyze<Propagate>::hasLearntClause() || lcReadPos < lcPos;
}

template <typename Propagate>
inline const LearntClause<typename Propagate::Database>& MultiUipAnalyze<Propagate>::getLearntClause()
{
   assert(hasLearntClause());
   return
         (FirstUipAnalyze<Propagate>::hasLearntClause()) ?
               FirstUipAnalyze<Propagate>::getLearntClause() : lcs[lcReadPos++];
}

template <typename Propagate>
inline constexpr bool MultiUipAnalyze<Propagate>::learntsMultipleClauses() const
{
   return true;
}


template <typename Propagate>
inline int MultiUipAnalyze<Propagate>::analyzeMuip(
                                                   int index,
                                                   int const conflictLevel,
                                                   LearntClause<Database>& lc)
{
   ImplicationGraph<Database> & ig = FUA::ig;
   assert((ig.reason(ig.getTrailLit(index)) != Database::npos()));
   int pathC = 0;
   unsigned numResolvents = 0, numBinResolvents = 0;
   Lit const prevUip = ig.getTrailLit(index);
   vec<Lit> & outC = lc.c;
   auto const onClauseVisit = [&](Lit const p, CRef const ref)
   {  ++numResolvents; numBinResolvents += FUA::ca[ref].size() == 2; return true;};
   auto const onLitVisit = [&](Lit const l)
   {
      Var const v = l.var();
      if (ig.level(v) >= conflictLevel)
      {
         ++pathC;
         return false;
      } else
      {
         outC.push(l);
         return true;
      }
   };
   auto const runFurther = [&](Lit const p)
   {
      return --pathC > 0;
   };

   // Generate conflict clause:
   //
   ig.setSeen(prevUip.var());  // mark the prev uip so it wont be in the path
   outC.clear();
   outC.growTo(2);      // leave room for the asserting literal and previous uip
   const int resIdx = ig.visitImplications(index - 1, ig.reason(prevUip), conflictLevel,
                                                onClauseVisit, onLitVisit, runFurther);
   if (numResolvents > 1 && numResolvents > numBinResolvents)
   {
      outC[0] = prevUip;
      outC[1] = ~ig.getTrailLit(resIdx);

      lc.lbd = FUA::inConflictMinimize.run(outC, 2, false);  // 2 ... Uips will be skipped in implication check
      lc.isAsserting = outC.size() == 1 || outC[1].var() != ig.getTrailLit(resIdx).var();
      // prepare for asserting
      if (outC.size() > 1)
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
   } else
   {
      outC.clear();
   }
   ig.unsetSeen(prevUip.var());
   ig.clearSeen();
   return resIdx;
}
}

#endif /* ANALYZE_MULTIUIPANALYZE_H_ */
