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

#ifndef ANALYZE_LEVELAWAREANALYZE_H_
#define ANALYZE_LEVELAWAREANALYZE_H_

#include "analyze/MultiUipAnalyze.h"
#include "analyze/ConflictCore.h"
#include "mtl/AvgQueue.h"

namespace ctsat
{

template <typename Propagate>
class LevelAwareAnalyze : public MultiUipAnalyze<Propagate>
{
   typedef MultiUipAnalyze<Propagate> MUA;
   typedef FirstUipAnalyze<Propagate> FUA;
   typedef typename Propagate::Database Database;
   typedef typename Database::Lit Lit;
   typedef typename Database::Var Var;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;

 public:
   LevelAwareAnalyze(
                     SolverConfig const & config,
                     Database & ca,
                     ImplicationGraph<Database> & ig,
                     Propagate & propEngine);

   template <typename Callbacks>
   void run(CRef const confl, Callbacks const & in);

   bool hasLearntClause() const;
   LearntClause<Database> const & getLearntClause();

   constexpr bool learntsMultipleClauses() const;

   void notifyRestart();

 protected:
   bool const alwaysSwap;
   int const levelDiffEnforce;
   int numInitialConflicts;
   int nConflictsSinceRestart;
   ConflictCore<Database> cc;
   AvgQueue<int> levelQueue;

};

template <typename Propagate>
inline LevelAwareAnalyze<Propagate>::LevelAwareAnalyze(
                                                       const SolverConfig& config,
                                                       Database& ca,
                                                       ImplicationGraph<Database>& ig,
                                                       Propagate& propEngine)
      : MUA(config, ca, ig, propEngine),
        alwaysSwap(config.LAA_alwaySwap),
        levelDiffEnforce(config.LAA_levelDiffEnforce),
        numInitialConflicts(config.LAA_numInitialConflicts),
        nConflictsSinceRestart(0),
        cc(config, ca, ig),
        levelQueue(config.LAA_levelQueueSz)
{
}

template <typename Propagate>
void LevelAwareAnalyze<Propagate>::notifyRestart()
{
   levelQueue.clear();
}

template <typename Propagate>
template <typename Callbacks>
inline void LevelAwareAnalyze<Propagate>::run(const CRef confl, const Callbacks& in)
{
   int curIndex = FUA::analyzeFistUip(confl, in);
   FUA::ig.clearSeen();
   int const conflictLevel = FUA::ig.level(FUA::lc.c[0]);
   bool const addClauses = (levelDiffEnforce > 0)
      && (--numInitialConflicts >= 0
         || !levelQueue.full()
         || levelQueue.avg() - conflictLevel >= levelDiffEnforce);
   if (FUA::lc.c.size() > 1 && (addClauses || alwaysSwap))
   {
      cc.run(confl, conflictLevel, FUA::inConflictMinimize);
      FUA::ig.clearSeen();
      if (cc.hasClause)
      {
         if (cc.lc.isAsserting)
         {
            if (cc.lc.c.size() < FUA::lc.c.size()) // the cc clause is better than the fuip
               cc.lc.swap(FUA::lc);
            else if(cc.lc.c[0] == FUA::lc.c[0])
               cc.hasClause = false; // the chance they are equal is high, so remove the cc clause
         } else if (cc.lc.c.size() > FUA::lc.c.size())
            cc.hasClause = false; // longer clauses and not asserting ... for sure meaningless
         else
            cc.lc.lbd = cc.lc.c.size()-1; // since the clause is not asserting, the lbd value some kine of false
      }

      MUA::collectMuipLearnts(curIndex, conflictLevel);
      for (int i = 0; i < MUA::lcPos; ++i)
         if (MUA::lcs[i].isAsserting && MUA::lcs[i].c.size() < FUA::lc.c.size())
            FUA::lc.swap(MUA::lcs[i]);
      if (!addClauses) // check if we only wanted to find a better clause, if so we delete the additional clauses
      {
         MUA::lcReadPos = MUA::lcPos;
         cc.hasClause = false;
      }

   }
   levelQueue.push(conflictLevel);

}

template <typename Propagate>
inline bool LevelAwareAnalyze<Propagate>::hasLearntClause() const
{
   return MUA::hasLearntClause() || cc.hasLearntClause();
}

template <typename Propagate>
inline const LearntClause<typename Propagate::Database> & LevelAwareAnalyze<Propagate>::getLearntClause()
{
   assert(hasLearntClause());
   return (MUA::hasLearntClause()) ? MUA::getLearntClause() : cc.getLearntClause();
}

template <typename Propagate>
inline constexpr bool LevelAwareAnalyze<Propagate>::learntsMultipleClauses() const
{
   return true;
}
}

#endif /* ANALYZE_LEVELAWAREANALYZE_H_ */
