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

#ifndef ANALYZE_FIRSTUIPANALYZE_H_
#define ANALYZE_FIRSTUIPANALYZE_H_

#include "mtl/Vec.h"
#include "core/ImplicationGraph.h"
#include "initial/SolverConfig.h"

#include "minimize/BasicInprocess.h"

namespace ctsat
{
template<typename Database>
struct LearntClause
{
   typedef typename Database::Lit Lit;
   LearntClause(LearntClause const&) = delete;
   LearntClause(LearntClause &&) = delete;
   LearntClause& operator=(LearntClause const&) = delete;
   LearntClause& operator=(LearntClause &&) = delete;

   bool isAsserting;
   int lbd;
   vec<Lit> c;

   LearntClause()
         : isAsserting(false),
           lbd(-1)
   {
   }

   void swap(LearntClause & lc)
   {
      std::swap(isAsserting, lc.isAsserting);
      std::swap(lbd, lc.lbd);
      c.swap(lc.c);
   }
};

template <typename Propagate>
class FirstUipAnalyze
{
   typedef typename Propagate::Database Database;
   typedef typename Database::Lit Lit;
   typedef typename Database::Var Var;
   typedef typename Database::Clause Clause;
   typedef typename Database::CRef CRef;
   typedef typename Database::lbool lbool;

 public:
   FirstUipAnalyze(
                   SolverConfig const & config,
                   Database & ca,
                   ImplicationGraph<Database> & ig,
                   Propagate & propEngine);

   template <typename Callbacks>
   void run(CRef const confl, Callbacks const & in);

   bool hasLearntClause() const;
   LearntClause<Database> const & getLearntClause();

   constexpr bool learntsMultipleClauses() const;

   void notifyRestart()
   {

   }

 protected:
   bool hasClause;
   Database & ca;
   ImplicationGraph<Database> & ig;
   BasicInprocessMinimize<Propagate> inConflictMinimize;

   LearntClause<Database> lc;

   template <typename Callbacks>
   int analyzeFistUip(CRef confl, Callbacks const & cb);
};

template <typename Propagate>
inline FirstUipAnalyze<Propagate>::FirstUipAnalyze(
                                                   SolverConfig const & config,
                                                   Database & ca,
                                                   ImplicationGraph<Database> & ig,
                                                   Propagate & propEngine)
      : hasClause(false),
        ca(ca),
        ig(ig),
        inConflictMinimize(config, ig, propEngine)
{
}

template <typename Propagate>
template <typename Callbacks>
inline void FirstUipAnalyze<Propagate>::run(const CRef confl, Callbacks const & in)
{
   analyzeFistUip(confl, in);
   ig.clearSeen();
}

template <typename Propagate>
template <typename Callbacks>
inline int FirstUipAnalyze<Propagate>::analyzeFistUip(CRef confl, Callbacks const & cb)
{
   lc.c.clear();

   int pathC = 0;
   int const conflictLvl = ig.level(ca[confl][0]);
   vec<Lit> & outC = lc.c;
   auto const onClauseVisit = [&](Lit const p, CRef const ref)
   {  cb.onClauseVisit(ref);
      return true;};
   auto const onLitVisit = [&](Lit const l)
   {
      Var const v = l.var();
      cb.onVarVisit(v);
      if (ig.level(v) >= conflictLvl)
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
   outC.push();      // (leave room for the asserting literal)
   int const index = ig.visitImplications(ig.nAssigns() - 1, confl, conflictLvl, onClauseVisit,
                                          onLitVisit, runFurther);
   outC[0] = ~ig.getTrailLit(index);

   lc.lbd = inConflictMinimize.run(outC) - 1;

// prepare for asserting
   if (outC.size() > 1)
   {
      int max_i = 1;
      assert(ig.level(outC[1]) > 0);
      for (int i = 2; i < outC.size(); i++)
      {
         if (ig.level(outC[i]) > ig.level(outC[max_i]))
            max_i = i;
         assert(ig.level(outC[i]) > 0);
      }
      std::swap(outC[1], outC[max_i]);
   }
   cb.onClauseCreated(outC);

   hasClause = true;
   lc.isAsserting = true;
   return index;
}

template <typename Propagate>
inline bool FirstUipAnalyze<Propagate>::hasLearntClause() const
{
   return hasClause;
}

template <typename Propagate>
inline const LearntClause<typename Propagate::Database>& FirstUipAnalyze<Propagate>::getLearntClause()
{
   assert(hasClause);
   hasClause = false;
   return lc;
}

template <typename Propagate>
inline constexpr bool FirstUipAnalyze<Propagate>::learntsMultipleClauses() const
{
   return false;
}
}

#endif /* ANALYZE_FIRSTUIPANALYZE_H_ */
