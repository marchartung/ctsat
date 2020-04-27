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

#ifndef SOURCES_PARALLEL_CONFLICTEXCHANGER_H_
#define SOURCES_PARALLEL_CONFLICTEXCHANGER_H_

#include "exchange/SimpleClauseExchanger.h"
#include "core/ImplicationGraph.h"
#include "core/Statistic.h"
#include "initial/SolverConfig.h"
#include "mtl/Vec.h"

namespace ctsat
{

/**
 * Imports clauses in a inbetween buffer, but adds them to the Propagate enigne. When a clause was used in a conflict,
 * it will be available for the solver to import
 */
template <typename Database, typename Connector, typename PropEngine>
class ConflictExchange : public SimpleClauseExchanger<Database, Connector, PropEngine>
{
   typedef SimpleClauseExchanger<Database, Connector, PropEngine> Super;
   typedef typename Super::Clause Clause;
   typedef typename Super::Lit Lit;
   typedef typename Super::CRef CRef;
 protected:
   typedef typename Connector::size_type size_type;
 public:

   static bool attachesClauses()
   {
      return true;
   }

   //static_assert(false);
   ConflictExchange(
                    SolverConfig const & config,
                    Statistic & stat,
                    Database & db,
                    ImplicationGraph<Database> & ig,
                    Connector & conn,
                    PropEngine & propEngine);
   ~ConflictExchange();

   void clauseUsedInConflict(CRef const ref);
   void clauseLearnt(CRef const ref);

   void conflictFound();

   bool hasImportClauses() const;
   std::tuple<bool, CRef> getImportClause();
   void fetchClauses();
   void relocAll(Database& to);

 protected:
   bool allowNonChronoTrail;
   bool onlyExportWhenMin;
   int minAttachLevel;
   int survivedFetchClauses;
   uint64_t numConflictsTillDelete;
   uint64_t lastcleanUp;
   uint64_t lastFetch;
   PropEngine & propEngine;
   vec<CRef> nonConflClauses;
   vec<CRef> conflClauses;

   void exportClause(Clause & c);

   void cleanNonConflictClauses();
};
template <typename Database, typename Connector, typename PropEngine>
void ConflictExchange<Database, Connector, PropEngine>::exportClause(Clause & c)
{
   c.setExport(0);
   Super::exportClause(c);
}

template <typename Database, typename Connector, typename PropEngine>
inline ConflictExchange<Database, Connector, PropEngine>::ConflictExchange(
                                                                           const SolverConfig& config,
                                                                           Statistic& stat,
                                                                           Database& db,
                                                                           ImplicationGraph<Database>& ig,
                                                                           Connector& conn,
                                                                           PropEngine& propEngine)
      : Super(config, stat, db, ig, conn, propEngine),
        allowNonChronoTrail(config.chrono > -1),
        onlyExportWhenMin(config.onlyExportWhenMin),
        minAttachLevel(Super::ig.nVars()),
        survivedFetchClauses(0),
        numConflictsTillDelete(config.numConflictsToDelete),
        lastcleanUp(0),
        lastFetch(0),
        propEngine(propEngine)
{
}

template <typename Database, typename Connector, typename PropEngine>
inline ConflictExchange<Database, Connector, PropEngine>::~ConflictExchange()
{
}

template <typename Database, typename Connector, typename PropEngine>
inline void ConflictExchange<Database, Connector, PropEngine>::clauseUsedInConflict(CRef const ref)
{
   Clause& c = Super::db[ref];
   if (c.getExport() != 0)
   {
      if (c.getLearnt() == 2)
      {
         if (c.getExport() == 2)
         {
            c.set_lbd(Super::ig.computeLBD(c));
            conflClauses.push(ref);
            ++Super::stat.nPromoted;
            --Super::stat.nHoldBackImported;
         }
         c.setExport(1);
      } else
      {
         if (c.getExport() == 2)
            c.setExport(1);
         else if ((!onlyExportWhenMin || c.simplified() ) && c.getExport() == 1
            && c.size() <= Super::max_export_sz
            && c.lbd() <= Super::max_export_lbd)
            exportClause(c);
      }
   }
}

template <typename Database, typename Connector, typename PropEngine>
void ConflictExchange<Database, Connector, PropEngine>::conflictFound()
{
   if (Super::stat.conflicts - lastcleanUp > numConflictsTillDelete / 2)
      cleanNonConflictClauses();
   if (Super::stat.conflicts - lastFetch > 100
      || (survivedFetchClauses > 50 && Super::ig.decisionLevel() <= minAttachLevel))
      fetchClauses();
}

template <typename Database, typename Connector, typename PropEngine>
inline void ConflictExchange<Database, Connector, PropEngine>::clauseLearnt(CRef const ref)
{
   Clause & c = Super::db[ref];
   if (!onlyExportWhenMin && c.lbd() < 3 && c.size() <= Super::max_export_sz)
      exportClause(c);
   else
   {
      c.setExport(2);  // marker for conflict usage
      c.touched() = Super::stat.conflicts;
   }
}

template <typename Database, typename Connector, typename PropEngine>
inline bool ConflictExchange<Database, Connector, PropEngine>::hasImportClauses() const
{
   return conflClauses.size() > 0;
}

template <typename Database, typename Connector, typename PropEngine>
inline std::tuple<bool, typename ConflictExchange<Database, Connector, PropEngine>::CRef> ConflictExchange<
      Database, Connector, PropEngine>::getImportClause()
{
   std::tuple<bool, CRef> res = std::make_tuple(true, Database::npos());
   if (conflClauses.size() > 0)
   {
      std::get<1>(res) = conflClauses.last();
      conflClauses.pop();
   }
   return res;
}

template <typename Database, typename Connector, typename PropEngine>
inline void ConflictExchange<Database, Connector, PropEngine>::fetchClauses()
{
   Super::fetchClauses();
   int i = 0, j = 0, level;
   vec<CRef> & clauses = Super::clauses;
   minAttachLevel = Super::ig.nVars();
   for (; i < clauses.size(); ++i)
   {
      CRef const ref = clauses[i];

      if ((level = propEngine.attachLevel(ref)) > Super::ig.decisionLevel())
      {
         Clause & c = Super::db[ref];
         c.setExport(2);
         nonConflClauses.push(ref);
         ++Super::stat.nHoldBackImported;

         c.setLearnt(2);
         propEngine.safeAttachClause(ref);
      } else
      {
         clauses[j++] = clauses[i];
         minAttachLevel = (level < minAttachLevel) ? level : minAttachLevel;
      }
   }
   clauses.shrink(i - j);
   survivedFetchClauses = j;
   lastFetch = Super::stat.conflicts;
}

template <typename Database, typename Connector, typename PropEngine>
inline void ConflictExchange<Database, Connector, PropEngine>::relocAll(Database& to)
{
   Super::relocAll(to);
   for (int i = 0; i < conflClauses.size(); ++i)
      Super::db.reloc(conflClauses[i], to);
   for (int i = 0; i < nonConflClauses.size(); ++i)
      Super::db.reloc(nonConflClauses[i], to);
}

template <typename Database, typename Connector, typename PropEngine>
inline void ctsat::ConflictExchange<Database, Connector, PropEngine>::cleanNonConflictClauses()
{
   const uint32_t nConfl = Super::stat.conflicts;
   int i = 0, j = 0;
   for (; i < nonConflClauses.size(); ++i)
   {
      Clause & c = Super::db[nonConflClauses[i]];
      if (c.getLearnt() == 2 && c.getExport() == 2)
      {
         uint32_t const age = nConfl - c.touched();
         if (Super::ig.locked(c)
            || age < numConflictsTillDelete
            || (c.lbd() < 4 && age < (5-c.lbd()) * numConflictsTillDelete))
            nonConflClauses[j++] = nonConflClauses[i];
         else
         {
            --Super::stat.nHoldBackImported;
            Super::db.remove(nonConflClauses[i]);
            propEngine.detachClause(nonConflClauses[i]);  // no actual remove, wait for garbage collect
         }
      }
   }

   nonConflClauses.shrink(i - j);
   lastcleanUp = Super::stat.conflicts;
}

}

#endif /* SOURCES_PARALLEL_CONFLICTEXCHANGER_H_ */
