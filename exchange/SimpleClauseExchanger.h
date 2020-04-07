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

#ifndef SOURCES_PARALLEL_SIMPLECLAUSEEXCHANGER_H_
#define SOURCES_PARALLEL_SIMPLECLAUSEEXCHANGER_H_

#include "exchange/NoClauseExchanger.h"
#include "core/ImplicationGraph.h"
#include "core/Statistic.h"
#include "initial/SolverConfig.h"
#include "initial/Inputs.h"
#include "mtl/Vec.h"

namespace CTSat
{

template <typename Database>
struct ExportClause
{
   typedef typename Database::Clause Clause;
   typedef typename Database::Lit Lit;

   ExportClause() = delete;
   ExportClause(ExportClause<Database> &&) = delete;
   ExportClause<Database>& operator=(ExportClause<Database> const&) = delete;
   ExportClause<Database>& operator=(ExportClause<Database> &&) = delete;

   uint16_t id;
   uint16_t lbd;
   int const sz;
   Lit clause[1];
   ExportClause(Clause const & c, unsigned const lbd, unsigned const id);
   ExportClause(ExportClause<Database> const & in);
   ExportClause(Lit const & l, unsigned const id);
   static uint64_t nbytes(Clause const & c);
   static uint64_t nbytes(Lit const & l);
   static uint64_t nbytes(ExportClause<Database> const & c);
   int size() const
   {
      return sz;
   }
   Lit const & operator[](int const idx) const
   {
      return clause[idx];
   }

   void set(Clause & c, bool const minimizeClause) const
   {
      c.set_lbd(lbd);
      c.setSimplified(!minimizeClause);
   }


};

template <typename Database, typename Connector, typename PropEngine>
class SimpleClauseExchanger : public NoClauseExchanger<Database, Connector, PropEngine>
{
 protected:
   typedef typename Database::Clause Clause;
   typedef typename Database::Lit Lit;
   typedef typename Database::CRef CRef;
   typedef ExportClause<Database> ExClause;
   typedef NoClauseExchanger<Database, Connector, PropEngine> Super;
   typedef typename Connector::size_type size_type;
 public:

   SimpleClauseExchanger(
                         SolverConfig const & config,
                         Statistic & stat,
                         Database & db,
                         ImplicationGraph<Database> & ig,
                         Connector & conn,
                         PropEngine & propEngine);
   ~SimpleClauseExchanger();

   void clauseLearnt(CRef const ref);
   void unitLearnt(Lit const l);

   bool shouldFetch() const;
   bool hasImportClauses() const;
   void fetchClauses();
   std::tuple<bool,CRef> getImportClause();
   Lit getImportUnit();

   void relocAll(Database& to);

 protected:
   bool const minimize_import_cl;
   int max_export_lbd;
   int max_export_sz;
   const unsigned id;
   size_type curReadPos;
   Database & db;
   ImplicationGraph<Database> & ig;
   vec<Lit> units;
   vec<CRef> clauses;
   vec<Lit> tmpClause;

   bool prepClause(vec<Lit> & preped, ExClause const & c)
   {
      // 1. copy clause to preped
      // 2. sort it so the highest level are in front
      // 3 remove set lits
      // 4. return true, when clause is sat
      preped.clear();
      for (int i = 0; i < c.size(); ++i)
      {
         lbool const val = ig.value(c[i]);
         if (val.isUndef() || ig.level(c[i].var()) > 0)
            preped.push(c[i]);
         else if (val.isTrue())
            return true;
      }
      return false;
   }

   template <typename ... Args>
   void exportClause(uint64_t const nBytes, Args ... args)
   {
      Super::conn.template exchange<ExClause, Args...>(nBytes, args...);
      ++Super::stat.nSendClauses;
      if (Super::conn.shouldImport(curReadPos))
         fetchClauses();
   }

   void exportClause(Clause const & c);

};
template <typename Database, typename Connector, typename PropEngine>
inline bool SimpleClauseExchanger<Database, Connector, PropEngine>::shouldFetch() const
{
   return Super::conn.shouldImport(curReadPos);
}

template <typename Database, typename Connector, typename PropEngine>
void SimpleClauseExchanger<Database, Connector, PropEngine>::exportClause(Clause const & c)
{
   exportClause<Clause const &, unsigned const, unsigned const>(ExportClause<Database>::nbytes(c),
                                                                c, c.lbd(), id);
}

template <typename Database, typename Connector, typename PropEngine>
void SimpleClauseExchanger<Database, Connector, PropEngine>::fetchClauses()
{
   while (Super::conn.isValid(curReadPos))
   {
      ExClause const & importCl = Super::conn.template get<ExClause>(curReadPos);
      if (importCl.id != id && !prepClause(tmpClause, importCl))
      {
         if (tmpClause.size() < 2)
         {
            if (tmpClause.size() == 0)
            {
               Super::ok = false;
               return;
            }
            units.push(tmpClause[0]);
         } else
         {
            CRef const ref = db.alloc(tmpClause, true);
            importCl.set(db[ref], minimize_import_cl);
            clauses.push(ref);
         }
         ++Super::stat.nReceivedClauses;
      }
      curReadPos = Super::conn.next(curReadPos);
   }
}

template <typename Database, typename Connector, typename PropEngine>
void SimpleClauseExchanger<Database, Connector, PropEngine>::relocAll(Database& to)
{
   for (int i = 0; i < clauses.size(); ++i)
      db.reloc(clauses[i], to);
}

template <typename Database, typename Connector, typename PropEngine>
inline void SimpleClauseExchanger<Database, Connector, PropEngine>::clauseLearnt(CRef const ref)
{
   Clause & c = db[ref];
   if (c.lbd() <= max_export_lbd && c.size() < max_export_sz)
      exportClause(c);
}

template <typename Database, typename Connector, typename PropEngine>
inline void SimpleClauseExchanger<Database, Connector, PropEngine>::unitLearnt(const Lit l)
{
   exportClause<Lit const &, unsigned const>(ExClause::nbytes(l), l, id);
}

template <typename Database>
ExportClause<Database>::ExportClause(Clause const & c, unsigned const lbd, unsigned const id)
      : id(id),
        lbd(lbd),
        sz(c.size())
{
   for (int i = 0; i < sz; ++i)
      clause[i] = c[i];
}

template <typename Database>
ExportClause<Database>::ExportClause(ExportClause<Database> const & c)
      : id(c.id),
        lbd(c.lbd),
        sz(c.size())
{
   for (int i = 0; i < sz; ++i)
      clause[i] = c[i];
}
template <typename Database>
ExportClause<Database>::ExportClause(Lit const & l, unsigned const id)
      : id(id),
        lbd(0),
        sz(1)
{
   clause[0] = l;
}

template <typename Database, typename Connector, typename PropEngine>
inline SimpleClauseExchanger<Database, Connector, PropEngine>::SimpleClauseExchanger(
                                                                                     SolverConfig const & config,
                                                                                     Statistic & stat,
                                                                                     Database& db,
                                                                                     ImplicationGraph<
                                                                                           Database>& ig,
                                                                                     Connector & conn,
                                                                                     PropEngine & propEngine)
      : Super(config, stat, db, ig, conn, propEngine),
        minimize_import_cl(config.minimize_import_cl),
        max_export_lbd(config.max_export_lbd),
        max_export_sz(config.max_export_sz),
        id(conn.getUniqueId()),
        curReadPos(0),
        db(db),
        ig(ig)
{
}

template <typename Database, typename Connector, typename PropEngine>
inline SimpleClauseExchanger<Database, Connector, PropEngine>::~SimpleClauseExchanger()
{
}

template <typename Database>
inline uint64_t ExportClause<Database>::nbytes(const Clause& c)
{
   return sizeof(ExportClause<Database> ) + (c.size() - 1) * sizeof(Lit);
}

template <typename Database>
inline uint64_t ExportClause<Database>::nbytes(const ExportClause<Database>& c)
{
   return sizeof(ExportClause<Database> ) + (c.size() - 1) * sizeof(Lit);
}

template <typename Database>
inline uint64_t ExportClause<Database>::nbytes(const Lit& l)
{
   return sizeof(ExportClause<Database> );
}

template <typename Database, typename Connector, typename PropEngine>
inline std::tuple<bool,typename SimpleClauseExchanger<Database, Connector, PropEngine>::CRef> SimpleClauseExchanger<
      Database, Connector, PropEngine>::getImportClause()
{
   CRef res = Database::npos();
   if (clauses.size() > 0)
   {
      res = clauses.last();
      clauses.pop();
   }
   return std::make_tuple(false,res);
}

template <typename Database, typename Connector, typename PropEngine>
inline bool SimpleClauseExchanger<Database, Connector, PropEngine>::hasImportClauses() const
{
   return clauses.size() > 0 || units.size() > 0;
}

template <typename Database, typename Connector, typename PropEngine>
inline typename SimpleClauseExchanger<Database, Connector, PropEngine>::Lit SimpleClauseExchanger<
      Database, Connector, PropEngine>::getImportUnit()
{
   Lit res = Lit::Undef();
   if (units.size() > 0)
   {
      res = units.last();
      units.pop();
   }
   return res;
}
}

#endif /* SOURCES_PARALLEL_SIMPLECLAUSEEXCHANGER_H_ */
