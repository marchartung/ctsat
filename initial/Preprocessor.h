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

#ifndef Minisat_SimpSolver_h
#define Minisat_SimpSolver_h

#include "database/BasicTypes.h"
#include "mtl/Queue.h"
#include "mtl/Heap.h"
#include "mtl/OccLists.h"
#include "core/ImplicationGraph.h"
#include "core/Statistic.h"
#include "branch/Branch.h"
#include "propagate/MiniSatPropagate.h"
#include "initial/SatInstance.h"
#include "initial/SolverConfig.h"
#include "utils/DratPrint.h"
#include "utils/Random.h"

#include "utils/Exceptions.h"
#include "utils/ParseUtils.h"
#include <zlib.h>

#include "EliminatedClauseDatabase.h"

namespace ctsat
{

//=================================================================================================

class Preprocessor
{
   typedef ClauseAllocator::Lit Lit;
   typedef ClauseAllocator::Var Var;
   typedef ClauseAllocator::Clause Clause;
   typedef ClauseAllocator::CRef CRef;
   typedef ClauseAllocator::lbool lbool;

   Preprocessor() = delete;
   Preprocessor(Preprocessor&&) = delete;
   Preprocessor(Preprocessor const&) = delete;
   Preprocessor& operator=(Preprocessor&&) = delete;
   Preprocessor& operator=(Preprocessor const &) = delete;
 public:
   // Constructor/Destructor:
   //
   Preprocessor(SolverConfig const & config);
   ~Preprocessor();

   SatInstance getInstance(std::string const & filename);

   void clear();

   // Mode of operation:
   //

 protected:

   // Helper structures:
   //
   struct ElimLt
   {
      const vec<int>& n_occ;
      explicit ElimLt(const vec<int>& no)
            : n_occ(no)
      {
      }

      // TODO: are 64-bit operations here noticably bad on 32-bit platforms? Could use a saturating
      // 32-bit implementation instead then, but this will have to do for now.
      uint64_t cost(Var const x) const
      {
         return (uint64_t) n_occ[Lit(x, false).toInt()] * (uint64_t) n_occ[(Lit(x, true)).toInt()];
      }
      bool operator()(Var const x, Var const y) const
      {
         return cost(x) < cost(y);
      }
   };

   struct ClauseDeleted
   {
      const ClauseAllocator& ca;
      explicit ClauseDeleted(const ClauseAllocator& _ca)
            : ca(_ca)
      {
      }
      bool operator()(const CRef& cr) const
      {
         return ca[cr].mark() == 1;
      }
   };

   bool elim;          // Perform variable elimination.
   bool ok;

   int verb;
   int grow;  // Allow a variable elimination step to grow by a number of clauses (default to zero).
   int clause_lim;  // Variables are not eliminated if it produces a resolvent with a length above this limit.
   int subsumption_lim;  // Do not check if subsumption against a clause larger than this. -1 means no limit.
   int nRedundantAssignedLits;
   int bwdsub_assigns;
   int n_touched;
   int eliminated_vars;
   CRef bwdsub_tmpunit;
   double simp_garbage_frac;  // A different limit for when to issue a GC during simplification (Also see 'garbage_frac').

   // Statistics:
   //
   EliminatedClauseDatabase elimDb;

   vec<char> touched;
   vec<CRef> clauses;
   OccLists<Var, vec<CRef>, ClauseDeleted> occurs;
   vec<int> n_occ;
   Heap<ElimLt> elim_heap;
   Queue<CRef> subsumption_queue;
   vec<char> eliminated;
   vec<Lit> add_tmp;


   Statistic stat; // currently only for branch interface
   SolveMode smode; // currently only for branch interface
   Random randEngine; // currently only for branch interface

   DratPrint<Lit> drat;
   ClauseAllocator ca;
   ImplicationGraph<ClauseAllocator> ig;
   Branch<ClauseAllocator> branch;
   MinisatPropagate<ClauseAllocator> propEngine;

   // Temporaries:
   //
   vec<Lit> add_oc;

   bool readInstance(std::string const & filename);
   bool readInstance_main(StreamBuffer & in);
   void readInstance_clause(StreamBuffer & in,vec<Lit> & lits);

   bool isEliminated(Var v) const;

   bool hasInterrupt() const
   {
      //FIXME use connector for interrupt
      return false;
   }

   bool eliminate();  // Perform variable elimination based simplification.
   bool eliminate_();

   bool eliminateMinimalNiver();

   int nVars() const
   {
      return ig.nVars();
   }
   int nClauses() const
   {
      return clauses.size();
   }
   int nFreeVars() const;

   int & verbosity()
   {
      return verb;
   }
   const int & verbosity() const
   {
      return verb;
   }

   // Memory managment:
   //
   void garbageCollect(bool const finalGarbage = false);

   // Main internal methods:
   //
   Var newVar();
   bool addEmptyClause();                // Add the empty clause to the solver.
   bool addClause(vec<Lit>& ps, bool const initial = false);
   bool substitute(Var v, Lit x);  // Replace all occurences of v with x (may cause a contradiction).

   void updateElimHeap(Var v);
   void gatherTouchedClauses();
   bool merge(const Clause& _ps, const Clause& _qs, Var v, vec<Lit>& out_clause);
   bool merge(const Clause& _ps, const Clause& _qs, Var v, int& size);
   bool backwardSubsumptionCheck(bool verbose = false);
   bool eliminateVar(Var const v);
   bool eliminateVarMinimal(Var const v);

   bool removeRedundant(bool const removeFalseLits = false);

   void removeClause(CRef const cr);
   bool strengthenClause(CRef cr, Lit l);
   void relocAll(ClauseAllocator& to, bool const finalGarbage = false);

   bool setOk(bool const b);
   bool isOk() const;
};

//=================================================================================================
// Implementation of inline methods:

inline bool Preprocessor::isOk() const
{
   return ok;
}

inline bool Preprocessor::setOk(bool const b)
{
   return ok = b;
}

inline bool Preprocessor::isEliminated(Var v) const
{
   return eliminated[v];
}
inline void Preprocessor::updateElimHeap(Var v)
{
   assert(elim);
   // if (!frozen[v] && !isEliminated(v) && value(v) == l_Undef)
   if (elim_heap.inHeap(v) || (!isEliminated(v) && ig.value(v).isUndef()))
      elim_heap.update(v);
}

inline bool Preprocessor::addEmptyClause()
{
   add_tmp.clear();
   return addClause(add_tmp);
}

//=================================================================================================
}

#endif
