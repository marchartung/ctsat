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

#ifndef SOURCES_INITIAL_SATINSTANCE_H_
#define SOURCES_INITIAL_SATINSTANCE_H_

#include "mtl/Vec.h"
#include "utils/DratPrint.h"

#include <cassert>
#include <algorithm>

#include "initial/EliminatedClauseDatabase.h"
#include "database/MinisatAllocatorDb.h"

namespace ctsat
{
class SatInstance
{

   SatInstance(SatInstance const &) = delete;
   SatInstance & operator=(SatInstance const &) = delete;

 public:
   typedef ClauseAllocator Database;
   typedef Database::Lit Lit;
   typedef Database::Var Var;
   typedef Database::Clause Clause;
   typedef Database::CRef CRef;

   SatInstance();
   SatInstance(SatInstance && in);
   SatInstance(
               vec<lbool> && model,
               vec<CRef> && clauses,
               Database && db,
               EliminatedClauseDatabase && elimDb,
               DratPrint<Lit> && drat);

   SatInstance(void const * data, uint64_t const nBytes);

   SatInstance & operator=(SatInstance&& in);

   bool operator==(SatInstance const & in) const;

   bool isOk() const;

   bool isClean() const;

   vec<lbool> model;
   vec<CRef> clauses;
   Database ca;
   EliminatedClauseDatabase elimDb;
   DratPrint<Lit> drat;

};
}

#endif /* SOURCES_INITIAL_SATINSTANCE_H_ */
