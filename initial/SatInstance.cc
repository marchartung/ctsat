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

#include "initial/SatInstance.h"

namespace CTSat
{

SatInstance::SatInstance()
{
}

SatInstance::SatInstance(SatInstance&& in)
{
   *this = std::move(in);
}

SatInstance::SatInstance(
                         vec<lbool> && model,
                         vec<CRef> && clauses,
                         Database&& db,
                         EliminatedClauseDatabase&& elimDb,
                         DratPrint<Lit> && drat)
      : model(std::move(model)),
        clauses(std::move(clauses)),
        ca(std::move(db)),
        elimDb(std::move(elimDb)),
        drat(std::move(drat))
{

}

SatInstance::SatInstance(void const * data, uint64_t const nBytes)
      : ca(data, nBytes)
{
   // read in clauses from stream, catch max var, collect clause refs, set model to highest var
   Var maxVar = 0;
   CRef pos = 0;
   while (pos != Database::npos())
   {
      clauses.push(pos);
      Clause const & c = ca[pos];
      assert(!c.learnt());
      for (int i = 0; i < c.size(); ++i)
         if (c[i].var() > maxVar)
            maxVar = c[i].var();
      pos = ca.next(pos);
   }
   model.growTo(maxVar + 1, lbool::Undef());
}

SatInstance& SatInstance::operator =(SatInstance&& in)
{
   model = std::move(in.model);
   clauses = std::move(in.clauses);
   ca = std::move(in.ca);
   elimDb = std::move(in.elimDb);
   return *this;
}

bool SatInstance::isOk() const
{
   return model.size() > 0;
}

bool SatInstance::isClean() const
{
   if (ca.wasted() > 0)
      return false;
   for (int i = 0; i < clauses.size(); ++i)
   {
      Clause const & c = ca[clauses[i]];
      for (int j = 0; j < c.size(); ++j)
         if (!model[c[j].var()].isUndef())
            return false;
   }
   return true;
}

bool SatInstance::operator==(SatInstance const & in) const
{
   if (in.clauses.size() != clauses.size())
      return false;
   for (int i = 0; i < clauses.size(); ++i)
   {
      Clause const & c1 = ca[clauses[i]];
      Clause const & c2 = in.ca[clauses[i]];
      if (c1.size() != c2.size())
         return false;
      for (int j = 0; j < c1.size(); ++j)
      {
         if (c1[j] != c2[j])
            return false;

         if(c1[j].var() >= model.size())
            return false;
         if(c1[j].var() >= in.model.size())
            return false;
      }
   }
   for(int i=0;i<clauses.size();++i)
   {
      bool found = false;
      for(int j=0;j<in.clauses.size();++j)
         if(clauses[i] == in.clauses[j])
         {
            found = true;
            break;
         }
      if(!found)
         return false;
   }

   return true;
}

}
