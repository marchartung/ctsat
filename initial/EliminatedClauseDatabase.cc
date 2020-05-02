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

#include "EliminatedClauseDatabase.h"

namespace ctsat
{

void EliminatedClauseDatabase::printModel(vec<lbool> & model) const
{
   printf("v ");
   for (int i = 0; i < model.size(); i++)
      if (!model[i].isUndef())
         printf("%s%s%d", (i == 0) ? "" : " ", (model[i].isTrue()) ? "" : "-", i + 1);
   printf(" 0\n");
}

void EliminatedClauseDatabase::addElimClause(Var const v, const Clause& c)
{
   {
      int first = elimclauses.size();
      int v_pos = -1;

// Copy clause to elimclauses-vector. Remember position where the
// variable 'v' occurs:
      for (int i = 0; i < c.size(); i++)
      {
         elimclauses.push(c[i].toInt());
         if (c[i].var() == v)
            v_pos = i + first;
      }
      assert(v_pos != -1);

// Swap the first literal with the 'v' literal, so that the literal
// containing 'v' will occur first in the clause:
      uint32_t tmp = elimclauses[v_pos];
      elimclauses[v_pos] = elimclauses[first];
      elimclauses[first] = tmp;

// Store the length of the clause last:
      elimclauses.push(c.size());
   }
}

void EliminatedClauseDatabase::extendModel(vec<lbool> & model) const
{
   Lit x;
   int i, j;
   for (i = elimclauses.size() - 1; i > 0; i -= j)
   {
      for (j = elimclauses[i--]; j > 1; j--, i--)
      {
         Lit const l = Lit::toLit(elimclauses[i]);
         lbool const val = model[l.var()];
         assert(!val.isUndef());
         if (l.sign() == val.sign())
            goto next;
      }

      x = Lit::toLit(elimclauses[i]);
      //assert(model[x.var()].isUndef() || model[x.var()] == lbool(!x.sign()));
      model[x.var()] = lbool(!x.sign());
      next: ;
   }
}

}
