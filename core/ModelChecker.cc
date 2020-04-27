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

#include "core/ModelChecker.h"

namespace ctsat
{
bool ModelChecker::checkSat(vec<lbool> const & model, std::string const filename)
{
   gzFile f = gzopen(filename.c_str(), "rb");
   if (f == NULL)
      printf("c ERROR! Could not open file: %s\n", filename.c_str()), throw InputException();
   StreamBuffer in(f);
   bool const res = checkSat(model, in);
   gzclose(f);
   return res;
}

void ModelChecker::printSatisfiedClauses(vec<lbool> const & model, std::string const filename)
{
   gzFile f = gzopen(filename.c_str(), "rb");
   if (f == NULL)
      printf("c ERROR! Could not open file: %s\n", filename.c_str()), throw InputException();
   StreamBuffer in(f);
   printSatisfiedClauses(model, in);
   gzclose(f);
}

void ModelChecker::printUndefClauses(vec<lbool> const & model, std::string const filename)
{
   gzFile f = gzopen(filename.c_str(), "rb");
   if (f == NULL)
      printf("c ERROR! Could not open file: %s\n", filename.c_str()), throw InputException();
   StreamBuffer in(f);
   printUndefClauses(model, in);
   gzclose(f);
}

void ModelChecker::printClause(vec<lbool> const & model, vec<Lit> & c)
{
   std::cout << "c [ ";
   for (int i = 0; i < c.size(); ++i)
   {
      lbool const mval = model[c[i].var()];
      lbool val = lbool::Undef();
      if (!mval.isUndef())
         val = (mval.sign() == c[i].sign()) ? lbool::True() : lbool::False();

      std::cout << ((c[i].sign()) ? "-" : "") << c[i].var() << "("
                << ((val.isUndef()) ? "u)" : ((val.isTrue()) ? "t)" : "f)")) << " ";
   }
   std::cout << "]\n";
}

void ModelChecker::printUndefClauses(vec<lbool> const & model, StreamBuffer & in)
{
   vec<Lit> lits;
//         int vars = 0;
//   int clauses = 0;
   int cnt = 0;
   for (;;)
   {
      skipWhitespace(in);
      if (*in == EOF)
         break;
      else if (*in == 'p')
      {
         if (eagerMatch(in, "p cnf"))
         {
//                  vars = parseInt(in);
            parseInt(in);
            /*clauses = */parseInt(in);
         } else
         {
            printf("c PARSE ERROR! Unexpected char: %c\n", *in), throw InputException();
         }
      } else if (*in == 'c' || *in == 'p')
         skipLine(in);
      else
      {
         cnt++;
         int parsed_lit;
         bool isSat = false, isUndef = false;
         for (;;)
         {
            parsed_lit = parseInt(in);
            if (parsed_lit == 0)
               break;  //{printf("\n"); break;}
            int const var = abs(parsed_lit) - 1;
            // printf("%d ", parsed_lit);
            if ((parsed_lit > 0 && model[var].isTrue()) || (parsed_lit < 0 && model[var].isFalse()))
               isSat = true;
            else if (model[var].isUndef())
               isUndef = true;
         }
         if (!isSat && isUndef)
         {
            printf("c clause %d is undefined\n", cnt);
            // break;
         }
      }
   }
}

bool ModelChecker::printSatisfiedClauses(vec<lbool> const & model, StreamBuffer & in)
{
   vec<Lit> lits;
//         int vars = 0;
//   int clauses = 0;
   int cnt = 0;
   bool ok = true;
   for (;;)
   {
      skipWhitespace(in);
      if (*in == EOF)
         break;
      else if (*in == 'p')
      {
         if (eagerMatch(in, "p cnf"))
         {
//                  vars = parseInt(in);
            parseInt(in);
            /* clauses =*/parseInt(in);
         } else
         {
            printf("c PARSE ERROR! Unexpected char: %c\n", *in), throw InputException();
         }
      } else if (*in == 'c' || *in == 'p')
         skipLine(in);
      else
      {
         cnt++;
         int parsed_lit, var;
         bool isSat = false;
         for (;;)
         {
            parsed_lit = parseInt(in);
            if (parsed_lit == 0)
               break;  //{printf("\n"); break;}
            var = abs(parsed_lit) - 1;
            // printf("%d ", parsed_lit);
            if ((parsed_lit > 0 && model[var].isTrue()) || (parsed_lit < 0 && model[var].isFalse()))
               isSat = true;
         }
         if (isSat)
         {
            printf("c clause %d is satisfied\n", cnt);
            ok = false;
            // break;
         }
      }
   }
   return ok;
}

bool ModelChecker::checkSat(vec<lbool> const & model, StreamBuffer & in)
{
   vec<Lit> lits;
//      int vars = 0;
//   int clauses = 0;
   int cnt = 0;
   bool ok = true;
   for (;;)
   {
      skipWhitespace(in);
      if (*in == EOF)
         break;
      else if (*in == 'p')
      {
         if (eagerMatch(in, "p cnf"))
         {
//               vars = parseInt(in);
            parseInt(in);
            /*clauses = */parseInt(in);
         } else
         {
            printf("c PARSE ERROR! Unexpected char: %c\n", *in), throw InputException();
         }
      } else if (*in == 'c' || *in == 'p')
         skipLine(in);
      else
      {
         cnt++;
         int parsed_lit, var;
         bool isSat = false;
         for (;;)
         {
            parsed_lit = parseInt(in);
            if (parsed_lit == 0)
            {
               break;  //{printf("\n"); break;}
            }
            var = abs(parsed_lit) - 1;
            lits.push(Lit(var, parsed_lit < 0));
            assert(!model[var].isUndef());
            // printf("%d ", parsed_lit);
            if ((parsed_lit > 0 && model[var].isTrue()) || (parsed_lit < 0 && model[var].isFalse()))
               isSat = true;
         }
         if (!isSat)
         {
            printf("c clause %d is not satisfied\n", cnt);
            printClause(model, lits);
            ok = false;
            // break;
         }
         lits.clear();
      }
   }
   return ok;
}
}

