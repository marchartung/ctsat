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

#include <errno.h>

#include <zlib.h>

#include "core/SolverRunner.h"
#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/ResourceLimits.h"
#include "utils/Options.h"
#include "initial/Inputs.h"
#include "database/BasicTypes.h"
using namespace CTSat;


//=================================================================================================
// Main:

int main(int argc, char** argv)
{
   try
   {
      setUsageHelp(
            "USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
      printf("c This is MapleLCMDistChronoBT.\n");

      // Extra options:
      //

      Inputs::setArgs(&argc, &argv);
      parseOptions(argc, argv, true);

      // Use signal handlers that forcibly quit until the solver will be able to respond to
      // interrupts:

      if (Inputs::cpu_lim != INT32_MAX)
         ResourceLimits::setCpuLimit(Inputs::cpu_lim);

      if (Inputs::mem_lim != INT32_MAX)
         ResourceLimits::setMemLimit(Inputs::mem_lim);

      lbool ret = lbool::Undef();
      if (argc == 1)
      {
         printf("No dimacs passed. Abort\n");
      }
      else
      {
         ret = SolverRunner::run(SolverConfig::getInputConfig());
      }


      return (ret.isTrue() ? 10 : ret.isFalse() ? 20 : 0);
   } catch (OutOfMemoryException&)
   {
      printf("c ===============================================================================\n");
      printf("c Out of memory\n");
      printf("s UNKNOWN\n");
      exit(0);
   } catch (InputException&)
   {
      printf("c ===============================================================================\n");
      printf("c Input error\n");
      printf("s UNKNOWN\n");
      exit(0);
   }

}
