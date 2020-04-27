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

#ifndef SOURCES_UTILS_RESOURCELIMITS_H_
#define SOURCES_UTILS_RESOURCELIMITS_H_

#include <cstdint>

namespace ctsat
{

class ResourceLimits
{
 public:
   static void setMemLimit(uint64_t const mblim)
   {
      rlim_t new_mem_lim = (rlim_t) mblim * 1024 * 1024;
      rlimit rl;
      getrlimit(RLIMIT_AS, &rl);
      if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max)
      {
         rl.rlim_cur = new_mem_lim;
         if (setrlimit(RLIMIT_AS, &rl) == -1)
            printf("c WARNING! Could not set resource limit: Virtual memory.\n");
      }
   }
   static void setCpuLimit(uint64_t const lim)
   {
      rlimit rl;
      getrlimit(RLIMIT_CPU, &rl);
      if (rl.rlim_max == RLIM_INFINITY || (rlim_t) lim < rl.rlim_max)
      {
         rl.rlim_cur = lim;
         if (setrlimit(RLIMIT_CPU, &rl) == -1)
            printf("c WARNING! Could not set resource limit: CPU-time.\n");
      }

   }
};


}



#endif /* SOURCES_UTILS_RESOURCELIMITS_H_ */
