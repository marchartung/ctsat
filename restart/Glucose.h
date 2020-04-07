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

#ifndef SOURCES_RESTART_GLUCOSE_H_
#define SOURCES_RESTART_GLUCOSE_H_

#include "mtl/AvgQueue.h"
#include "initial/SolverConfig.h"
#include "core/SolveMode.h"
#include "core/Statistic.h"

namespace CTSat
{

class GlucoseRestart
{
 public:

   GlucoseRestart(SolverConfig const & config, SolveMode & smode, Statistic & stat);

   void notifyRestart();
   void notifyConflictFound();
   void notifyConflictResolved(int const lbd);
   bool shouldRestart();

 private:
   bool cached;
   AvgQueue<int> lbd_queue;   // For computing moving averages of recent LBD values.
   SolveMode & smode;
   Statistic & stat;

   inline static double luby(double y, int x)
   {

      // Find the finite subsequence that contains index 'x', and the
      // size of that subsequence:
      int size, seq;
      for (size = 1, seq = 0; size < x + 1; seq++, size = 2 * size + 1)
         ;

      while (size - 1 != x)
      {
         size = (size - 1) >> 1;
         seq--;
         x = x % size;
      }

      return pow(y, seq);
   }
};

inline GlucoseRestart::GlucoseRestart(const SolverConfig& config, SolveMode& smode, Statistic& stat)
      : cached(false),
        lbd_queue(config.lbd_queue),
        smode(smode),
        stat(stat)
{
}

inline void GlucoseRestart::notifyRestart()
{
   cached = false;
   lbd_queue.clear();
}

inline void GlucoseRestart::notifyConflictResolved(const int lbd)
{

   cached = false;
   lbd_queue.push(lbd);
   stat.global_lbd_sum += (lbd > 50 ? 50 : lbd);
}

inline void GlucoseRestart::notifyConflictFound()
{
}

inline bool GlucoseRestart::shouldRestart()
{
   bool restart = false;
   if (!cached)
   {
      restart = lbd_queue.full() && (lbd_queue.avg() * 0.8 > stat.global_lbd_sum / stat.conflicts);
      cached = true;
   }
   return restart;
}

}

#endif /* SOURCES_RESTART_GLUCOSE_H_ */
