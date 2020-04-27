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

#ifndef UTILS_LOGGER_H_
#define UTILS_LOGGER_H_

#include <pthread.h>
#include <fstream>
#include <string>
#include <cassert>
#include <iostream>

#include "utils/Exceptions.h"
#include "initial/Inputs.h"

namespace ctsat
{

class Logger
{
 public:
   static thread_local Logger instance;
   inline Logger()
         : active(false)
   {

   }

   inline void init(int const id)
   {
      if (Inputs::log)
      {
         pthread_t pid = pthread_self();
         std::string const filename = std::string(Inputs::logDirectory)
            + std::to_string(id)
            + "_"
            + std::to_string(pid)
            + ".log";
         outs.open(filename, std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
         if (!outs.good())
         {
            std::cout << "Error: Could not open log file" << std::endl;
            throw InputException();
         }
         active = true;
         log("logger initialized");
      }
   }
   inline void deinit()
   {
      if (active)
      {
         log("logger shut down");
         assert(outs.good());
         outs.close();
         active = false;
      }
   }

   inline void log(std::string const & str)
   {
      if (active)
      {
         assert(outs.good());
         outs << str << std::endl;
      }
   }

 private:
   bool active;
   std::ofstream outs;
};


}
#ifdef LOG_CTSAT_ENABLE

#define LOG_INIT(id) ctsat::Logger::instance.init(id);
#define LOG_DEINIT ctsat::Logger::instance.deinit();
#define LOG(str) ctsat::Logger::instance.log(str);

#else
inline void checkConsistent()
{
   if(ctsat::Inputs::log)
      std::cout << "Warning: Logging was not build. Use -DLOG_CTSAT_ENABLE as compile flag" << std::endl;
}

#define LOG_INIT(id) checkConsistent();
#define LOG_DEINIT
#define LOG(str)
#endif

#endif /* UTILS_LOGGER_H_ */
