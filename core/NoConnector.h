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

#ifndef SOURCES_PARALLEL_NOCONNECTOR_H_
#define SOURCES_PARALLEL_NOCONNECTOR_H_

#include <cassert>
#include <cstdint>
#include <atomic>
#include <vector>
#include <unistd.h>
#include <signal.h>

#include "mtl/Vec.h"
#include "database/BasicTypes.h"

namespace CTSat
{

class NoConnector
{
   NoConnector(NoConnector const &) = delete;
   NoConnector & operator=(NoConnector const &) = delete;
   NoConnector(NoConnector &&) = delete;
   NoConnector & operator=(NoConnector &&) = delete;

 public:
   typedef uint32_t size_type;

   static void signalAbort(int signum);

   NoConnector();

   ~NoConnector();

   void notifyThreadInitialized();
   void waitInitialize(int const nThreads);
   int nInitialzedThreads();

   void notifyThreadStart();
   void notifyThreadEnd();
   int nRunningThreads();

   void allowMemFree();

   void abort();
   bool isAborted() const;
   bool setFinished(lbool const res);
   bool isFinished() const;

   lbool getResult() const;

   void checkBufferState(size_type const pos)
   {
   }

   bool isValid(size_type const & pos) const;
   size_type next(size_type const & pos) const;

   vec<lbool> & getModel();
   void commitModel(vec<lbool> && m);

   template <typename B>
   void commitModel(vec<B> const & min)
   {
      // only translates lbool or copies it
      vec<lbool> m(min.size(), lbool::Undef());
      for (int i = 0; i < min.size(); ++i)
         if (min[i].isFalse())
            m[i] = lbool::False();
         else if (min[i].isTrue())
            m[i] = lbool::True();
      commitModel(std::move(m));
   }

   template <typename ClauseType>
   inline ClauseType const & get(size_type const & pos) const
   {
      assert(false);
      ClauseType const * res = nullptr;
      return *res;
   }

   template <typename ClauseType, typename ... Args>
   inline size_type exchange(uint64_t const & nBytes, Args ... args)
   {
      return 0;
   }

   unsigned getUniqueId()
   {
      return idCounter.fetch_add(1);
   }

   bool shouldImport(size_type const pos) const
   {
      return false;
   }

   void sleep()
   {
      usleep(uSleepTime);
   }

 protected:
   static const unsigned uSleepTime = 50;
   static std::atomic<unsigned> idCounter;
   static std::vector<NoConnector*> abortConnectors;

   std::atomic<bool> modelCommited;
   std::atomic<bool> allowedToFreeMem;
   std::atomic<int> result;
   std::atomic<int> numRunningThreads;
   std::atomic<int> numInitializedThreads;
   vec<lbool> model;

   struct Result
   {

      static constexpr int Undef()
      {
         return 0;
      }
      static constexpr int Abort()
      {
         return 1;
      }
      static constexpr int False()
      {
         return 2;
      }
      static constexpr int True()
      {
         return 3;
      }

      static int getInt(lbool const in)
      {
         int res = Undef();
         if (in.isTrue())
            res = True();
         else if (in.isFalse())
            res = False();
         return res;
      }

      static lbool getLbool(int const in)
      {
         switch (in)
         {
            case True():
               return lbool::True();
            case False():
               return lbool::False();
            default:
               return lbool::Undef();
         }
      }
   };
};

inline void NoConnector::signalAbort(int signum)
{
   for (size_t i = 0; i < abortConnectors.size(); ++i)
      abortConnectors[i]->abort();
}

inline NoConnector::NoConnector()
      : modelCommited(false),
        allowedToFreeMem(false),
        result(Result::Undef()),
        numRunningThreads(0),
        numInitializedThreads(0)
{
   if (abortConnectors.empty())
   {
      signal(SIGINT, NoConnector::signalAbort);
      signal(SIGXCPU, NoConnector::signalAbort);

   }
   abortConnectors.push_back(this);
}

inline NoConnector::~NoConnector()
{
}

inline void NoConnector::waitInitialize(int const nThreads)
{
   while (nThreads > nInitialzedThreads())
      usleep(uSleepTime);
}

inline void NoConnector::notifyThreadStart()
{
   numRunningThreads.fetch_add(1);
}

inline void NoConnector::notifyThreadEnd()
{
   numRunningThreads.fetch_sub(1);
   while (!allowedToFreeMem)
      usleep(uSleepTime);
}

inline void NoConnector::allowMemFree()
{
   allowedToFreeMem = true;
}

inline void NoConnector::notifyThreadInitialized()
{
   numInitializedThreads.fetch_add(1);
}

inline int NoConnector::nRunningThreads()
{
   return numRunningThreads;
}

inline int NoConnector::nInitialzedThreads()
{
   return numInitializedThreads;
}

inline bool NoConnector::isAborted() const
{
   return Result::Abort() == result;
}

inline void NoConnector::commitModel(vec<lbool> && m)
{
   assert(model.size() == 0);
   model = std::move(m);
   assert(!modelCommited);
   modelCommited = true;
}

inline void NoConnector::abort()
{
   int assumed = Result::Undef();
   result.compare_exchange_strong(assumed, Result::Abort());
}

inline bool NoConnector::setFinished(lbool const res)
{
   int assumed = Result::Undef();
   return result.compare_exchange_strong(assumed, Result::getInt(res));
}
inline bool NoConnector::isFinished() const
{
   return result != Result::Undef();
}

inline lbool NoConnector::getResult() const
{
   return Result::getLbool(result.load());
}

inline bool NoConnector::isValid(size_type const & pos) const
{
   return false;
}

inline NoConnector::size_type NoConnector::next(size_type const & pos) const
{
   return 0;
}

inline vec<lbool> & NoConnector::getModel()
{
   assert(getResult().isTrue());
   while (!modelCommited)
      usleep(uSleepTime);
   return model;
}

}

#endif /* SOURCES_PARALLEL_NOCONNECTOR_H_ */
