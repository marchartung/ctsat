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

#include "utils/CPUBind.h"
#include <thread>
#include <sstream>
#include <cassert>
#include <fstream>
#include <string>
#include <sstream>
#include <sched.h>
#include <iostream>

namespace ctsat
{

NumaAwareSet NumaAwareSet::instance;

unsigned findNumberAfter(const std::string & pattern, const std::string & str)
{
   std::string::size_type pos = str.find(pattern);
   if (pos == std::string::npos)
      return static_cast<unsigned>(-1);
   pos += pattern.size();
   std::stringstream ss;

   for (; pos < str.size(); ++pos)
      if (std::isdigit(str[pos]))
         break;
   if (pos >= str.size())
      return static_cast<unsigned>(-1);

   for (; pos < str.size(); ++pos)
      if (std::isdigit(str[pos]))
         ss << str[pos];
      else
         break;
   unsigned res;
   ss >> res;
   return res;
}

std::string NumaAwareSet::getProcString() const
{
   std::ifstream f("/proc/cpuinfo", std::ios_base::in | std::ios_base::binary);
   std::stringstream ss;
   char c;
   while (f.good())
   {
      f.get(c);
      if (c != '\n')
         ss << c;
   }
   f.close();
   return ss.str();
}

std::vector<std::string> splitBefore(const std::string & pattern, const std::string & str)
{
   std::vector<std::string> res;
   std::string::size_type pos = str.find(pattern), next;
   if (pos > 0)
   {
      if (pos == std::string::npos)
         return res;
      res.push_back(str.substr(0, pos));
   }
   while ((next = str.find(pattern, pos + pattern.size())) != std::string::npos)
   {
      res.push_back(str.substr(pos, next));
      pos = next;
   }
   if (pos < str.size())
      res.push_back(str.substr(pos, str.size()));
   return res;
}

void NumaAwareSet::initializeFallback()
{
   idealMap.resize(std::thread::hardware_concurrency());
   for (size_t i = 0; i < idealMap.size(); ++i)
      idealMap[i] = i;
   numVCores = idealMap.size();
   numCores = idealMap.size();
   numNumaNodes = 1;
}
struct CpuInfo
{
   unsigned threadId;
   unsigned vcore;
   unsigned coreId;
   unsigned numaId;

   CpuInfo(const std::string & s)
         : threadId(findNumberAfter("processor", s)),
           vcore(0),
           coreId(findNumberAfter("core id", s)),
           numaId(findNumberAfter("physical id", s))
   {
   }
};

void NumaAwareSet::initialize()
{
   std::string file = getProcString();
   std::vector<std::string> split = splitBefore("processor", file);
   bool correct = false;

   if (split.size() > 1)
   {

      unsigned maxVCores = 0;
      std::vector<CpuInfo> cinfos;
      std::vector<unsigned> numDiffCores;
      {
         std::vector<std::vector<unsigned>> coreIdMap;
         std::vector<std::vector<unsigned>> vCoreNum;
         for (std::string const & str : split)
         {
            cinfos.emplace_back(str);
            CpuInfo & i = cinfos.back();
            unsigned const numa = i.numaId;
            unsigned const core = i.coreId;
            if (numa >= coreIdMap.size())
            {
               coreIdMap.resize(numa + 1);
               vCoreNum.resize(numa + 1);
               numDiffCores.resize(numa + 1, 0);
            }
            if (core >= vCoreNum[numa].size())
            {
               vCoreNum[numa].resize(core + 1, 0);
               coreIdMap[numa].resize(core + 1, -1);
            }
            if (vCoreNum[numa][core] == 0)
            {
               coreIdMap[numa][core] = numDiffCores[numa]++;
            }
            // reasign core ids so they are sequential and add vcore numbering:
            i.coreId = coreIdMap[numa][core];
            i.vcore = vCoreNum[numa][core];  // numa and threadid remain

            ++vCoreNum[numa][core];
            maxVCores = (maxVCores < vCoreNum[numa][core]) ? vCoreNum[numa][core] : maxVCores;
         }
      }
      for (size_t i = 1; i < numDiffCores.size(); ++i)
         assert(numDiffCores[i] == numDiffCores[0] && "Heterogeneous processors not supported");
      idealMap.clear();
      idealMap.resize(cinfos.size(), cinfos.size() + 1);

      unsigned const nNuma = numDiffCores.size();
      unsigned const nCores = numDiffCores[0];
      unsigned coresPerVPartition = nCores * nNuma;

      for (size_t i = 0; i < idealMap.size(); ++i)
      {
         CpuInfo & info = cinfos[i];
         unsigned const pos = info.coreId * nNuma + info.numaId + info.vcore * coresPerVPartition;
         idealMap[pos] = info.threadId;
      }
      numCores = nCores * nNuma;
      numVCores = numCores * maxVCores;
      numNumaNodes = nNuma;
      isFallback = false;
      correct = true;
      for (auto const & m : idealMap)
         assert(m >= 0 && static_cast<unsigned>(m) < cinfos.size());
   }
   if (!correct)
      initializeFallback();
//   print();
}

void NumaAwareSet::print()
{
   std::cout << " virtual cores     : " << numVCores << std::endl;
   std::cout << " cores per node    : " << numCores / numNumaNodes << std::endl;
   std::cout << " numa nodes        : " << numNumaNodes << std::endl;
   std::cout << " threads per core  : " << numVCores / numCores << std::endl << std::endl;
}

NumaAwareSet::NumaAwareSet()
      : isFallback(true),
        numVCores(std::thread::hardware_concurrency()),
        numCores(std::thread::hardware_concurrency() / 2),
        numNumaNodes(1)
{
   initialize();
}

std::vector<int> const & NumaAwareSet::getIdealThreadCoreMapping() const
{
   return idealMap;
}

int NumaAwareSet::getNumCores() const
{
   return numCores;
}
int NumaAwareSet::getNumVCores() const
{
   return numVCores;
}
int NumaAwareSet::getNumNumaNodes() const
{
   return numNumaNodes;
}

void CPUBind::bindThread(const int& coreId)
{
#ifndef MACOS
   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(coreId, &cpuset);
   sched_setaffinity(0, sizeof(cpuset), &cpuset);
#endif
}
}

