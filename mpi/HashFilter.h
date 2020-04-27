/*
 * HashFilter.h
 *
 *  Created on: 23.04.2020
 *      Author: hartung
 */

#ifndef MINIMIZE_HASHFILTER_H_
#define MINIMIZE_HASHFILTER_H_

#include <unordered_map>
#include "mtl/Vec.h"
#include "mtl/Sort.h"
namespace ctsat
{

class HashFilter
{
 public:

   template <typename ClauseType>
   bool add(ClauseType const & c)
   {
      size_t const hash = generateHash(c);
      auto it = cHashes.find(hash);
      if(it == cHashes.end())
      {
         cHashes.insert(std::make_pair(hash,age));
         return false;
      }
      else
      {
         it->second = age;
         return true;
      }
   }

   void cleanUp(unsigned olderThan)
   {
      auto it = cHashes.begin();
      while(it != cHashes.end())
      {
         if(age - it->second > olderThan)
            it = cHashes.erase(it);
         else
            ++it;
      }
   }

   void increaseAge()
   {
      ++age;
   }

 private:
   unsigned age;
   vec<int> tmp;
   std::unordered_map<size_t, unsigned> cHashes;

   inline size_t getHashOfTmp() const
   {
      std::hash<int> hasher;
      size_t result = tmp.size();
      for (int i = 0; i < tmp.size(); ++i)
      {
         result = (result << 1) ^ hasher(tmp[i]);
      }
      return result;
   }

   template <typename ClauseType>
   inline size_t generateHash(ClauseType const & c)
   {
      tmp.clear();
      tmp.growTo(c.size());
      for (int i = 0; i < c.size(); ++i)
         tmp[i] = c[i].toInt();
      sort(tmp);
      return getHashOfTmp();
   }

};
}

#endif /* MINIMIZE_HASHFILTER_H_ */
