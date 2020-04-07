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

#ifndef SOURCES_DATABASE_MINISATALLOCATORDB_H_
#define SOURCES_DATABASE_MINISATALLOCATORDB_H_

#include "database/BasicTypes.h"
#include "mtl/Vec.h"
#include "mtl/Alloc.h"
#include "mtl/Alg.h"

namespace CTSat
{

//=================================================================================================
// Clause -- a simple class for representing a clause:

class Clause;
typedef RegionAllocator<uint32_t>::Ref CRef;

class Clause
{
   static const unsigned numBitsLbd = 23;
   struct Header
   {
      Header() = delete;
      Header(Header&&) = delete;
      Header& operator=(Header&&) = delete;
      Header& operator=(Header const&) = delete;

      unsigned has_extra :1;
      unsigned reloced :1;
      unsigned removable :1;
      unsigned simplified :1;
      unsigned learnt :2;
      unsigned mark :2;
      unsigned exportMark :2;
      unsigned lbd :numBitsLbd;
      unsigned size :32;
      //simplify
      Header(unsigned const size, bool const use_extra, bool const learnt)
            : has_extra(learnt | use_extra),
              reloced(0),
              removable(1),
              simplified(0),
              learnt(learnt),
              mark(0),
              exportMark(0),
              lbd(0),
              size(size)
      {
      }
      Header(Header const & h)
            : has_extra(h.has_extra),
              reloced(h.reloced),
              removable(h.removable),
              simplified(h.simplified),
              learnt(h.learnt),
              mark(h.mark),
              exportMark(h.exportMark),
              lbd(h.lbd),
              size(h.size)
      {
      }
   };
   Header header;

   union
   {
      Lit lit;
      float act;
      uint32_t abs;
      uint32_t touched;
      CRef rel;
   } data[0];

   friend class ClauseAllocator;

   // NOTE: This constructor cannot be used directly (doesn't allocate enough memory).
   template <class V>
   Clause(const V& ps, bool use_extra, bool learnt)
         : header(ps.size(), use_extra, learnt)
   {
      for (int i = 0; i < ps.size(); i++)
         data[i].lit = ps[i];

      if (header.has_extra)
      {
         if (header.learnt)
         {
            data[header.size].act = 0;
            data[header.size + 1].touched = 0;
         } else
            calcAbstraction();
      }
   }

   Clause(Clause const & c)
         : header(c.header)
   {
      for (int i = 0; i < c.size(); i++)
         data[i].lit = c[i];

      if (header.has_extra)
      {
         if (header.learnt)
         {
            data[header.size].act = 0;
            data[header.size + 1].touched = 0;
         } else
            calcAbstraction();
      }
   }

 public:
   void calcAbstraction()
   {
      assert(header.has_extra);
      uint32_t abstraction = 0;
      for (int i = 0; i < size(); i++)
         abstraction |= 1 << ((data[i].lit.var()) & 31);
      data[header.size].abs = abstraction;
   }

   int size() const;
   void shrink(int i);
   void pop();
   bool learnt() const;
   bool has_extra() const;
   uint32_t mark() const;
   void mark(uint32_t m);

   void setExport(unsigned const m)
   {
      header.exportMark = m;
   }

   unsigned getExport() const
   {
      return header.exportMark;
   }


   const Lit& last() const;

   bool reloced() const;
   CRef relocation() const;
   void relocate(CRef c);

   int lbd() const;
   void set_lbd(int lbd);
   bool removable() const;
   void removable(bool b);

   // NOTE: somewhat unsafe to change the clause in-place! Must manually call 'calcAbstraction' afterwards for
   //       subsumption operations to behave correctly.
   Lit& operator [](int i);
   Lit operator [](int i) const;
   operator const Lit*(void) const;

   Lit* getData()
   {
      return &(data[0].lit);
   }

   uint32_t const& touched() const;
   uint32_t& touched();
   float const & activity() const;
   float& activity();
   uint32_t const & abs() const;

   void setLearnt(unsigned const l)
   {
      assert(l < 4);
      header.learnt = l;
   }

   unsigned getLearnt() const
   {
      return header.learnt;
   }

   Lit subsumes(const Clause& other) const;
   void strengthen(Lit p);
   // simplify
   //
   void setSimplified(bool b);
   bool simplified();
};

inline int Clause::size() const
{
   return header.size;
}
inline void Clause::shrink(int i)
{
   assert(i <= size());
   if (header.has_extra)
      data[header.size - i] = data[header.size];
   header.size -= i;
}
inline void Clause::pop()
{
   shrink(1);
}
inline bool Clause::learnt() const
{
   return header.learnt == 1;
}
inline bool Clause::has_extra() const
{
   return header.has_extra;
}
inline uint32_t Clause::mark() const
{
   return header.mark;
}
inline void Clause::mark(uint32_t m)
{
   header.mark = m;
}

inline const Lit& Clause::last() const
{
   return data[header.size - 1].lit;
}

inline bool Clause::reloced() const
{
   return header.reloced;
}
inline CRef Clause::relocation() const
{
   return data[0].rel;
}
inline void Clause::relocate(CRef c)
{
   header.reloced = 1;
   data[0].rel = c;
}

inline int Clause::lbd() const
{
   return header.lbd;
}
inline void Clause::set_lbd(int lbd)
{
   static const int maxLbd = 1 << (numBitsLbd - 1);
   header.lbd = (lbd > maxLbd) ? maxLbd : lbd;
}
inline bool Clause::removable() const
{
   return header.removable;
}
inline void Clause::removable(bool b)
{
   header.removable = b;
}

// NOTE: somewhat unsafe to change the clause in-place! Must manually call 'calcAbstraction' afterwards for
//       subsumption operations to behave correctly.
inline Lit& Clause::operator [](int i)
{
   return data[i].lit;
}
inline Lit Clause::operator [](int i) const
{
   return data[i].lit;
}
inline Clause::operator const Lit*(void) const
{
   return (Lit*) data;
}

inline uint32_t& Clause::touched()
{
   assert(header.has_extra && header.learnt);
   return data[header.size + 1].touched;
}

inline uint32_t const & Clause::touched() const
{
   assert(header.has_extra && header.learnt);
   return data[header.size + 1].touched;
}

inline float const& Clause::activity() const
{
   assert(header.has_extra);
   return data[header.size].act;
}

inline float& Clause::activity()
{
   assert(header.has_extra);
   return data[header.size].act;
}
inline uint32_t const & Clause::abs() const
{
   assert(header.has_extra);
   return data[header.size].abs;
}
// simplify
//
inline void Clause::setSimplified(bool b)
{
   header.simplified = b;
}
inline bool Clause::simplified()
{
   return header.simplified;
}

//=================================================================================================
// ClauseAllocator -- a simple class for allocating memory for clauses:

const CRef CRef_Undef = RegionAllocator<uint32_t>::Ref_Undef;
class ClauseAllocator : public RegionAllocator<uint32_t>
{
   static int clauseWord32Size(int size, int extras)
   {
      return (sizeof(Clause) + (sizeof(Lit) * (size + extras))) / sizeof(uint32_t);
   }
 public:
   typedef uint32_t base_type;
   bool extra_clause_field;

   typedef CTSat::Lit Lit;
   typedef CTSat::CRef CRef;
   typedef CTSat::Var Var;
   typedef CTSat::lbool lbool;
   typedef CTSat::Clause Clause;

   inline static CRef npos()
   {
      return CRef_Undef;
   }

   ClauseAllocator(uint32_t start_cap)
         : RegionAllocator<uint32_t>(start_cap),
           extra_clause_field(false)
   {
   }

   ClauseAllocator(ClauseAllocator const & ca)
         : RegionAllocator<uint32_t>(ca),
           extra_clause_field(ca.extra_clause_field)
   {

   }

   ClauseAllocator(void const * data, uint64_t const nBytes)
         : RegionAllocator<uint32_t>(data, nBytes),
           extra_clause_field(false)
   {
   }

   ClauseAllocator(ClauseAllocator && ca)
         : RegionAllocator<uint32_t>(std::move(ca)),
           extra_clause_field(ca.extra_clause_field)
   {
   }

   ClauseAllocator()
         : RegionAllocator<uint32_t>(1024 * 1024),
           extra_clause_field(false)
   {
   }

   ClauseAllocator & operator=(ClauseAllocator const & ca)
   {
      extra_clause_field = ca.extra_clause_field;
      *reinterpret_cast<RegionAllocator<uint32_t>*>(this) = ca;
      return *this;
   }

   ClauseAllocator & operator=(ClauseAllocator && ca)
   {
      std::swap(extra_clause_field, ca.extra_clause_field);
      *reinterpret_cast<RegionAllocator<uint32_t>*>(this) = std::move(ca);
      return *this;
   }

   void moveTo(ClauseAllocator& to)
   {
      to.extra_clause_field = extra_clause_field;
      RegionAllocator<uint32_t>::moveTo(to);
   }

   template <class Lits>
   CRef alloc(const Lits& ps, bool learnt)
   {
      assert(sizeof(Lit) == sizeof(uint32_t));
      assert(sizeof(float) == sizeof(uint32_t));
      int extras = learnt ? 2 : (int) extra_clause_field;
      CRef cid = RegionAllocator<uint32_t>::alloc(clauseWord32Size(ps.size(), extras));
      new (lea(cid)) Clause(ps, extra_clause_field, learnt);
      return cid;
   }

   inline CRef alloc(Clause const & c)
   {
      int extras = (c.getLearnt() > 0) ? 2 : (int) extra_clause_field;
      CRef cid = RegionAllocator<uint32_t>::alloc(clauseWord32Size(c.size(), extras));
      new (lea(cid)) Clause(c);
      return cid;
   }

   inline CRef next(CRef const ref) const
   {
      Clause const & c = operator[](ref);
      CRef res = ref
         + clauseWord32Size(c.size(), (c.getLearnt() != 0) ? 2 : (int) extra_clause_field);
      if (res >= size())
         res = CRef_Undef;
      return res;
   }

   // Deref, Load Effective Address (LEA), Inverse of LEA (AEL):
   inline Clause& operator[](Ref r)
   {
      return (Clause&) RegionAllocator<uint32_t>::operator[](r);
   }
   inline const Clause& operator[](Ref r) const
   {
      return (Clause&) RegionAllocator<uint32_t>::operator[](r);
   }
   Clause* lea(Ref r)
   {
      return (Clause*) RegionAllocator<uint32_t>::lea(r);
   }
   const Clause* lea(Ref r) const
   {
      return (Clause*) RegionAllocator<uint32_t>::lea(r);
   }
   Ref ael(const Clause* t)
   {
      return RegionAllocator<uint32_t>::ael((uint32_t*) t);
   }

   void free(CRef cid)
   {
      Clause& c = operator[](cid);
      int extras = c.learnt() ? 2 : (int) c.has_extra();
      RegionAllocator<uint32_t>::free(clauseWord32Size(c.size(), extras));
   }

   void reloc(CRef& cr, ClauseAllocator& to)
   {
      Clause& c = operator[](cr);
      if (!c.reloced())
      {
         cr = to.alloc(c);
         c.relocate(cr);
      }
      cr = c.relocation();
   }
};

inline std::ostream& operator<<(std::ostream& out, const Clause& cls)
{
   for (int i = 0; i < cls.size(); ++i)
   {
      out << cls[i] << " ";
   }

   return out;
}

/*_________________________________________________________________________________________________
 |
 |  subsumes : (other : const Clause&)  ->  Lit
 |
 |  Description:
 |       Checks if clause subsumes 'other', and at the same time, if it can be used to simplify 'other'
 |       by subsumption resolution.
 |
 |    Result:
 |       lit_Error  - No subsumption or simplification
 |       lit_Undef  - Clause subsumes 'other'
 |       p          - The literal p can be deleted from 'other'
 |________________________________________________________________________________________________@*/
inline Lit Clause::subsumes(const Clause& other) const
{
   //if (other.size() < size() || (extra.abst & ~other.extra.abst) != 0)
   //if (other.size() < size() || (!learnt() && !other.learnt() && (extra.abst & ~other.extra.abst) != 0))
   assert(!learnt());
   assert(!other.learnt());
   assert(has_extra());
   assert(other.has_extra());
   if (other.size() < size()
      || (abs() & ~other.abs()) != 0)
      return Lit::Error();

   Lit ret = Lit::Undef();
   const Lit* c = (const Lit*) (*this);
   const Lit* d = (const Lit*) other;

   for (int i = 0; i < size(); i++)
   {
      // search for c[i] or ~c[i]
      for (int j = 0; j < other.size(); j++)
         if (c[i] == d[j])
            goto ok;
         else if (ret == Lit::Undef() && c[i] == ~d[j])
         {
            ret = c[i];
            goto ok;
         }

      // did not find it
      return Lit::Error();
      ok: ;
   }

   return ret;
}

inline void Clause::strengthen(Lit p)
{
   remove(*this, p);
   calcAbstraction();
}

}

#endif /* SOURCES_DATABASE_MINISATALLOCATORDB_H_ */
