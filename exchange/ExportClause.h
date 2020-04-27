/*
 * ExportClause.h
 *
 *  Created on: 26.04.2020
 *      Author: hartung
 */

#ifndef EXCHANGE_EXPORTCLAUSE_H_
#define EXCHANGE_EXPORTCLAUSE_H_


#include <cstdint>

namespace ctsat
{
template <typename Database>
struct ExportClause
{
   typedef typename Database::Clause Clause;
   typedef typename Database::Lit Lit;

   ExportClause() = delete;
   ExportClause(ExportClause<Database> &&) = delete;
   ExportClause<Database>& operator=(ExportClause<Database> const&) = delete;
   ExportClause<Database>& operator=(ExportClause<Database> &&) = delete;

   uint16_t id;
   uint16_t lbd;
   int const sz;
   Lit clause[1];
   ExportClause(Clause const & c, unsigned const lbd, unsigned const id);
   ExportClause(ExportClause<Database> const & in);
   ExportClause(Lit const & l, unsigned const id);
   static uint64_t nbytes(Clause const & c);
   static uint64_t nbytes(Lit const & l);
   static uint64_t nbytes(ExportClause<Database> const & c);
   int size() const
   {
      return sz;
   }
   Lit const & operator[](int const idx) const
   {
      return clause[idx];
   }

   void set(Clause & c, bool const minimizeClause) const
   {
      c.set_lbd(lbd);
      c.setSimplified(!minimizeClause);
   }
};

}




#endif /* EXCHANGE_EXPORTCLAUSE_H_ */
