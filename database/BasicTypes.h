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

#ifndef SOURCES_DATABASE_BASICTYPES_H_
#define SOURCES_DATABASE_BASICTYPES_H_

#include <iostream>
#include <cstdint>

namespace ctsat {

//=================================================================================================
// Variables, literals, lifted booleans, clauses:

// NOTE! Variables are just integers. No abstraction here. They should be chosen from 0..N,
// so that they can be used as array indices.

typedef int Var;
#define var_Undef (-1)
//=================================================================================================
// Lifted booleans:
//
// NOTE: this implementation is optimized for the case when comparisons between values are mostly
//       between one variable and one constant. Some care had to be taken to make sure that gcc
//       does enough constant propagation to produce sensible code, and this appears to be somewhat
//       fragile unfortunately.


class lbool {
	uint8_t value;

public:

	static constexpr lbool Undef();
	static constexpr lbool True();
	static constexpr lbool False();

	bool isTrue() const;
   bool isFalse() const;
   bool isUndef() const;

   bool sign() const;

	explicit constexpr lbool(uint8_t const & v);
	constexpr lbool();
	explicit constexpr lbool(bool const & x);

	lbool & operator=(lbool const & in);

	bool constexpr operator ==(lbool const & b) const;
	bool constexpr operator !=(lbool const & b) const;
	lbool constexpr operator ^(bool const & b) const;
	lbool constexpr operator &&(lbool const & b) const;
	lbool constexpr operator ||(lbool const & b) const;

	int constexpr toInt() const;
	static constexpr lbool toLbool(int const v);
};


inline lbool & lbool::operator=(lbool const & in)
{
   value = in.value;
   return *this;
}


inline bool lbool::isTrue() const
{
   return *this == lbool::True();
}
inline bool lbool::isFalse() const
{
   return *this == lbool::False();
}
inline bool lbool::isUndef() const
{
   return *this == lbool::Undef();
}
inline bool lbool::sign() const
{
   return isFalse();
}

inline constexpr lbool lbool::Undef() {
	return lbool((uint8_t)2);
}

inline constexpr lbool lbool::True() {
	return lbool((uint8_t)0);
}
constexpr lbool lbool::False() {
	return lbool((uint8_t)1);
}

inline constexpr lbool::lbool(uint8_t const & v) :
		value(v) {
}

inline constexpr lbool::lbool() :
		value(0) {
}
inline constexpr lbool::lbool(bool const & x) :
		value(!x) {
}

inline bool constexpr lbool::operator ==(lbool const & b) const {
	return ((b.value & 2) & (value & 2)) | (!(b.value & 2) & (value == b.value));
}
inline bool constexpr lbool::operator !=(lbool const & b) const {
	return !(*this == b);
}
inline lbool constexpr lbool::operator ^(bool const & b) const {
	return lbool((uint8_t)(value ^ (uint8_t) b));
}

inline lbool constexpr lbool::operator &&(lbool const & b) const {
//   constexpr uint8_t sel = ((value << 1) | (b.value << 3));
//   constexpr uint8_t v = (0xF7F755F4 >> sel) & 3;
	return lbool(static_cast<uint8_t>((0xF7F755F4 >> ((value << 1) | (b.value << 3))) & 3));
}

inline lbool constexpr lbool::operator ||(lbool const & b) const {
//      constexpr uint8_t sel = (value << 1) | (b.value << 3);
//      constexpr uint8_t v = (0xFCFCF400 >> sel) & 3;
	return lbool(
			static_cast<uint8_t>((0xFCFCF400
					>> static_cast<uint8_t>((value << 1) | (b.value << 3))) & 3));
}
inline constexpr int lbool::toInt() const {
	return value;
}

inline constexpr lbool lbool::toLbool(int const v) {
	return lbool((uint8_t) v);
}

struct Lit {
	int x;

	static constexpr int const undefVal = -1;
	static constexpr int const errorVal = -2;

	constexpr Lit();

	// Use this as a constructor:
	constexpr Lit(Var var, bool sign);

	constexpr Lit operator ~() const;
	constexpr Lit operator ^(bool b) const;

	constexpr bool operator ==(Lit p) const;
	constexpr bool operator !=(Lit p) const;
	constexpr bool operator <(Lit p) const;

	constexpr int toInt() const;

	constexpr bool sign() const;

	constexpr Var var() const;

	constexpr lbool value(lbool const lb) const;

	static constexpr Lit Undef();

	static constexpr Lit Error();

	static constexpr Lit toLit(int i);
private:
	constexpr Lit(int const x) :
			x(x) {
	}
};

inline constexpr Lit Lit::toLit(int i) {
	return Lit(i);
}

inline constexpr Lit::Lit() :
		x(undefVal) {
}

// Use this as a constructor:
inline constexpr Lit::Lit(Var var, bool sign) :
		x(var + var + (int) sign) {
}

inline constexpr Lit Lit::operator ~() const {
	return Lit(x ^ 1);
}
inline constexpr Lit Lit::operator ^(bool b) const {
	return Lit(x ^ (unsigned int) b);
}

inline constexpr bool Lit::operator ==(Lit p) const {
	return x == p.x;
}
inline constexpr bool Lit::operator !=(Lit p) const {
	return x != p.x;
}
inline constexpr bool Lit::operator <(Lit p) const {
	return x < p.x;
}  // '<' makes p, ~p adjacent in the ordering.

inline constexpr lbool Lit::value(lbool const lb) const {
	return lb ^ sign();
}
inline constexpr Var Lit::var() const {
	return x >> 1;
}

inline constexpr Lit Lit::Undef() {
	return Lit(undefVal);
}
inline constexpr Lit Lit::Error() {
	return Lit(errorVal);
}

inline constexpr int Lit::toInt() const {
	return x;
}
inline constexpr bool Lit::sign() const {
	return x & 1;
}

inline std::ostream& operator<<(std::ostream& out, const Lit& val) {
	out << ((val.sign()) ? -(val.var()) : (val.var())) << std::flush;
	return out;
}

// Needed for OccLists
template<typename Type>
inline int toInt(Type const & in) {
	return in;
}
template<>
inline int toInt<Lit>(Lit const & in) {
	return in.toInt();
}

}

#endif /* SOURCES_DATABASE_BASICTYPES_H_ */
