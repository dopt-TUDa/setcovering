// -*- C++ -*-
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is modified from polymake                                    */
/*                                                                           */
/*       polymake                                                            */
/*       Copyright (c) 1997-2005 Ewgenij Gawrilow, Michael Joswig            */
/*                               (Technische Universitaet Berlin, Germany)   */
/*       www.math.tu-berlin.de/polymake                                      */
/*       mailto: polymake@math.tu-berlin.de                                  */
/*       polymake is distributed under the GNU General Public License,       */
/*       see file COPYING in the polymake distribution.                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   Bitset.h
 * @brief  Bitset class based on GMP
 * @author Ewgenij Gawrilow, Marc Pfetsch
 */

#ifndef BITSET_H
#define BITSET_H

#include <gmp.h>
#include <ostream>
#include <iterator>
#include <cstddef>

// from polymake: "inlines_int.h"inline
inline int log2_round(unsigned long x, int round)
{
   if (x <= 1UL) 
      return 0;

   int step = (sizeof(unsigned long) * 8)/2;
   int log2 = step;
   unsigned long pow2 = 1UL << step;

   while (x != pow2) 
   {
      if (x < pow2) 
      {
         if (!(step >>= 1)) 
	    return log2-1 + round;
         pow2 >>= step;
         log2 -= step;
      } 
      else 
      {
         if (!(step >>= 1)) 
	    return log2 + round;
         pow2 <<= step;
         log2 += step;
      }
   }
   return log2;
}

inline int log2_floor(unsigned long x) { return log2_round(x,0); }



class Bitset;

class Bitset_iterator
{
   friend class Bitset;
public:
   typedef std::forward_iterator_tag iterator_category;
   typedef int value_type;
   typedef int reference;
   typedef const int* pointer;
   typedef ptrdiff_t difference_type;
   typedef Bitset_iterator iterator;
   typedef Bitset_iterator const_iterator;
protected:
   mpz_srcptr bits;
   int cur;

   static const int bits_per_limb
#if __GNU_MP_VERSION >= 4
#if (__GNU_MP_VERSION >= 5) || (__GNU_MP_VERSION_MINOR > 0)
     =GMP_NUMB_BITS
#else
     =__GMP_BITS_PER_MP_LIMB
#endif
#endif
   ;

   Bitset_iterator(mpz_srcptr bits_arg, int cur_arg) : bits(bits_arg), cur(cur_arg) { }

public:
   Bitset_iterator() : bits(0), cur(0) { }

   explicit Bitset_iterator(mpz_srcptr bits_arg) : bits(bits_arg) { rewind(); }

   reference operator* () const { return cur; }
   pointer operator-> () const { return &cur; }

   iterator& operator++ ()
   {
      ++cur;
      if (!at_end()) cur=mpz_scan1(bits,cur);
      return *this;
   }

   const iterator operator++ (int) { iterator copy=*this; (void)operator++(); return copy; }

   bool operator== (const iterator& it) const { return cur==it.cur; }
   bool operator!= (const iterator& it) const { return !operator==(it); }

   bool at_end() const
   {
      int n = cur/bits_per_limb;
      int diff = n+1 - mpz_size(bits);
      return diff>0 ||
	 (diff==0 && !(mpz_getlimbn(bits, n) & (mp_limb_t)(-1) << (cur%bits_per_limb)));
   }

   void rewind()
   {
      cur=0;
      if (mpz_size(bits)) cur=mpz_scan1(bits,cur);
   }
};


/** Dense set of integers
    Data representation is a bit string as implemented in the GNU Multiprecision Library (mpz_t).
*/
class Bitset
{
protected:
   mpz_t rep;

   void fill1s(size_t n);
public:
   typedef int value_type;
   typedef const int const_reference;
   typedef const_reference reference;
   typedef Bitset_iterator iterator;
   typedef iterator const_iterator;
   typedef Bitset persistent_type;

   /// An empty set, without preallocated storage.
   Bitset() { mpz_init(rep); }

   /** An empty set with preallocated storage for elements 0..n-1.
       It can dynamically grow beyond this limit if needed; to avoid performance penalties, however,
       you should specify here the highest element you really expect to occur. */
   explicit Bitset(int n, bool full=false)
   {
      mpz_init(rep);
      reserve(n);
      if (full && n>0) fill1s(n);
   }

   Bitset(const Bitset& s) { mpz_init_set(rep, s.rep); }

   ~Bitset() { mpz_clear(rep); }


protected:
   template <typename Iterator>
   void _assign(Iterator src)
   {
      for (; !src.at_end(); ++src)
	 mpz_setbit(rep, *src);
   }

   template <typename Iterator>
   void _assign(Iterator src, Iterator src_end)
   {
      for (; src!=src_end; ++src)
	 mpz_setbit(rep, *src);
   }

public:
   Bitset(void (*f)(mpz_ptr,unsigned long), mpz_srcptr a, unsigned long k)
   {
      mpz_init_set(rep,a);
      f(rep,k);
   }

   template <typename Arg>
   Bitset(void (*f)(mpz_ptr,Arg), Arg a)
   {
      mpz_init(rep);
      f(rep,a);
   }

   template <typename Arg1, typename Arg2>
   Bitset(void (*f)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b)
   {
      mpz_init(rep);
      f(rep,a,b);
   }

   template <typename Iterator>
   Bitset(Iterator src, Iterator src_end)
   {
      mpz_init(rep);
      _assign(src,src_end);
   }

   template <typename Iterator>
   Bitset(Iterator src)
   {
      mpz_init(rep);
      _assign(src);
   }

   Bitset& operator= (const Bitset& s)
   {
      mpz_set(rep, s.rep);
      return *this;
   }

   /// Reserve storage for $n$ elements
   void reserve(int n) 
   {
      if (n > rep[0]._mp_alloc*iterator::bits_per_limb)
#if __GNU_MP_VERSION < 4
	 _mpz_realloc(rep, (n+iterator::bits_per_limb-1)/iterator::bits_per_limb);
#else
         mpz_realloc2(rep,n);
#endif
   }

   /// Make the set empty.
   void clear()
   {
      mpz_set_ui(rep, 0);
   }

   void swap(Bitset& s) { mpz_swap(rep, s.rep); }

   friend void relocate(Bitset* from, Bitset* to)
   {
      to->rep[0] = from->rep[0];
   }

   bool empty() const
   {
      return !mpz_sgn(rep);
   }

   /** Count the elements.
       CAUTION: depending on the hardware, can take O(n) time! */
   int size() const
   {
      return mpz_popcount(rep);
   }

   bool contains(int i) const
   {
      return mpz_tstbit(rep, i);
   }

   int front() const
   {
      return mpz_scan1(rep,0);
   }

   int back() const
   {
      int n = mpz_size(rep)-1;
      return n*iterator::bits_per_limb + log2_floor(mpz_getlimbn(rep,n));
   }

   /** Insert an element.
       This is the quickest way to manipulate single elements. */
   Bitset& operator+= (int i)
   {
      mpz_setbit(rep, i);
      return *this;
   }

   /** Remove an element if it existed.
       This is the quickest way to manipulate single elements. */
   Bitset& operator-= (int i)
   {
      mpz_clrbit(rep, i);
      return *this;
   }

   Bitset& operator+= (const Bitset& s)
   {
      mpz_ior(rep, rep, s.rep);
      return *this;
   }

private:
   static void difference(mpz_ptr dst, mpz_srcptr src1, mpz_srcptr src2);
public:
   /// difference
   Bitset& operator-= (const Bitset& s)
   {
      difference(rep, rep, s.rep);
      return *this;
   }

   /// intersection
   Bitset& operator*= (const Bitset& s)
   {
      mpz_and(rep, rep, s.rep);
      return *this;
   }

   Bitset& operator^= (int i)
   {
      if (mpz_tstbit(rep, i))
	 mpz_clrbit(rep, i);
      else
	 mpz_setbit(rep, i);
      return *this;
   }

   Bitset& operator^= (const Bitset& s)
   {
      mpz_xor(rep, rep, s.rep);
      return *this;
   }

   iterator begin() const { return iterator(rep); }
   iterator end() const { return iterator(rep, empty() ? 0 : back()+1); }

   iterator insert(int i)
   {
      *this += i;
      return iterator(rep, i);
   }

   void push_back(int i) { *this += i; }
   void push_front(int i) { *this += i; }
   void erase(int i) { *this -= i; }

   iterator insert(const iterator&, int i) { return insert(i); }

   void erase(const iterator& where) { *this -= *where; }

   void pop_front() { *this -= front(); }

   void pop_back() { *this -= back(); }

   friend
   Bitset operator+ (const Bitset& s1, const Bitset& s2)
   {
      return Bitset(mpz_ior, s1.rep, s2.rep);
   }

   /// alias for operator+
   friend
   Bitset operator| (const Bitset& s1, const Bitset& s2)
   {
      return Bitset(mpz_ior, s1.rep, s2.rep);
   }

   friend
   Bitset operator+ (const Bitset& s1, int i2)
   {
      return Bitset(mpz_setbit, s1.rep, (unsigned long)i2);
   }

   friend
   Bitset operator+ (int i1, const Bitset& s2)
   {
      return Bitset(mpz_setbit, s2.rep, (unsigned long)i1);
   }

   friend
   Bitset operator- (const Bitset& s1, const Bitset& s2)
   {
      return Bitset(difference, s1.rep, s2.rep);
   }

   friend
   Bitset operator- (const Bitset& s1, int i2)
   {
      return Bitset(mpz_clrbit, s1.rep, (unsigned long)i2);
   }

   friend
   Bitset operator* (const Bitset& s1, const Bitset& s2)
   {
      return Bitset(mpz_and, s1.rep, s2.rep);
   }

   /// alias for operator*
   friend
   Bitset operator& (const Bitset& s1, const Bitset& s2)
   {
      return Bitset(mpz_and, s1.rep, s2.rep);
   }

   friend
   Bitset operator* (const Bitset& s1, int i2)
   {
      return s1.contains(i2) ? Bitset(mpz_setbit, (unsigned long)i2) : Bitset();
   }

   friend
   Bitset operator* (int i1, const Bitset& s2)
   {
      return s2.contains(i1) ? Bitset(mpz_setbit, (unsigned long)i1) : Bitset();
   }

   friend
   Bitset operator^ (const Bitset& s1, const Bitset& s2)
   {
      return Bitset(mpz_xor, s1.rep, s2.rep);
   }

   friend
   Bitset operator^ (const Bitset& s1, int i2)
   {
      return Bitset(s1.contains(i2) ? &mpz_clrbit : &mpz_setbit, s1.rep, (unsigned long)i2);
   }

   friend
   Bitset operator^ (int i1, const Bitset& s2)
   {
      return Bitset(s2.contains(i1) ? &mpz_clrbit : &mpz_setbit, s2.rep, (unsigned long)i1);
   }

   friend
   bool operator== (const Bitset& s1, const Bitset& s2)
   {
      return !mpz_cmp(s1.rep, s2.rep);
   }

   friend
   bool operator!= (const Bitset& s1, const Bitset& s2)
   {
      return mpz_cmp(s1.rep, s2.rep);
   }

   friend
   bool operator< (const Bitset& s1, const Bitset& s2)
   {
      return mpz_cmp(s1.rep, s2.rep) < 0;
   }

   friend
   std::ostream& operator<< (std::ostream& os, const Bitset& B)
   {
      os << "{ ";
      const_iterator it, iend = B.end();
      for (it = B.begin(); it != iend; ++it)
	 os << *it << " ";
      os << "}";
      return os;
   }

   mpz_srcptr get_rep() const { return rep; }
};


// Test relative inclusion
int incl(const Bitset& s1, const Bitset& s2);


namespace std {
   inline void swap(Bitset& s1, Bitset& s2) { s1.swap(s2); }
}



// define hash functions
namespace std_ext {
   template <typename _Key> struct hash;

   // hash function for GMP integers
   template <> struct hash<MP_INT> {
   protected:
      size_t _do(mpz_srcptr a) const {
	 size_t result=0;
	 for (int i=0, n=mpz_size(a); i<n; ++i)
	    (result <<= 1) ^= mpz_getlimbn(a, i);
	 return result;
      }
   public:
      size_t operator() (const MP_INT& a) const {
	 return _do(&a);
      }
   };

   // hash function for Bitsets
   template <> struct hash<Bitset> : hash<MP_INT> {
      size_t operator() (const Bitset& s) const {
         return _do(s.get_rep());
      }
   };
}


#endif // BITSET_H
