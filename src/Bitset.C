/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is modified from a file from polymake                        */
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

/**@file   Bitset.C
 * @brief  Implementation part of Bitset class based on GMP
 * @author Ewgenij Gawrilow, Marc Pfetsch
 *
 * The main operation is to replace the standard memory calls used by GMP
 * to STL allocators.
 */

#include "Bitset.h"


#if __GNU_MP_VERSION < 4
const int Bitset_iterator::bits_per_limb = mp_bits_per_limb;
#endif

void Bitset::difference(mpz_ptr dst, mpz_srcptr src1, mpz_srcptr src2)
{
   const mp_limb_t *s2=src2->_mp_d;
   if (dst==src1) {
      // assignment
      mp_limb_t *d=dst->_mp_d;

      if (dst->_mp_size <= src2->_mp_size) {
	 // lhs can get shortened
	 mp_limb_t *end=d+dst->_mp_size, *non0=d;
	 while (d<end)
	    if ((*d++ &= ~(*s2++))) non0=d;
	 dst->_mp_size=non0-dst->_mp_d;
      } else {
	 const mp_limb_t *end=s2+src2->_mp_size;
	 while (s2<end)
	    (*d++) &= ~(*s2++);
      }
   } else {
      (void) _mpz_realloc(dst, src1->_mp_size);
      mp_limb_t *d=dst->_mp_d;
      const mp_limb_t *s1=src1->_mp_d;

      if (src1->_mp_size <= src2->_mp_size) {
	 const mp_limb_t *end=s1+src1->_mp_size;
	 mp_limb_t *non0=d;
	 while (s1<end)
	    if ((*d++ = (*s1++) & ~(*s2++))) non0=d;
	 dst->_mp_size=non0-dst->_mp_d;
      } else {
	 const mp_limb_t *end2=s2+src2->_mp_size, *end=d+(dst->_mp_size=src1->_mp_size);
	 while (s2<end2)
	    *d++ = (*s1++) & ~(*s2++);
	 while (d<end)
	    *d++ = *s1++;
      }
   }
}

void Bitset::fill1s(size_t n)
{
   mpz_ptr dst=rep;
   const size_t n_limbs=dst->_mp_size=(n+iterator::bits_per_limb-1)/iterator::bits_per_limb;
   mp_limb_t *l=dst->_mp_d, *l_top=l+n_limbs-1;
   mp_limb_t ones=mp_limb_t(mp_limb_signed_t(-1));	// take care of sign expansion
#if GMP_NAIL_BITS
   ones &= GMP_NUMB_MASK;
#endif
   for (; l<l_top; ++l)
      *l=ones;
   *l=ones>>(n_limbs*iterator::bits_per_limb-n);
}

enum cmp_value { cmp_lt=-1, cmp_eq=0, cmp_gt=1 }; /*line !e830*/

cmp_value sign(int x)
{
   const unsigned int bits=sizeof(int)*8-1;
   return cmp_value((x >> bits) + (((unsigned int)x + ((unsigned int)(-1) >> 1)) >> bits));
}
cmp_value sign(long x)
{
   const unsigned int bits=sizeof(long)*8-1;
   return cmp_value((x >> bits) + (((unsigned long)x + ((unsigned long)(-1L) >> 1)) >> bits));
}

int incl(const Bitset& s1, const Bitset& s2)
{
   mpz_srcptr rep1=s1.get_rep(), rep2=s2.get_rep();
   int size1=mpz_size(rep1), size2=mpz_size(rep2), result=(int) sign(size1-size2);
   for (const mp_limb_t *e1=rep1->_mp_d, *e2=rep2->_mp_d, *stop=e1+std::min(size1,size2);
	e1 != stop;
	++e1, ++e2) {
      const mp_limb_t intersect=*e1 & *e2;
      if (*e1 != intersect) {
	 if (result<0) return 2;
	 result=1;
      }
      if (*e2 != intersect) {
	 if (result>0) return 2;
	 result=-1;
      }
   }
   return result;
}
