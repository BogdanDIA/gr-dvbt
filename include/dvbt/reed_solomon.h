/* -*- c++ -*- */
/* 
 * Copyright 2013 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#ifndef INCLUDED_DVBT_REED_SOLOMON_H
#define INCLUDED_DVBT_REED_SOLOMON_H

#include <dvbt/api.h>

namespace gr {
  namespace dvbt {

    /*!
     * \brief <+description+>
     *
     */
    class reed_solomon
    {
    /*!
     * \brief Reed Solomon encoder.
     * \ingroup dvbt
     * \param p characteristic of GF(p^m) \n
     * \param m we use GF(p^m) \n
     * \param gfpoly Generator Polynomial \n
     * \param n length of codeword of RS coder \n
     * \param k length of information sequence of RS encoder \n
     * \param t number of corrected errors \n
     * \param s shortened length \n
     * \param blocks number of blocks to process at once\n
     */

    private:
      int d_lambda;
      int d_gfpoly;
      int d_p;
      int d_m;
      int d_n;
      int d_k;
      int d_t;
      unsigned char *d_gf_exp;
      unsigned char *d_gf_log;
      unsigned char *d_l;
      unsigned char *d_g;

      int d_s;
      int d_blocks;
      unsigned char *d_syn;

      int gf_add(int a, int b);
      int gf_mul(int a, int b);
      int gf_div(int a, int b);
      int gf_exp(int a);
      int gf_pow(int a, int power);
      int gf_lpow(int power);

      void gf_init(int p, int m, int gfpoly);
      void gf_uninit();
      void rs_init(int lambda, int n, int k, int t);
      void rs_uninit();

    public:
      /*!
       * ETSI EN 300 744 Clause 4.3.2 \n
       * RS(N=204,K=239,T=8)
       */
      int rs_encode(const unsigned char *data, unsigned char *parity);
      int rs_decode(unsigned char *data, unsigned char *eras, const int no_eras);

      reed_solomon(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks);
      ~reed_solomon();
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_REED_SOLOMON_H */

