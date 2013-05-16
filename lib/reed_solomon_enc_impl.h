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

#ifndef INCLUDED_DVBT_REED_SOLOMON_ENC_IMPL_H
#define INCLUDED_DVBT_REED_SOLOMON_ENC_IMPL_H

#include <dvbt/reed_solomon_enc.h>

namespace gr {
  namespace dvbt {
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
    class reed_solomon_impl : public reed_solomon
    {
    private:
      int d_lambda;
      int d_gfpoly;
      int d_p;
      int d_m;
      int d_n;
      int d_k;
      int d_t;
      int d_s;
      int d_blocks;

    public:
      reed_solomon_impl(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks);
      ~reed_solomon_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

	/*!
         * ETSI EN 300 744 Clause 4.3.2 \n
         * RS(N=204,K=239,T=8)
         */
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_REED_SOLOMON_ENC_IMPL_H */

