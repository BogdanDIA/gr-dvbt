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

#ifndef INCLUDED_DVBT_REED_SOLOMON_IMPL_H
#define INCLUDED_DVBT_REED_SOLOMON_IMPL_H

#include <dvbt/reed_solomon.h>

namespace gr {
  namespace dvbt {

    class reed_solomon_impl : public reed_solomon
    {
    private:
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

      int gf_add(int a, int b);
      int gf_mul(int a, int b);
      int gf_exp(int a);
      void gf_init(int p, int m, int gfpoly);
      void gf_uninit();
      void rs_init(int lambda, int n, int k, int t);
      void rs_uninit();
      int rs_encode(const unsigned char *data_in, unsigned char *data_out);

    public:
      reed_solomon_impl(int p, int m, int gfpoly, int n, int k, int t);
      ~reed_solomon_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_REED_SOLOMON_IMPL_H */

