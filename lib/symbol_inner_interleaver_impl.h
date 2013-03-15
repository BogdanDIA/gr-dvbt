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

#ifndef INCLUDED_DVBT_SYMBOL_INNER_INTERLEAVER_IMPL_H
#define INCLUDED_DVBT_SYMBOL_INNER_INTERLEAVER_IMPL_H

#include <dvbt/symbol_inner_interleaver.h>

namespace gr {
  namespace dvbt {

    class symbol_inner_interleaver_impl : public symbol_inner_interleaver
    {
    private:
      const dvbt_config config;

      int d_ninput;
      int d_noutput;

      unsigned int d_h[2048];

      //Keeps the stmbol index
      unsigned int d_symbol_index;

      void generate_H();
      unsigned int H(unsigned int q);
      unsigned int calculate_R(unsigned int i);

    public:
      symbol_inner_interleaver_impl(int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy);
      ~symbol_inner_interleaver_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_SYMBOL_INNER_INTERLEAVER_IMPL_H */

