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

#ifndef INCLUDED_DVBT_INNER_CODER_IMPL_H
#define INCLUDED_DVBT_INNER_CODER_IMPL_H

#include <dvbt/inner_coder.h>
#include <dvbt/dvbt_config.h>

namespace gr {
  namespace dvbt {

    class inner_coder_impl : public inner_coder
    {
    private:
      const dvbt_config config;

      unsigned char d_reg;

      //counts the bits in the bytes
      //in input stream
      unsigned char d_bitcount;

      // Code rate k/n
      unsigned char d_k;
      unsigned char d_n;
      // Cosntellation with m
      unsigned char d_m;



      void generate_codeword(unsigned char in, unsigned char &x, unsigned char &y);
      unsigned char generate_punctured_code(dvbt_code_rate_t coderate, unsigned char c);

    public:
      inner_coder_impl(int ninput, int noutput, dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, dvbt_code_rate_t coderate);
      ~inner_coder_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_INNER_CODER_IMPL_H */

