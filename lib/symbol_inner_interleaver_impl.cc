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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gr_io_signature.h>
#include "symbol_inner_interleaver_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    void
    symbol_inner_interleaver_impl::generate_H()
    {
      unsigned int Mmax = 2048;
      unsigned int Nmax = 1512;
      unsigned int Nr = 11;
      unsigned int q = 0;
      unsigned int hq = 0;

      for (int i = 0; i < Mmax; i++)
      {
        d_h[q] = ((i % 2) << (Nr - 1)) + calculate_R(i);
        if (d_h[q] < Nmax) q++;
      }
    }

    unsigned int
    symbol_inner_interleaver_impl::H(unsigned int q)
    {
      return d_h[q];
    }

    unsigned int
    symbol_inner_interleaver_impl::calculate_R(unsigned int i)
    {
      const int Nr = 11;
      unsigned int reg = 0;

      unsigned char bit_perm[Nr - 1] = {4, 3, 9, 6, 2, 8, 1, 5, 7, 0};

      if (i == 0)
        reg = 0;
      else if (i == 1)
        reg = 0;
      else if (reg == 2)
        reg = 1;
      else
      {
        reg = 1;
        for (int k = 3; k <= i; k++)
        {
          unsigned char new_bit = (reg ^ (reg >> 3)) & 0x1;
          reg = ((reg >> 1) | (new_bit << 9)) & 0x3ff;
        }
      }

      unsigned int newreg = 0;

      for (int k = 0; k < (Nr - 1); k++)
      {
        unsigned char bit = (reg >> k) & 0x1;
        newreg = newreg | (bit << bit_perm[k]); 
      }

      return newreg;
    }

    symbol_inner_interleaver::sptr
    symbol_inner_interleaver::make(int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy)
    {
      return gnuradio::get_initial_sptr (new symbol_inner_interleaver_impl(ninput, \
		    noutput, constellation, hierarchy));
    }

    /*
     * The private constructor
     */
    symbol_inner_interleaver_impl::symbol_inner_interleaver_impl(int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy)
      : gr_block("symbol_inner_interleaver",
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * ninput),
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * noutput)),
      config(constellation, hierarchy),
      d_symbol_index(0)
    {
      generate_H();
    }

    /*
     * Our virtual destructor.
     */
    symbol_inner_interleaver_impl::~symbol_inner_interleaver_impl()
    {
    }

    void
    symbol_inner_interleaver_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    symbol_inner_interleaver_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        /*
         * Clause 4.3.4.2
         * One block is 12groupsx126datawords=1512datawords
         * Input format: 
         * 000000I0I1 - QPSK
         * 0000I0I1I2I3 - 16QAM
         * 00I0I1I2I3I4I5 - 64QAM
         * Output format:
         * 000000Y0Y1 - QPSK
         * 0000Y0Y1Y2Y3 - 16QAM
         * 00Y0Y1Y2Y3Y4Y5 - 64QAM
         */

        for (int k = 0; k < noutput_items; k++)
        {
          for (int q = 0; q < d_ninput; q++)
          {
            int blocks = k * d_ninput;

            if (d_symbol_index % 2)
              out[blocks + q] = in[blocks + H(q)];
            else
              out[blocks + H(q)] = in[blocks + q];
          }
          d_symbol_index++;
        }

        // Do <+signal processing+>
        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */

