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

    const char symbol_inner_interleaver_impl::d_bit_perm_2k[] = {4, 3, 9, 6, 2, 8, 1, 5, 7, 0};
    const char symbol_inner_interleaver_impl::d_bit_perm_8k[] = {7, 1, 4, 2, 9, 6, 8, 10, 0, 3, 11, 5};

    void
    symbol_inner_interleaver_impl::generate_H()
    {
      const int Mmax = d_fft_length;
      const int Nmax = d_payload_length;
      const int Nr = int(ceil(log2(d_fft_length)));
      int q = 0;
      int hq = 0;

      for (int i = 0; i < Mmax; i++)
      {
        d_h[q] = ((i % 2) << (Nr - 1)) + calculate_R(i);
        if (d_h[q] < Nmax) q++;
      }
    }

    int
    symbol_inner_interleaver_impl::H(int q)
    {
      return d_h[q];
    }

    int
    symbol_inner_interleaver_impl::calculate_R(int i)
    {
      const int Nr = int(ceil(log2(d_fft_length)));
      int reg = 0;

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
          char new_bit = 0;

          if (d_transmission_mode == gr::dvbt::T2k)
            new_bit = (reg ^ (reg >> 3)) & 1;
          else if (d_transmission_mode == gr::dvbt::T8k)
            new_bit = (reg ^ (reg >> 1) ^ (reg >> 4) ^ (reg >> 6)) & 1;
          else
            new_bit = (reg ^ (reg >> 3)) & 1;

          int mask = (1 << Nr) - 1;
          reg = ((reg >> 1) | (new_bit << (Nr - 2))) & mask;
        }
      }

      int newreg = 0;

      for (int k = 0; k < (Nr - 1); k++)
      {
        char bit = (reg >> k) & 1;
        newreg = newreg | (bit << d_bit_perm[k]); 
      }

      return newreg;
    }

    symbol_inner_interleaver::sptr
    symbol_inner_interleaver::make(int nsize, \
        dvbt_transmission_mode_t transmission, int direction)
    {
      return gnuradio::get_initial_sptr (new symbol_inner_interleaver_impl(nsize, \
		    transmission, direction));
    }

    /*
     * The private constructor
     */
    symbol_inner_interleaver_impl::symbol_inner_interleaver_impl(int nsize, \
        dvbt_transmission_mode_t transmission, int direction)
      : gr_block("symbol_inner_interleaver",
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * nsize),
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * nsize)),
      config(gr::dvbt::QAM16, gr::dvbt::NH, gr::dvbt::C1_2, gr::dvbt::C1_2, gr::dvbt::G1_32, transmission),
      d_nsize(nsize), d_direction(direction),
      d_fft_length(0), d_payload_length(0),
      d_symbol_index(0)
    {
      d_symbols_per_frame = config.d_symbols_per_frame;
      d_transmission_mode = config.d_transmission_mode;
      d_fft_length = config.d_fft_length;
      d_payload_length = config.d_payload_length;
      d_direction = direction;

      // Verify if transmission mode matches with size of block
      assert(d_payload_length != d_nsize);

      // Allocate memory for h vector
      d_h = new int[d_fft_length];
      if (d_h == NULL)
      {
        std::cout << "cannot allocate memory" << std::endl;
      }

      // Setup bit permutation vectors
      if (d_transmission_mode == gr::dvbt::T2k)
        d_bit_perm = d_bit_perm_2k;
      else if (d_transmission_mode == gr::dvbt::T8k)
        d_bit_perm = d_bit_perm_8k;
      else
        d_bit_perm = d_bit_perm_2k;

      // Generate the h function
      generate_H();
    }

    /*
     * Our virtual destructor.
     */
    symbol_inner_interleaver_impl::~symbol_inner_interleaver_impl()
    {
      delete [] d_h;
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

        for (int k = 0; k < noutput_items; k++)
        {
          for (int q = 0; q < d_nsize; q++)
          {
            int blocks = k * d_nsize;

            if (d_symbol_index % 2)
              if (d_direction)
                out[blocks + q] = in[blocks + H(q)];
              else
                out[blocks + H(q)] = in[blocks + q];
            else
              if (d_direction)
                out[blocks + H(q)] = in[blocks + q];
              else
                out[blocks + q] = in[blocks + H(q)];
          }
          d_symbol_index = (++d_symbol_index) % d_symbols_per_frame;
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

