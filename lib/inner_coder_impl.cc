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
#include "inner_coder_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    void
    inner_coder_impl::generate_codeword(unsigned char in, unsigned char &x, unsigned char &y)
    {
      //insert input bit
      d_reg |= ((in & 0x1) << 7);
      d_reg  = d_reg >> 1;

      //generate output G1=171(OCT)
      x = ((d_reg >> 6) ^ (d_reg >> 5) ^ (d_reg >> 4) ^ \
                        (d_reg >> 3) ^ d_reg) & 0x1;
      //generate output G2=133(OCT)
      y = ((d_reg >> 6) ^ (d_reg >> 4) ^ (d_reg >> 3) ^ \
                        (d_reg >> 1) ^ d_reg) & 0x1;
    }

    //TODO - do this based on puncturing matrix
    /*
     * Input e.g rate 2/3:
     * 000000x0x1
     * Output e.g. rate 2/3
     * 00000c0c1c2
     */
    unsigned char
    inner_coder_impl::generate_punctured_code(dvbt_code_rate_t coderate, unsigned char in)
    {
      unsigned char x, y;
      unsigned char c = 0;

      switch(coderate)
      {
        //X1Y1
        case gr::dvbt::C1_2:
          generate_codeword(in, x, y);
          c = (x << 1) | y;
          break;
        //X1Y1Y2
        case gr::dvbt::C2_3:
          generate_codeword(in >> 1, x, y);
          c = (x << 1) | y;
          generate_codeword(in, x, y);
          c = (c << 1) | y;
          break;
        //X1Y1Y2X3
        case gr::dvbt::C3_4:
          generate_codeword(in >> 2, x, y);
          c = (x << 1) | y;
          generate_codeword(in >> 1, x, y);
          c = (c << 1) | y;
          generate_codeword(in, x, y);
          c = (c << 1) | x;
          break;
        //X1Y1Y2X3Y4X5
        case gr::dvbt::C5_6:
          generate_codeword(in >> 4, x, y);
          c = (x << 1) | y;
          generate_codeword(in >> 3, x, y);
          c = (c << 1) | y;
          generate_codeword(in >> 2, x, y);
          c = (c << 1) | x;
          generate_codeword(in >> 1, x, y);
          c = (c << 1) | y;
          generate_codeword(in, x, y);
          c = (c << 1) | x;
          break;
        //X1Y1Y2X3Y4X5Y6X7
        case gr::dvbt::C7_8:
          generate_codeword(in >> 6, x, y);
          c = (x << 1) | y;
          generate_codeword(in >> 5, x, y);
          c = (c << 1) | y;
          generate_codeword(in >> 4, x, y);
          c = (c << 1) | x;
          generate_codeword(in >> 3, x, y);
          c = (c << 1) | y;
          generate_codeword(in >> 2, x, y);
          c = (c << 1) | x;
          generate_codeword(in >> 1, x, y);
          c = (c << 1) | y;
          generate_codeword(in, x, y);
          c = (c << 1) | x;
          break;
        default:
          generate_codeword(in, x, y);
          c = (x << 1) | y;
          break;
      }

      return c;
    }

    inner_coder::sptr
    inner_coder::make(int ninput, int noutput, dvbt_constellation_t constellation, \
        dvbt_hierarchy_t hierarchy, dvbt_code_rate_t coderate)
    {
      return gnuradio::get_initial_sptr (new inner_coder_impl(ninput, noutput, constellation, \
        hierarchy, coderate));
    }

    /*
     * The private constructor
     */
    inner_coder_impl::inner_coder_impl(int ninput, int noutput, dvbt_constellation_t constellation, \
        dvbt_hierarchy_t hierarchy, dvbt_code_rate_t coderate)
      : gr_block("inner_coder",
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * ninput),
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * noutput)),
      config(ninput, noutput, constellation, hierarchy, coderate, coderate),
      d_reg(0),
      d_bitcount(0)
    {
      d_k = config.d_k;
      d_n = config.d_n;
      d_m = config.d_m;
    }

    /*
     * Our virtual destructor.
     */
    inner_coder_impl::~inner_coder_impl()
    {
    }

    void
    inner_coder_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    inner_coder_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        //TODO - allocate these dynamically 
        unsigned char in_buff[d_n * d_m];
        unsigned char out_buff[d_n * d_m];

        /*
         * Clause 4.3.3
         * Mother convolutional code with rate 1/2
         * k=1, n=2, K=6
         * Generator polinomial G1=171(OCT), G2=133(OCT)
         * Punctured to obtain rates of 2/3, 3/4, 5/6, 7/8
         *
         * TODO - comment in and out data format
         */

        int in_count = 0;
        int out_count = 0;
        int krem = d_k;

        while (out_count < (noutput_items * config.d_noutput))
        {
          int ncount = 0;
          //keep buffer into 64bits
          long long ll = 0;
          int y = 0;

          while(ncount < (d_n * d_m))
          {
            int x = 0;
            int c = 0;
            int bitrem = 8;
            int byte = in[in_count];

            //Do this till finish the byte
            while(bitrem >= krem)
            {
              x = ((byte >> (bitrem - krem)) & ((1 << d_k) - 1));
              c = generate_punctured_code(config.d_code_rate_HP, y | x);
              ll = (ll << d_n) | c;
              bitrem -= krem;
              krem = d_k;
              ncount += d_n;
            }
            //Put remaining bits 
            if (bitrem > 0)
            {
              krem = d_k - bitrem;
              y = (byte & ((1 << bitrem) - 1)) << krem;
            }
            in_count++;
          }

          int count = ncount / d_m;

          while(count--)
            out[out_count++] = (ll >> (d_m * count)) & ((1 << d_m) - 1);
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

