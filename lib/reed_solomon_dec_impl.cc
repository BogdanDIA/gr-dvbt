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
#include "reed_solomon_dec_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    reed_solomon_dec::sptr
    reed_solomon_dec::make(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks)
    {
      return gnuradio::get_initial_sptr (new reed_solomon_dec_impl(p, m, gfpoly, n, k, t, s, blocks));
    }

    /*
     * The private constructor
     */
    reed_solomon_dec_impl::reed_solomon_dec_impl(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks)
      : gr_block("reed_solomon_dec",
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * blocks * (n - s)),
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * blocks * (k - s))),
      d_p(p), d_m(m), d_gfpoly(gfpoly), d_n(n), d_k(k), d_t(t), d_s(s), d_blocks(blocks),
      d_rs(p, m, gfpoly, n, k, t, s, blocks)
    {
    
    }

    /*
     * Our virtual destructor.
     */
    reed_solomon_dec_impl::~reed_solomon_dec_impl()
    {
    }

    void
    reed_solomon_dec_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    reed_solomon_dec_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        // We receive only nonzero data
        int in_bsize = d_n - d_s;
        int out_bsize = d_k - d_s;

        unsigned char d_in[d_n];

        int in_count = 0;
        int out_count = 0;

        for (int i = 0; i < (d_blocks * noutput_items); i++)
        {
          //TODO - zero copy?
          // Set first d_s symbols to zero
          memset(&d_in[0], 0, d_s);
          // Then copy actual data
          memcpy(&d_in[d_s], &in[i * in_bsize], in_bsize);

#if 0
          for (int j = 0; j < 8; j++)
          {
            unsigned char c = (unsigned char)rand();
#define min(a, b) (a) < (b) ? (a) : (b)
            c = min(c, d_n - d_s);
            d_in[d_s + c] = 0;
            printf("ccc: %i\n", c);
          }
#endif

#if 0
          int start = 134;

          for (int j = (start - 1); j >= (start - 1 - 7); j--)
            d_in[j] = 0;
#endif

          d_rs.rs_decode(d_in, NULL, 0);

          memcpy(&out[i * out_bsize], &d_in[d_s], out_bsize);
        }

        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */

