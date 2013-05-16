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
#include "reed_solomon_enc_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    reed_solomon_enc::sptr
    reed_solomon_enc::make(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks)
    {
      return gnuradio::get_initial_sptr (new reed_solomon_enc_impl(p, m, gfpoly, n, k, t, s, blocks));
    }

    /*
     * The private constructor
     */
    reed_solomon_enc_impl::reed_solomon_enc_impl(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks)
      : gr_block("reed_solomon",
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * blocks * (k - s)),
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * blocks * (n - s))),
      d_p(p), d_m(m), d_gfpoly(gfpoly), d_n(n), d_k(k), d_t(t), d_s(s), d_blocks(blocks),
      d_rs(p, m, gfpoly, n, k, t, s, blocks)
    {

    }

    /*
     * Our virtual destructor.
     */
    reed_solomon_enc_impl::~reed_solomon_enc_impl()
    {
    }

    void
    reed_solomon_enc_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    reed_solomon_enc_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];
 
        int in_bsize = d_k - d_s;
        int out_bsize = d_n - d_s;

        unsigned char parity[2 * d_t];

        int in_count = 0;
        int out_count = 0;

        for (int i = 0; i < (d_blocks * noutput_items); i++)
        {
          //TODO - zero copy?
          //memcpy(&d_in[d_s], &in[i * in_bsize], in_bsize);

          d_rs.rs_encode(&in[i * in_bsize], parity);

          memcpy(&out[i * out_bsize], &in[i * in_bsize], in_bsize);
          memcpy(&out[i * out_bsize + in_bsize], parity, 2 * d_t);
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

