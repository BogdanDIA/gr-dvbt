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
#include "test_impl.h"

#include <iostream>

namespace gr {
  namespace dvbt {

    test::sptr
    test::make(int itemsize, int ninput, int noutput)
    {
      return gnuradio::get_initial_sptr (new test_impl(itemsize, ninput, noutput));
    }

    /*
     * The private constructor
     */
    test_impl::test_impl(int itemsize, int ninput, int noutput)
      : gr_block("test",
		      gr_make_io_signature(1, 1, itemsize * ninput),
		      gr_make_io_signature(1, 1, itemsize * noutput))
    {
      d_itemsize = itemsize;
      d_ninput_len = itemsize * ninput;
      d_noutput_len = itemsize * noutput;

      d_prefix_len = itemsize * (1 + (int)((noutput - ninput) / 2));
      d_suffix_len = d_noutput_len - d_ninput_len - d_prefix_len;
    }

    /*
     * Our virtual destructor.
     */
    test_impl::~test_impl()
    {
    }

    void
    test_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
      ninput_items_required[0] = noutput_items;
    }

    int
    test_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const char *in = (const char *) input_items[0];
        char *out = (char *) output_items[0];

        int icnt = 0;
        int ocnt = 0;
        int blocks = 0;
        
        
        while (blocks < noutput_items)
        {
          // add prefix
          memset(&out[ocnt], 0, d_prefix_len);
          ocnt += d_prefix_len;

          // add useful data
          memcpy(&out[ocnt], &in[icnt], d_ninput_len);
          icnt += d_ninput_len;
          ocnt += d_ninput_len;
          
          // add suffix
          memset(&out[ocnt], 0, d_suffix_len);
          ocnt += d_suffix_len;

          blocks++;
        }

        // Do <+signal processing+>
        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return (noutput_items);
    }

  } /* namespace dvbt */
} /* namespace gr */

