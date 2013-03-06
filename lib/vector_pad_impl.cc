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
#include "vector_pad_impl.h"

namespace gr {
  namespace dvbt {

    vector_pad::sptr
    vector_pad::make(int itemsize, int ninput, int noutput)
    {
      return gnuradio::get_initial_sptr (new vector_pad_impl(itemsize, ninput, noutput));
    }

    /*
     * The private constructor
     */
    vector_pad_impl::vector_pad_impl(int itemsize, int ninput, int noutput)
      : gr_block("vector_pad",
		      gr_make_io_signature(1, 1, itemsize * ninput),
		      gr_make_io_signature(1, 1, itemsize * noutput))
    {
      d_itemsize = itemsize;
      d_ninput = ninput;
      d_noutput = noutput;

      d_prefix_len = d_itemsize * (1 + (int)((d_noutput - d_ninput) / 2));
      d_ninput_len = d_itemsize * d_ninput;
      d_suffix_len = d_itemsize * d_noutput - d_ninput_len - d_prefix_len;
    }

    /*
     * Our virtual destructor.
     */
    vector_pad_impl::~vector_pad_impl()
    {
    }

    void
    vector_pad_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
      ninput_items_required[0] = noutput_items;
    }

    int
    vector_pad_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const char *in = (const char *) input_items[0];
      char *out = (char *) output_items[0];

      int icnt = 0;
      int ocnt = 0;
      int blocks = 0;
      
      
      for (int i = 0; i < noutput_items; i++)
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

