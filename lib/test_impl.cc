/* -*- c++ -*- */
/* 
 * Copyright 2013 <Bogdan Diaconescu, yo3iiu@yo3iiu.ro>.
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

#include <gnuradio/io_signature.h>
#include "test_impl.h"

#include <iostream>
#include <stdio.h>

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
      : block("test",
          io_signature::make(1, 1, itemsize * ninput),
          io_signature::make(1, 1, itemsize)),
      d_ninput(ninput)
    {
      printf("in_bsize: %i\n", input_signature()->sizeof_stream_item (0));
      printf("out_bsize: %i\n", output_signature()->sizeof_stream_item (0));
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

      int ibsize = input_signature()->sizeof_stream_item (0);
      int obsize = output_signature()->sizeof_stream_item (0);

      int size = std::min (ninput_items[0] * ibsize, noutput_items * obsize);

      memcpy(out, in, size); 

      // Do <+signal processing+>
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (size / ibsize);
      printf("ninput_items: %i\n", ninput_items[0]);
      printf("noutput_items: %i\n", noutput_items);

      // Tell runtime system how many output items we produced.
      return (noutput_items);
    }

  } /* namespace dvbt */
} /* namespace gr */

