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
#include "convolutional_deinterleaver_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    convolutional_deinterleaver::sptr
    convolutional_deinterleaver::make(int nsize, int I, int M)
    {
      return gnuradio::get_initial_sptr (new convolutional_deinterleaver_impl(nsize, I, M));
    }

    /*
     * The private constructor
     */
    convolutional_deinterleaver_impl::convolutional_deinterleaver_impl(int blocks, int I, int M)
      : gr_sync_decimator("convolutional_deinterleaver",
		      gr_make_io_signature(1, 1, sizeof (unsigned char)),
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * I * blocks), I * blocks),
      d_blocks(blocks), d_I(I), d_M(M)
    {
      //The positions are shift registers (FIFOs)
      //of lenght i*M, except the first position
      for (int i = (d_I - 1); i >= 0; i--)
        d_shift.push_back(new std::deque<unsigned char>(d_M * i, 0));
    }

    /*
     * Our virtual destructor.
     */
    convolutional_deinterleaver_impl::~convolutional_deinterleaver_impl()
    {
      for (int i = 0; i < d_shift.size(); i++)
      {
        delete d_shift.back();
        d_shift.pop_back();
      }
    }

    int
    convolutional_deinterleaver_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        int is_sync = 0;

        for (int i = 0; i < (noutput_items * d_blocks); i++)
        {
          //Process one block of I symbols
       	  int j = 0;

          while (j < d_shift.size())
          {
            int c = in[(d_I * i) + j];

            // Look for a SYNC or a /SYNC
            // If we find a SYNC or a /SYNC
            // route it to the branch with the null delay

#if 0
            if (c == 0x47 || c == 0xb8)
            {
              is_sync = 1;
              j = (d_I - 1);
            }

            if (is_sync)
#endif
            {
              d_shift[j]->push_back(c);
              out[(d_I * i) + j] = d_shift[j]->front();
              d_shift[j]->pop_front();
            }

            j++;
          }
        }

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */

