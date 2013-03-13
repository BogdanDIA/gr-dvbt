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
#include "convolutional_interleaver_impl.h"
#include <deque>

namespace gr {
  namespace dvbt {

    convolutional_interleaver::sptr
    convolutional_interleaver::make(int nsize, int I, int M)
    {
      return gnuradio::get_initial_sptr (new convolutional_interleaver_impl(nsize, I, M));
    }

    /*
     * The private constructor
     */
    convolutional_interleaver_impl::convolutional_interleaver_impl(int blocks, int I, int M)
      : gr_block("convolutional_interleaver",
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * I * blocks),
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * I * blocks)),
      d_blocks(blocks), d_I(I), d_M(M)
    {
      //Zero elements in the first position
      d_shift.push_back(new std::deque<unsigned char>(0, 0));

      //The rest of positions are shift registers (FIFOs)
      //of lenght i*M
      for (int i = 1; i < d_I; i++)
      {
        d_shift.push_back(new std::deque<unsigned char>(d_M * i, 0));
      }
    }

    /*
     * Our virtual destructor.
     */
    convolutional_interleaver_impl::~convolutional_interleaver_impl()
    {
      for (int i = 1; i < d_shift.size(); i++)
      {
        delete d_shift.back();
        d_shift.pop_back();
      }
    }

    void
    convolutional_interleaver_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    convolutional_interleaver_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        /*
         * Clause 4.3.1
         * Forney (Ramsey type III) convolutional interleaver
         * Data input: Blocks of 204 bytes
         * Data output: Blocks of 204 bytes
         */
        for (int i = 0; i < (d_blocks * noutput_items); i++)
        {
          //Process one block
          for (int j = 0; j < d_shift.size(); j++)
          {
            d_shift[j]->push_front(in[(d_I * i) + j]);
            out[(d_I * i) + j] = d_shift[j]->back();
            d_shift[j]->pop_back();
          }
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

