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
#include "convolutional_deinterleaver_impl.h"
#include <stdio.h>

//#define DEBUG 1

#ifdef DEBUG
#define PRINTF(a...) printf(a)
#else
#define PRINTF(a...)
#endif


namespace gr {
  namespace dvbt {

    const int convolutional_deinterleaver_impl::d_SYNC = 0x47;
    const int convolutional_deinterleaver_impl::d_NSYNC = 0xB8;
    const int convolutional_deinterleaver_impl::d_MUX_PKT = 8;

    convolutional_deinterleaver::sptr
    convolutional_deinterleaver::make(int nsize, int I, int M)
    {
      return gnuradio::get_initial_sptr (new convolutional_deinterleaver_impl(nsize, I, M));
    }

    /*
     * The private constructor
     */
    convolutional_deinterleaver_impl::convolutional_deinterleaver_impl(int blocks, int I, int M)
      : block("convolutional_deinterleaver",
          io_signature::make(1, 1, sizeof (unsigned char)),
          io_signature::make(1, 1, sizeof (unsigned char) * I * blocks)),
      d_blocks(blocks), d_I(I), d_M(M)
    {
      set_relative_rate(1.0 / I * d_blocks);
      set_output_multiple(2);
      //The positions are shift registers (FIFOs)
      //of lenght i*M
      for (int i = (d_I - 1); i >= 0; i--)
        d_shift.push_back(new std::deque<unsigned char>(d_M * i, 0));

      // There are 8 mux packets
      assert(d_blocks / d_m == d_MUX_PKT);
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

    void
    convolutional_deinterleaver_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      int ninputs = ninput_items_required.size ();

      for (int i = 0; i < ninputs; i++)
        ninput_items_required[i] = noutput_items * d_I * d_blocks;
    }


    int
    convolutional_deinterleaver_impl::general_work(int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        int to_out = noutput_items;

        PRINTF("INTERLEAVER: noutput_items: %i\n", noutput_items);

        for (int count = 0, i = 0; i < to_out; i++)
        {
          for (int mux_pkt = 0; mux_pkt < d_MUX_PKT; mux_pkt++)
          {
            // This is actually the interleaver
            for (int k = 0; k < (d_M * d_I); k++)
            {
              d_shift[k % d_I]->push_back(in[count]);
              out[count++] = d_shift[k % d_I]->front();
              d_shift[k % d_I]->pop_front();
            }
          }
        }

        PRINTF("INTERLEAVER: to_out: %i\n", to_out);

        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each(d_I * d_blocks * to_out);

        // Tell runtime system how many output items we produced.
        return (to_out);
    }

  } /* namespace dvbt */
} /* namespace gr */

