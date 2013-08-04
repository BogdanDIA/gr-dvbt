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

#include <gr_io_signature.h>
#include "convolutional_deinterleaver_impl.h"
#include <stdio.h>

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
      : gr_sync_decimator("convolutional_deinterleaver",
		      gr_make_io_signature(1, 1, sizeof (unsigned char)),
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * I * blocks), I * blocks),
      d_blocks(blocks), d_I(I), d_M(M), d_index(0), d_offset(0)
    {
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

    int
    convolutional_deinterleaver_impl::work(int noutput_items,
			  gr_vector_const_void_star &input_items,
			  gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        int count = 0;
        int to_out = 0;

        PRINTF("INTERLEAVER: d_offset: %i, d_index: %i, noutput_items: %i\n", d_offset, d_index, noutput_items);

        // Output only (noutput_items - 1) items in order to 
        // have room for an initial offset
        to_out = noutput_items - 1;

        for (int i = 0; i < to_out; i++)
        {
          int sync = d_NSYNC;
          int mux_pkt = 0;

          while (mux_pkt < d_MUX_PKT)
          {
            if (in[d_index + count] != sync)
            {
              PRINTF("INTERLEAVER: error: Missing sync: %x on input[%i] %x!\n", sync, d_index + count, in[d_index + count]);

              // Go to the next input element
              d_index++;
              // Reset output counter
              count = 0;
              // Restart the search for NSYNC on the next index
              sync = d_NSYNC;
              // Reset MUX packet to 0
              mux_pkt = 0;

              if (d_index == d_blocks * d_I)
              {
                d_index = 0;
                break;
              }

              continue;
            }
            else
              PRINTF("INTERLEAVER: sync found: %x on input[%i] %x, mux_pkt: %i\n", sync, d_index + count, in[d_index + count], mux_pkt);

            // This is actually the actual interleaver
            for (int k = 0; k < (d_M * d_I); k++)
            {
              d_shift[k % d_I]->push_back(in[d_index + count]);
              out[count++] = d_shift[k % d_I]->front();
              d_shift[k % d_I]->pop_front();
            }

            sync = d_SYNC;

            mux_pkt++;
          }
        }

        d_offset += to_out * d_blocks * d_I;

        PRINTF("INTERLEAVER: to_out: %i\n", to_out);

        // Tell runtime system how many output items we produced.
        return (to_out);
    }

  } /* namespace dvbt */
} /* namespace gr */

