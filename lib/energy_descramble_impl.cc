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
#include "energy_descramble_impl.h"
#include <stdio.h>

//#define DEBUG 1

#ifdef DEBUG
#define PRINTF(a...) printf(a)
#else
#define PRINTF(a...)
#endif

namespace gr {
  namespace dvbt {

    const int energy_descramble_impl::d_nblocks = 8;
    const int energy_descramble_impl::d_bsize = 188;
    const int energy_descramble_impl::d_SYNC = 0x47;
    const int energy_descramble_impl::d_NSYNC = 0xB8;
    const int energy_descramble_impl::d_MUX_PKT = 8;

    void
    energy_descramble_impl::init_prbs()
    {
      d_reg = 0xa9;
    }

    int
    energy_descramble_impl::clock_prbs(int clocks)
    {
      int res = 0;
      int feedback = 0;

      for(int i = 0; i < clocks; i++)
      {
        feedback = ((d_reg >> (14 - 1)) ^ (d_reg >> (15 - 1))) & 0x1;
        d_reg = ((d_reg << 1) | feedback) & 0x7fff;

        res = (res << 1) | feedback;
      }

      return res;
    }

    energy_descramble::sptr
    energy_descramble::make(int nblocks)
    {
      return gnuradio::get_initial_sptr (new energy_descramble_impl(nblocks));
    }

    /*
     * The private constructor
     */
    energy_descramble_impl::energy_descramble_impl(int nblocks)
      : block("energy_descramble",
          io_signature::make(1, 1, sizeof (unsigned char) * d_nblocks * d_bsize),
          io_signature::make(1, 1, sizeof (unsigned char))),
      d_index (0)
    {
      set_relative_rate((double) (d_nblocks * d_bsize)); 
      set_output_multiple(2 * d_nblocks * d_bsize);

      PRINTF("d_nblocks: %i, d_bsize: %i\n", d_nblocks, d_bsize);
    }

    /*
     * Our virtual destructor.
     */
    energy_descramble_impl::~energy_descramble_impl()
    {
    }

    void
    energy_descramble_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = 10 + (noutput_items / (d_nblocks * d_bsize));
    }

    int
    energy_descramble_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        int count = 0;
        int to_consume = 0;
        int to_out = 0;

        PRINTF("ENERGY: d_offset: %i, d_index: %i, noutput_items: %i\n", d_offset, d_index, noutput_items);

        // Calculate the input items to consume
        // in order to produce noutput_items
        // It is important to set_output_multiple() first
        int isize = (noutput_items / (d_nblocks * d_bsize)) - 1;

        // Output only (noutput_items - 1) items in order to 
        // have room for an initial offset
        to_out = isize * d_nblocks * d_bsize;
        to_consume = isize;

        for (int i = 0; i < isize; i++)
        {
          int sync = d_NSYNC;
          int mux_pkt = 0;

          init_prbs();

          while (mux_pkt < d_MUX_PKT)
          {
            if (in[d_index + count] != sync)
            {
              PRINTF("ENERGY: error: Missing sync: %x on input[%i] %x!\n", sync, d_index + count, in[d_index + count]);

              // Go to the next input element
              d_index++;
              // Reset output counter
              count = 0;
              // Restart the search for NSYNC on the next index
              sync = d_NSYNC;
              // Reset MUX packet to 0
              mux_pkt = 0;

              if (d_index == d_nblocks * d_bsize)
              {
                d_index = 0;
                break;
              }

              continue;
            }
            else
              PRINTF("ENERGY: sync found: %x on input[%i] %x, mux_pkt: %i\n", sync, d_index + count, in[d_index + count], mux_pkt);

            out[count++] = d_SYNC;
            // PRBS clocking starts right after NSYNC

            for (int k = 1; k < d_bsize; k++)
              out[count++] = in[d_index + count] ^ clock_prbs(d_nblocks);

            // For subsequent blocks PRBS is clocked also on SYNC
            // but its output is not used
            sync = d_SYNC;
            clock_prbs(d_nblocks);

            mux_pkt++;
          }
        }

        // For debug purposes
        d_offset += to_consume * d_nblocks * d_bsize;

        PRINTF("ENERGY: to_consume: %i, to_out: %i\n", to_consume, to_out);

        // Tell runtime how many input items we consumed
        consume_each(to_consume);

        // Tell runtime system how many output items we produced.
        return (to_out);
    }

  } /* namespace dvbt */
} /* namespace gr */

