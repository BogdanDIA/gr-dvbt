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
          io_signature::make(1, 1, sizeof (unsigned char)))
    {
      set_relative_rate((double) (d_nblocks * d_bsize)); 
      set_output_multiple(2 * d_nblocks * d_bsize);

      PRINTF("ENERGY: d_nblocks: %i, d_bsize: %i\n", d_nblocks, d_bsize);
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

        int to_consume = noutput_items / (d_nblocks * d_bsize);
        int to_out = noutput_items;

        PRINTF("ENERGY: noutput_items: %i\n", noutput_items);

        for (int count = 0, i = 0; i < to_consume; i++)
        {
          init_prbs();

          for (int mux_pkt = 0; mux_pkt < d_MUX_PKT; mux_pkt++)
          {
            out[count++] = d_SYNC;
            // PRBS clocking starts right after NSYNC

            for (int k = 1; k < d_bsize; k++)
              out[count++] = in[count] ^ clock_prbs(d_nblocks);

            // For subsequent blocks PRBS is clocked also on SYNC
            // but its output is not used
            clock_prbs(d_nblocks);
          }
        }

        PRINTF("ENERGY: to_consume: %i, to_out: %i\n", to_consume, to_out);

        // Tell runtime how many input items we consumed
        consume_each(to_consume);

        // Tell runtime system how many output items we produced.
        return (to_out);
    }

  } /* namespace dvbt */
} /* namespace gr */

