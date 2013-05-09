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
#include "energy_descramble_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    const int energy_descramble_impl::d_nblocks = 8;
    const int energy_descramble_impl::d_bsize = 188;
    const int energy_descramble_impl::d_SYNC = 0x47;
    const int energy_descramble_impl::d_NSYNC = 0xB8;

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
      : gr_block("energy_descramble",
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * d_nblocks * d_bsize),
		      gr_make_io_signature(1, 1, sizeof (unsigned char)))
    {
      set_relative_rate((double) (d_nblocks * d_bsize)); 
      set_output_multiple(d_nblocks * d_bsize);
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
        ninput_items_required[0] = noutput_items / (d_nblocks * d_bsize);
    }

    int
    energy_descramble_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        int index = 0;
        int count = 0;
        int to_consume = 0;
        int to_out = 0;
        int is_sync = 0;

        // Search for NSYNC byte
        while(is_sync == 0 && index < d_bsize)
          if (in[index] == d_NSYNC)
            is_sync = 1;
          else
            index++;

        // Calculate the input items to consume
        // in order to produce noutput_items
        // It is important to set_output_multiple first
        int isize = noutput_items / (d_nblocks * d_bsize);

        // If we found a NSYNC byte
        if (is_sync)
        {
          for (int i = 0; i < isize; i++)
          {
            init_prbs();
	          int sync = d_NSYNC;

            for (int j = 0; j < d_nblocks; j++)
            {
              if (in[index + count] != sync)
                printf("error: Missing sync: %i on input!\n", sync);

              out[count++] = d_SYNC;
	            // PRBS clocking starts right after NSYNC

              for (int k = 1; k < d_bsize; k++)
                out[count++] = in[index + count] ^ clock_prbs(d_nblocks);

              // For subsequent blocks PRBS is clocked also on SYNC
              // but its output is not used
              sync = d_SYNC;
              clock_prbs(d_nblocks);
            }
          }

          to_consume = isize;
          to_out = noutput_items;
        }
        else
        {
          to_consume = 1;
          to_out = 0;
        }

        // Tell runtime how many input items we consumed
        consume_each(to_consume);

        // Tell runtime system how many output items we produced.
        return to_out;
    }

  } /* namespace dvbt */
} /* namespace gr */

