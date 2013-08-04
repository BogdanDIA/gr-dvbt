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
#include "energy_dispersal_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    const int energy_dispersal_impl::d_nblocks = 8;
    const int energy_dispersal_impl::d_bsize = 188;
    const int energy_dispersal_impl::d_SYNC = 0x47;
    const int energy_dispersal_impl::d_NSYNC = 0xB8;

    void
    energy_dispersal_impl::init_prbs()
    {
      d_reg = 0xa9;
    }

    int
    energy_dispersal_impl::clock_prbs(int clocks)
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

    energy_dispersal::sptr
    energy_dispersal::make(int nblocks)
    {
      return gnuradio::get_initial_sptr (new energy_dispersal_impl(nblocks));
    }

    /*
     * The private constructor
     */
    energy_dispersal_impl::energy_dispersal_impl(int nblocks)
      : gr_block("energy_dispersal",
		      gr_make_io_signature(1, 1, sizeof(unsigned char)),
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * d_nblocks * d_bsize))
    {
      set_relative_rate(1.0/(double) (d_nblocks * d_bsize)); 
    }

    /*
     * Our virtual destructor.
     */
    energy_dispersal_impl::~energy_dispersal_impl()
    {
    }

    void
    energy_dispersal_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      // Add one block size for SYNC search
      ninput_items_required[0] = d_nblocks * (d_bsize + 1) * noutput_items;
    }

    int
    energy_dispersal_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        int index = 0;
        int count = 0;
        int ret = 0;
        int is_sync = 0;

        // Search for SYNC byte
        while(is_sync == 0 && index < d_bsize)
          if (in[index] == d_SYNC)
            is_sync = 1;
          else
            index++;

        // If we found a SYNC byte
        if (is_sync)
        {
          for (int i = 0; i < noutput_items; i++)
          {
            init_prbs();

            int sync = d_NSYNC;

            for (int j = 0; j < d_nblocks; j++)
            {
              if (in[index + count] != d_SYNC)
                printf("error: Malformed MPEG-TS!\n");

              out[count++] = sync;

              for (int k = 1; k < d_bsize; k++)
                out[count++] = in[index + count] ^ clock_prbs(d_nblocks);

              sync = d_SYNC;
              clock_prbs(d_nblocks);
            }
          }
          consume_each(index + d_nblocks * d_bsize * noutput_items);
          ret = noutput_items;
        }
        else
        {
          consume_each(index);
          ret = 0;
        }

        // Tell runtime system how many output items we produced.
        return ret;
    }
  } /* namespace dvbt */
} /* namespace gr */

