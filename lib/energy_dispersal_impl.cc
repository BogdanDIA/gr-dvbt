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
#include "energy_dispersal_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    const int energy_dispersal_impl::d_nblocks = 8;
    const int energy_dispersal_impl::d_bsize = 188;
    const int energy_dispersal_impl::d_sync = 0x47;
    const int energy_dispersal_impl::d_nsync = 0xB8;

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
		      gr_make_io_signature(1, 1, d_nblocks * d_bsize * sizeof(unsigned char)),
		      gr_make_io_signature(1, 1, d_nblocks * d_bsize * sizeof(unsigned char)))
    {}

    /*
     * Our virtual destructor.
     */
    energy_dispersal_impl::~energy_dispersal_impl()
    {
    }

    void
    energy_dispersal_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    energy_dispersal_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        int count = 0;

        //TODO - look for sync first

        for (int i = 0; i < noutput_items; i++)
        {
          init_prbs();

          int sync = d_nsync;

          for (int j = 0; j < d_nblocks; j++)
          {
            out[count++] = sync;

            for (int k = 0; k < (d_bsize - 1); k++)
              out[count++] = in[count] ^ clock_prbs(d_nblocks);

            sync = d_sync;
            clock_prbs(d_nblocks);
          }
        }

        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */

