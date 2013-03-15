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
#include "bit_inner_interleaver_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    const unsigned char bit_inner_interleaver_impl::d_bsize = 126;

    unsigned char
    bit_inner_interleaver_impl::H(unsigned char e, unsigned char w)
    {
      int rez = 0;

      switch (e)
      {
        case 0:
          rez = w; break;
        case 1:
          rez = (w + 63) % 126; break;
        case 2:
          rez = (w + 105) % 126; break;
        case 3:
          rez = (w + 42) % 126; break;
        case 4:
          rez = (w + 21) % 126; break;
        case 5:
          rez = (w + 84) % 126; break;
        default:
          break;
      }

      return rez;
    }

    bit_inner_interleaver::sptr
    bit_inner_interleaver::make(int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy)
    {
      return gnuradio::get_initial_sptr (new bit_inner_interleaver_impl(ninput, noutput, \
        constellation, hierarchy));
    }

    /*
     * The private constructor
     */
    bit_inner_interleaver_impl::bit_inner_interleaver_impl(int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy)
      : gr_block("bit_inner_interleaver",
		      gr_make_io_signature(1, 2, sizeof(unsigned char) * ninput),
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * noutput)),
      config(constellation, hierarchy),
      d_ninput(ninput), d_noutput(noutput)
    {
      d_v = config.d_m;

      d_perm = (unsigned char *)new unsigned char[d_v * d_bsize];
      if (d_perm == NULL)
      {
        std::cout << "Error allocating d_perm" << std::endl;
        return;
      }

      //Init permutation table
      for (int i = 0; i <  d_bsize * d_v; i++)
      {
        if (config.d_hierarchy == gr::dvbt::NH)
          d_perm[i] = ((i % d_v) / (d_v / 2)) + 2 * (i % (d_v / 2));
        else
        {
          d_perm[i] = (i % (d_v - 2)) / ((d_v - 2) / 2) + 2 * (i % ((d_v -2)/2)) + 2;
        }
      }


      if (ninput % d_bsize)
        std::cout << "Input size must be multiple of block size: " \
          << "ninput: " << ninput << "bsize: " << d_bsize << std::endl;

      if (noutput % d_bsize)
        std::cout << "Output size must be multiple of block size: " \
          << "noutput: " << ninput << "bsize: " << d_bsize << std::endl;

      if (ninput != noutput)
        std::cout << "Input size must be ninput=noutput: "  \
          << "ninput: " << ninput << "noutput: " << noutput << std::endl;
    }

    /*
     * Our virtual destructor.
     */
    bit_inner_interleaver_impl::~bit_inner_interleaver_impl()
    {
      delete [] d_perm;
    }

    void
    bit_inner_interleaver_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
        ninput_items_required[1] = noutput_items;
    }

    int
    bit_inner_interleaver_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *inh = (const unsigned char *) input_items[0];
        const unsigned char *inl = (const unsigned char *) input_items[1];
        unsigned char *out = (unsigned char *) output_items[0];


        unsigned char d_b[d_v][d_bsize];
        int bmax = noutput_items * d_ninput / d_bsize;

        /*
         * Clause 4.3.4.1
         *
         * Data Input format Non-Hierarchical:
         * 000000X0X1 - QPSK
         * 0000X0X1X2X3 - QAM16
         * 00X0X1X2X3X4X5 - QAM64
         *
         * Data Input format Hierarchical:
         * 0000000X0 - H-QPSK
         * 0000000X1 - L-QPSK
         * 000000X0X1 - H-QAM16
         * 000000X0X1 - L-QAM16
         * 000000X0X1 - H-QAM64
         * 0000X0X1X2X3 - L-QAM64
         *
         * Data Output format:
         * 000000B0B1 - QPSK
         * 0000B0B1B2B3 - QAM16
         * 00B0B1B2B3B4B5 - QAM64
         * bit interleaver block size is 126
         */

         for (int bcount = 0; bcount < bmax; bcount++)
         {
          for (int i = 0; i < d_bsize; i++)
          {
            if (config.d_hierarchy == gr::dvbt::NH)
            {
              char c = inh[bcount * d_bsize + i];

              for (int k = 0; k < d_v; k++)
                d_b[d_perm[(d_v * i) + k]][((d_v * i) + k) / d_v] = (c >> (d_v - k - 1)) & 0x1;
            }
            else
            {
              char ch = inh[(bcount * d_bsize * (d_v - 2)) + i];
              char cl = inl[(bcount * d_bsize * (d_v - 2)) + i];
              for (int k = 0; k < (d_v / 2); k++)
              {
                //TODO - error - fix this!
                //High priority input
                d_b[(d_v * i + k) % 2][i / 2] = (ch >> ((d_v / 2) - k - 1)) & 0x1;
                //Low priority input
                d_b[d_perm[d_v * i + k]][(d_v * i + k) / (d_v - 2)] = (cl >> ((d_v / 2) - k - 1)) % 0x1;
              }
            }
          }

          //Take one bit from each interleaver
          //and format the output

          for (int w = 0; w < d_bsize; w++)
          {
            unsigned char val = 0;

            for (int e = 0; e < d_v; e++)
              val = (val << 1) | d_b[e][H(e, w)];

            out[(bcount * d_bsize) + w] = val;
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

