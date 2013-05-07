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
#include <complex>
#include "dvbt_map_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {
 
    dvbt_map::sptr
    dvbt_map::make(int nsize, dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, \
        dvbt_transmission_mode_t transmission, float gain)
    {
      return gnuradio::get_initial_sptr (new dvbt_map_impl(nsize, constellation, hierarchy, transmission, gain));
    }

    /*
     * The private constructor
     */
    dvbt_map_impl::dvbt_map_impl(int nsize, dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, dvbt_transmission_mode_t transmission, float gain)
      : gr_block("dvbt_map",
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * nsize),
		      gr_make_io_signature(1, 1, sizeof (gr_complex) * nsize)),
      config(constellation, hierarchy, gr::dvbt::C1_2, gr::dvbt::C1_2, gr::dvbt::G1_32, transmission),
      d_nsize(nsize),
      d_gain(gain)
    {
      //TODO - clean up here
      unsigned char c_gain[4] = {0x9b, 0xe8, 0xa1, 0xbe};
      float * pf = (float *)&c_gain;
      d_gain = (float)(-1) / float(*pf);

      //Get parameters from config object
      d_transmission_mode = config.d_transmission_mode;
      d_alpha = config.d_alpha;
      d_bits_per_symbol = config.d_m;

      d_qaxis_points = d_bits_per_symbol - 2;
      d_qaxis_steps = d_qaxis_points - 1;
    }

    /*
     * Our virtual destructor.
     */
    dvbt_map_impl::~dvbt_map_impl()
    {
    }

    unsigned int
    dvbt_map_impl::bin_to_gray(unsigned int val)
    {
      return (val >> 1) ^ val;
    }

    void
    dvbt_map_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items;
    }

    int
    dvbt_map_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];
        char q0, q1, x, y;
        float v_x, v_y;
        
        for (int i = 0; i < noutput_items; i++)
        {
          for (int k = 0; k < d_nsize; k++)
          {
            unsigned char bits = in[k + i * d_nsize];

            //TODO - use VOLK for multiplication
            switch (config.d_constellation)
            {
              case (gr::dvbt::QPSK):
              {
                q0 = (bits >> 1) & 0x1; q1 = bits & 0x1;
                v_x = (float)(d_alpha) / d_gain; v_y = (float)(d_alpha) / d_gain;
              }; break;
              case (gr::dvbt::QAM16):
              {
                q0 = (bits >> 3) & 0x1; q1 = (bits >> 2) & 0x1;
                x = (bits >> 1) & 0x1; y = bits & 0x1;

                v_x = (float)(d_alpha + 2 * (d_qaxis_steps - x)) / d_gain;
                v_y = (float)(d_alpha + 2 * (d_qaxis_steps - y)) / d_gain;
              }; break;
              case (gr::dvbt::QAM64):
              {
                q0 = (bits >> 5) & 0x1; q1 = (bits >> 4) & 0x1;

                x = ((bits >> 2) | (bits >> 1)) & 0x3; y = ((bits >> 1) | bits) & 0x3;

                v_x = (float)(d_alpha + 2 * (d_qaxis_steps - bin_to_gray(x))) / d_gain;
                v_y = (float)(d_alpha + 2 * (d_qaxis_steps - bin_to_gray(y))) / d_gain;
              } break;
              default:
              {
                //Defaults to QPSK
                q0 = (bits >> 1) & 0x1; q1 = bits & 0x1;
                v_x = (float)(d_alpha) / d_gain; v_y = (float)(d_alpha) / d_gain;
              }; break;
            }

            int sign0 = 1 - 2 * q0; 
            int sign1 = 1 - 2 * q1; 
            out[k + d_nsize * i] = gr_complex(sign0 * v_x, sign1 * v_y);
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

