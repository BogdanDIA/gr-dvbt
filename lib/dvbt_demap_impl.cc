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
#include <volk/volk.h>

#include <gr_io_signature.h>
#include <dvbt/dvbt_config.h>
#include "dvbt_demap_impl.h"
#include <gr_math.h>
#include <stdio.h>
#include <sys/time.h>

struct timeval tvs, tve;
struct timezone tzs, tze;

#define USE_VOLK 1

namespace gr {
  namespace dvbt {

    dvbt_demap::sptr
    dvbt_demap::make(int nsize, dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, \
        dvbt_transmission_mode_t transmission, float gain)
    {
      return gnuradio::get_initial_sptr (new dvbt_demap_impl(nsize, constellation, hierarchy, transmission, gain));
    }

    /*
     * The private constructor
     */
    dvbt_demap_impl::dvbt_demap_impl(int nsize, dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, \
        dvbt_transmission_mode_t transmission, float gain)
      : gr_block("dvbt_demap",
		      gr_make_io_signature(1, 1, sizeof (gr_complex) * nsize),
		      gr_make_io_signature(1, 1, sizeof (unsigned char) * nsize)),
      config(constellation, hierarchy, gr::dvbt::C1_2, gr::dvbt::C1_2, gr::dvbt::G1_32, transmission),
      d_nsize(nsize),
      d_constellation_size(0),
      d_step(0),
      d_alpha(0),
      d_gain(0.0)
    {
      //Get parameters from config object
      d_constellation_size = config.d_constellation_size;
      d_transmission_mode = config.d_transmission_mode;
      d_step = config.d_step;
      d_alpha = config.d_alpha;
      d_gain = gain / config.d_norm;

      printf("d_constellation_size: %i\n", d_constellation_size);
      printf("d_step: %i\n", d_step);
      printf("d_alpha: %i\n", d_alpha);
      printf("d_gain: %f\n", d_gain);

      d_constellation_points = new gr_complex[d_constellation_size];
      if (d_constellation_points == NULL)
      {
        std::cout << "cannot allocate memory" << std::endl;
      }

      d_constellation_bits = new int[d_constellation_size];
      if (d_constellation_bits == NULL)
      {
        std::cout << "cannot allocate memory" << std::endl;
      }

      make_constellation_points(d_constellation_size, d_step, d_alpha);
    }

    /*
     * Our virtual destructor.
     */
    dvbt_demap_impl::~dvbt_demap_impl()
    {
      delete [] d_constellation_points;
      delete [] d_constellation_bits;
    }

    void
    dvbt_demap_impl::make_constellation_points(int size, int step, int alpha)
    {
      //TODO - verify if QPSK works
      
      int bits_per_axis = log2(size) / 2;

      for (int i = 0; i < size; i++)
      {
        int x = i >> bits_per_axis;
        int y = i & ((1 << bits_per_axis) - 1);

        d_constellation_points[i] = gr_complex(alpha + step - 2 * x, alpha + step - 2 * y);
        d_constellation_bits[i] = (bin_to_gray(x) << bits_per_axis) + bin_to_gray(y);

        // ETSI EN 300 744 Clause 4.3.5
        // Actually the constellation is gray coded
        // but the bits on each axis are not taken in consecutive order
        // So we need to convert from b0b2b4b1b3b5->b0b1b2b3b4b5(QAM64)

        x = 0; y = 0;

        for (int j = 0; j < bits_per_axis; j++)
        {
          x += ((d_constellation_bits[i] >> (1 + 2 * j)) & 1) << j;
          y += ((d_constellation_bits[i] >> (2 * j)) & 1) << j;
        }

        d_constellation_bits[i] = (x << bits_per_axis) + y;

        printf("points: %f, %f, bits: %x\n", d_constellation_points[i].real(), d_constellation_points[i].imag(), \
            d_constellation_bits[i]);
      }
    }

    int
    dvbt_demap_impl::find_constellation_value(gr_complex val)
    {
      float min_dist = std::norm(val - d_constellation_points[0]);
      int min_index = 0;


#ifdef USE_VOLK
      float ff[d_constellation_size];

      volk_32fc_x2_square_dist_32f_a(&ff[0], &val, &d_constellation_points[0], d_constellation_size * 8);

      for (int i = 0; i < d_constellation_size; i++)
      {
        //printf("ff: %f\n", ff[i]);
        if (ff[i] < min_dist)
        {
          min_dist = ff[i];
          min_index = i;
        }
      }
#else
      for (int i = 0; i < d_constellation_size; i++)
      {
        float dist = std::norm(val - d_constellation_points[i]);
        //printf("dist: %f\n", dist);

        if (dist < min_dist)
        {
          min_dist = dist;
          min_index = i;
        }
      }
#endif

      return d_constellation_bits[min_index];
    }

    int
    dvbt_demap_impl::bin_to_gray(int val)
    {
      return (val >> 1) ^ val;
    }

    void
    dvbt_demap_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    dvbt_demap_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];

        // TODO - use DFE (Decission Feedback Equalizer)

        //gettimeofday(&tvs, &tzs);

#ifdef USE_VOLk
        // TODO - allocate this on heap
        gr_complex ff[noutput_items * d_nsize];
        gr_complex gg(d_gain, 0.0);

        volk_32fc_s32fc_multiply_32fc_a(&ff[0], &in[0], gg, noutput_items * d_nsize);


        for (int i = 0; i < (noutput_items * d_nsize); i++)
          out[i] = find_constellation_value(ff[i]);
#else
        for (int i = 0; i < (noutput_items * d_nsize); i++)
          out[i] = find_constellation_value(in[i] * d_gain);
#endif
        //gettimeofday(&tve, &tze);
        //printf("us: %f\n", (float) (tve.tv_usec - tvs.tv_usec) / (float) (noutput_items * d_nsize));

        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */

