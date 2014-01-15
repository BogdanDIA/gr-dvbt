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

#include <gnuradio/io_signature.h>
#include <dvbt/dvbt_config.h>
#include "dvbt_demap_impl.h"
#include <gnuradio/math.h>
#include <stdio.h>
#include <sys/time.h>

// For timing debug
static struct timeval tvs, tve;
static struct timezone tzs, tze;

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
      : block("dvbt_demap",
          io_signature::make(1, 1, sizeof (gr_complex) * nsize),
          io_signature::make(1, 1, sizeof (unsigned char) * nsize)),
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

      printf("DVBT demap, d_constellation_size: %i\n", d_constellation_size);
      printf("DVBT demap, d_step: %i\n", d_step);
      printf("DVBT demap, d_alpha: %i\n", d_alpha);
      printf("DVBT demap, d_gain: %f\n", d_gain);

      d_constellation_points = new gr_complex[d_constellation_size];
      if (d_constellation_points == NULL)
      {
        std::cout << "cannot allocate d_constellation_points" << std::endl;
      }

      d_constellation_bits = new int[d_constellation_size];
      if (d_constellation_bits == NULL)
      {
        std::cout << "cannot allocate d_constellation_bits" << std::endl;
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
      
      // The simetry of the constellation is used to calculate 
      // 16-QAM from QPSK and 64-QAM form 16-QAM

      int bits_per_axis = log2(size) / 2;
      int steps_per_axis = sqrt(size) / 2 - 1;

      for (int i = 0; i < size; i++)
      {
        // This is the quadrant made of the first two bits starting from MSB
        int q = i >> (2 * (bits_per_axis - 1)) & 3;
        // Sign for correctly calculate I and Q in each quadrant
        int sign0 = (q >> 1) ? -1 : 1; int sign1 = (q & 1) ? -1 : 1;

        int x = (i >> (bits_per_axis - 1)) & ((1 << (bits_per_axis - 1)) - 1);
        int y = i & ((1 << (bits_per_axis - 1)) - 1);

        int xval = alpha + (steps_per_axis - x) * step;
        int yval = alpha + (steps_per_axis - y) * step;

        d_constellation_points[i] = gr_complex(sign0 * xval, sign1 * yval);
        d_constellation_bits[i] = (bin_to_gray(x) << (bits_per_axis - 1)) + bin_to_gray(y);

        // ETSI EN 300 744 Clause 4.3.5
        // Actually the constellation is gray coded
        // but the bits on each axis are not taken in consecutive order
        // So we need to convert from b0b2b4b1b3b5->b0b1b2b3b4b5(QAM64)

        x = 0; y = 0;

        for (int j = 0; j < (bits_per_axis - 1); j++)
        {
          x += ((d_constellation_bits[i] >> (1 + 2 * j)) & 1) << j;
          y += ((d_constellation_bits[i] >> (2 * j)) & 1) << j;
        }

        d_constellation_bits[i] = (q << 2 * (bits_per_axis - 1)) + (x << (bits_per_axis - 1)) + y;

        printf("DVBT map, constellation points[%i]: %f, %f, bits: %x\n", \
            i, d_constellation_points[i].real(), d_constellation_points[i].imag(), d_constellation_bits[i]);
      }
    }

    int
    dvbt_demap_impl::find_constellation_value(gr_complex val)
    {
      float min_dist = std::norm(val - d_constellation_points[0]);
      int min_index = 0;

#ifdef USE_VOLK
      float ff[d_constellation_size];

      if (is_unaligned())
        volk_32fc_x2_square_dist_32f_u(&ff[0], &val, &d_constellation_points[0], d_constellation_size);
      else
        volk_32fc_x2_square_dist_32f_a(&ff[0], &val, &d_constellation_points[0], d_constellation_size);

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

#ifdef USE_VOLK
        // TODO - allocate this on heap
        gr_complex ff[noutput_items * d_nsize];
        gr_complex gg(d_gain, 0.0);

        if (is_unaligned())
          volk_32fc_s32fc_multiply_32fc_u(&ff[0], &in[0], gg, noutput_items * d_nsize);
        else
          volk_32fc_s32fc_multiply_32fc_a(&ff[0], &in[0], gg, noutput_items * d_nsize);

        for (int i = 0; i < (noutput_items * d_nsize); i++)
          out[i] = find_constellation_value(ff[i]);
#else
        for (int i = 0; i < (noutput_items * d_nsize); i++)
          out[i] = find_constellation_value(in[i] * d_gain);
#endif
        //gettimeofday(&tve, &tze);
        //printf("dvbt demap: us: %f\n", (float) (tve.tv_usec - tvs.tv_usec) / (float) (noutput_items * d_nsize));

        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */

