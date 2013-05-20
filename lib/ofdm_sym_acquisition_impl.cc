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
#include "ofdm_sym_acquisition_impl.h"
#include <complex.h>
#include <gr_math.h>
#include <gr_expj.h>
#include <stdio.h>

namespace gr {
  namespace dvbt {

    int ofdm_sym_acquisition_impl::peak_detect_init(float threshold_factor_rise, float threshold_factor_fall, int look_ahead, float alpha)
    {
      d_avg_alpha = alpha;
      d_avg = 0;

      return (0);
    }

    int ofdm_sym_acquisition_impl::peak_detect_process(float * datain, int datain_length, unsigned char * peak_pos, int &peak_pos_length)
    {
      int state = 0;
      float peak_val = -(float)INFINITY; int peak_ind = 0;

      peak_pos_length = 0;
      int i = 0;

      while(i < datain_length)
      {
        if (state == 0)
        {
          if (datain[i] > d_avg * d_threshold_factor_rise)
          {
            state = 1;
          }
          else
          {
            d_avg = d_avg_alpha * datain[i] + (1 - d_avg_alpha) * d_avg;
            i++;
          }
        } 
        else if (state == 1)
        {
          if (datain[i] > peak_val)
          {
            peak_val = datain[i];
            peak_ind = i;
            d_avg = d_avg_alpha * datain[i] + (1 - d_avg_alpha) * d_avg;
            i++;
          }
          else if (datain[i] > d_avg * d_threshold_factor_fall)
          {
            d_avg = (d_avg_alpha) * datain[i] + (1 - d_avg_alpha) * d_avg;
            i++;
          }
          else
          {
            peak_pos[peak_pos_length++] = peak_ind;
            state = 0;
            peak_val = - (float)INFINITY;
          }
        }
      }

      return (peak_pos_length);
    }


    ofdm_sym_acquisition::sptr
    ofdm_sym_acquisition::make(int blocks, int fft_length, int occupied_tones, int cp_length, float snr)
    {
      return gnuradio::get_initial_sptr (new ofdm_sym_acquisition_impl(blocks, fft_length, occupied_tones, cp_length, snr));
    }

    /*
     * The private constructor
     */
    ofdm_sym_acquisition_impl::ofdm_sym_acquisition_impl(int blocks, int fft_length, int occupied_tones, int cp_length, float snr)
      : gr_block("ofdm_sym_acquisition",
          gr_make_io_signature(1, 1, sizeof (gr_complex) * blocks),
          gr_make_io_signature2(2, 2, sizeof (gr_complex) * blocks * fft_length, sizeof (char) * blocks * fft_length)),
      d_blocks(blocks), d_fft_length(fft_length), d_cp_length(cp_length), d_snr(snr),
      d_search_max(2 * d_fft_length + d_cp_length)
    {
      set_relative_rate(1.0 / (double) (d_cp_length + d_fft_length));
      // TODO - set this to use data for just one output item
      //set_output_multiple(2);

      d_snr = pow(10, d_snr / 10.0);
      d_rho = d_snr / (d_snr + 1.0);

      printf("blocks: %i\n", blocks);
      printf("fft_length: %i\n", fft_length);
      printf("occupied_tones: %i\n", occupied_tones);
      printf("snr: %f\n", snr);

      d_gamma = new gr_complex[d_search_max];
      if (d_gamma == NULL)
        std::cout << "cannot allocate memory" << std::endl;

      d_lambda = new float[d_search_max];
      if (d_lambda == NULL)
        std::cout << "cannot allocate memory" << std::endl;

      peak_detect_init(0.2, 0.25, 30, 0.0005);

      d_index = 0;
    }

    /*
     * Our virtual destructor.
     */
    ofdm_sym_acquisition_impl::~ofdm_sym_acquisition_impl()
    {
      delete [] d_gamma;
      delete [] d_lambda;
    }

    void
    ofdm_sym_acquisition_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      int ninputs = ninput_items_required.size ();

      // make sure we receive at least (symbol_length + fft_length)
      for (int i = 0; i < ninputs; i++)
        ninput_items_required[i] = (2 * d_fft_length + d_cp_length) * noutput_items;
    }

    int
    ofdm_sym_acquisition_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];
        unsigned char *trigger = (unsigned char *) output_items[1];

        printf("ninput_items: %i\n", ninput_items[0]);
        printf("noutput_items: %i\n", noutput_items);

        // Array to store peak positions
        unsigned char peak_pos[d_search_max];
        int to_consume = 0;

        {
          for (int i = 0; i < (d_fft_length); i++)
          {
            float phi = 0.0;

            d_gamma[i] = 0;

            // Moving sum for calculating gamma and phi
            for (int j = 0; j < d_cp_length; j++)
            {
              gr_complex c1 = in[i + j];
              gr_complex c2 = in[i + j + d_fft_length];

              // Calculate gamma
              // Keep gamma for later
              d_gamma[i] += c1 * std::conj(c2);
              // Calculate phi
              phi += std::norm(c1) + std::norm(c2);
            }

            // Keep angle for later
            // Calculate arg(gamma - rho*phi)
            d_lambda[i] = abs(d_gamma[i]) - (phi * d_rho / 2.0);

            //printf("gamma[%i]: %f, %f\n", i, d_gamma[i].real(), d_gamma[i].imag());
            //printf("in[%i]: %f, %f\n", i, in[i].real(), in[i].imag());
            printf("lambda[%i]: %f\n", i, d_lambda[i]);
          }

          int peak_length = 0;

          // Find peaks of lambda
          peak_detect_process(&d_lambda[0], d_fft_length, &peak_pos[0], peak_length);

          if (peak_length)
          {
            for (int i = 0; i < peak_length; i++)
              printf("peak[%i]: %i\n", i, i);

            // Calculate frequency correction
            float peak_epsilon = (- 1.0 / M_2_PI) * std::arg(d_gamma[peak_pos[0]]);
            printf("angle: %f\n", peak_epsilon);

            trigger[0] = 1;

            // We found a CP starting at peak_pos[0]
            // Copy symbol data after CP
            for (int i = 0; i < d_fft_length; i++)
              out[i] = /*gr_expj(peak_epsilon * 2 * M_2_PI * i / (float)d_fft_length) */ in[i + d_cp_length + peak_pos[0]];

            // Then we consumed fft+cp length 
            //to_consume += d_fft_length + d_cp_length;
          }
          else
          {
            printf("miss\n");
            //to_consume += d_fft_length;
          }

          int addr1 = (d_index) * 8 * 2112;
          printf("work: index: %i, addr1: %x, %x, %x, %x, %x\n", d_index, addr1, \
              ((char *)in)[0], ((char *)in)[1], ((char *)in)[2], ((char *)in)[3]);

          int addr2 = (d_index) * 8 * 2112 + 8 * 2048;
          printf("work: index: %i, addr2: %x, %x, %x, %x, %x\n", d_index, addr2, \
              ((char *)in)[8 * 2048], ((char *)in)[1 + 8 * 2048], \
              ((char *)in)[2 + 8 * 2048], ((char *)in)[3 + 8 * 2048]);

          d_index++;
        }

        to_consume =  d_fft_length + d_cp_length;

        printf("to_consume: %i\n", to_consume);

        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each(to_consume);

        // Tell runtime system how many output items we produced.
        return (1); //noutput_items
    }

  } /* namespace dvbt */
} /* namespace gr */

