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
      d_threshold_factor_rise = threshold_factor_rise;
      d_threshold_factor_fall = threshold_factor_fall;
      d_avg = 0;

      return (0);
    }

    int ofdm_sym_acquisition_impl::peak_detect_process(const float * datain, const int datain_length, int * peak_pos)
    {
      int state = 0;
      float peak_val = -(float)INFINITY; int peak_index = 0; int peak_pos_length = 0;

      int i = 0;

      while(i < datain_length)
      {
        if (state == 0)
        {
          if (datain[i] > d_avg * d_threshold_factor_rise)
          {
            //printf("state 0, rise: datain: %f, d_avg: %f, d_threshold_factor_rise: %f\n", datain[i], d_avg, d_threshold_factor_rise);
            state = 1;
          }
          else
          {
            d_avg = d_avg_alpha * datain[i] + (1 - d_avg_alpha) * d_avg;
            //printf("state 0: %i, avg: %f\n", i, d_avg);
            i++;
          }
        } 
        else if (state == 1)
        {
          if (datain[i] > peak_val)
          {
            peak_val = datain[i];
            peak_index = i;
            //printf("state 1, index: %i, peak_val: %f, peak_index: %i, avg: %f, dtain: %f\n", i, peak_val, peak_index, d_avg, datain[i]);
            d_avg = d_avg_alpha * datain[i] + (1 - d_avg_alpha) * d_avg;
            i++;
          }
          else if (datain[i] > d_avg * d_threshold_factor_fall)
          {
            //printf("state 1, fall: %i, avg: %f, datain: %f\n", i, d_avg, datain[i]);
            d_avg = (d_avg_alpha) * datain[i] + (1 - d_avg_alpha) * d_avg;
            i++;
          }
          else
          {
            peak_pos[peak_pos_length] = peak_index;
            //printf("state 1, finish: %i, peak_pos[%i]: %i\n", i, peak_pos_length, peak_pos[peak_pos_length]);
            peak_pos_length++;
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
      d_index(0), d_phase(0.0), d_phaseinc(0)
    {
      set_relative_rate(1.0 / (double) (d_cp_length + d_fft_length));
      // TODO - set this to use data for just one output item
      //set_output_multiple(2);

      d_snr = pow(10, d_snr / 10.0);
      d_rho = d_snr / (d_snr + 1.0);

      printf("blocks: %i\n", blocks);
      printf("fft_length: %i\n", fft_length);
      printf("occupied_tones: %i\n", occupied_tones);
      printf("SNR: %f\n", d_snr);

      d_gamma = new gr_complex[d_fft_length];
      if (d_gamma == NULL)
        std::cout << "cannot allocate memory" << std::endl;

      d_lambda = new float[d_fft_length];
      if (d_lambda == NULL)
        std::cout << "cannot allocate memory" << std::endl;

      peak_detect_init(0.2, 0.25, 30, 0.0005);
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

    /*
     * ML Estimation of Time and Frequency Offset in OFDM systems
     * Jan-Jaap van de Beck
     */

    int
    ofdm_sym_acquisition_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];
        unsigned char *trigger = (unsigned char *) output_items[1];

        // Array to store peak positions
        int peak_pos[d_fft_length];
        int to_consume = 0;

        /********************************************************/
        // Calculate time delay and frequency correction
        // This looks like spgetti code but it is fast
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
          d_lambda[i] = std::abs(d_gamma[i]) - (phi * d_rho / 2.0);
        }

        // Find peaks of lambda
        int peak_length = peak_detect_process(&d_lambda[0], d_fft_length, &peak_pos[0]);

        // We found a CP starting at peak_pos[0]
        if (peak_length)
        {
          trigger[0] = 1;

          // Calculate frequency correction
          float peak_epsilon = std::arg(d_gamma[peak_pos[0]]);

          // Increment phase for data before CP
          d_phase += (float)d_phaseinc * (float)peak_pos[0];
          // Calculate new phase increment (it aplies for the cp and fft lengths)
          d_phaseinc = (float)(1) * peak_epsilon / (float)d_fft_length;
          // Increment phase for CP
          d_phase += (float)d_phaseinc * (float)peak_pos[0];
          // Increment the phase for the rest of data
          // We'll use this to update the final phase
          int copy_phase = d_phase + (float)d_phaseinc * (float)(d_fft_length - peak_pos[0]);

          // Derotate the data of size fft_length
          for (int i = 0; i < d_fft_length; i++)
          {
            // We are interested only in fft_length
            d_phase += d_phaseinc;

            while (d_phase > (float)M_PI)
              d_phase -= (float)M_2_PI;
            while (d_phase < (float)(-M_PI))
              d_phase += (float)M_2_PI;

            out[i] = gr_expj(d_phase) * in[peak_pos[0] + d_cp_length + i];
          }

          // We'll consume cp_length+fft_length
          to_consume = d_cp_length + d_fft_length;
          // Even though we copied the fft data we consume only cp_length+fft_length
          // So restore the phase corresponding to this block
          d_phase = copy_phase;
        }
        else
        {
          trigger[0] = 0;
          to_consume = d_fft_length;

          // TODO - get rid of trigger
          for (int i = 0; i < d_fft_length; i++)
          {
            // Set the phase for next symbol
            d_phase += d_phaseinc; 

            while (d_phase > (float)M_PI)
              d_phase -= (float)M_2_PI;
            while (d_phase < (float)(-M_PI))
              d_phase += (float)M_2_PI;

            // Set the output to zero
            out[i] = gr_complex(0.0, 0.0);
          }
        }

        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each(to_consume);

        // Tell runtime system how many output items we produced.
        return (1);
    }

  } /* namespace dvbt */
} /* namespace gr */

