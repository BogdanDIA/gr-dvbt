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
#include <complex>
#include <gr_math.h>
#include <gr_expj.h>
#include <stdio.h>
#include <volk/volk.h>
#include <gr_fxpt.h>

void print_float(float f)
{
  for (int j = 0; j < 4; j++)
    printf("%x", 0xff & ((char *)&f)[j]);
}

void print_double(float d)
{
  for (int j = 0; j < 8; j++)
    printf("%x", 0xff & ((char *)&d)[j]);
}


namespace gr {
  namespace dvbt {

    

    int 
    ofdm_sym_acquisition_impl::peak_detect_init(float threshold_factor_rise, float threshold_factor_fall, int look_ahead, float alpha)
    {
      d_avg_alpha = alpha;
      d_threshold_factor_rise = threshold_factor_rise;
      d_threshold_factor_fall = threshold_factor_fall;
      d_avg = 0;

      return (0);
    }

    int 
    ofdm_sym_acquisition_impl::peak_detect_process(const float * datain, const int datain_length, int * peak_pos)
    {
#if 0
      int ret = 0;
      float peak = datain[0];

      for (int i = 0; i < datain_length; i++)
        if (datain[i] > peak)
        {
          peak = datain[i];
          *peak_pos = i;
        }

      if (*peak_pos == datain_length - 1)
        ret = 0;
      else
        ret = 1;

      return (ret);
#endif
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

    int
    ofdm_sym_acquisition_impl::ml_sync(const gr_complex * in, int lookup_start, int lookup_stop, int * cp_pos, gr_complex * derot, int * to_consume, int * to_out)
    {

      assert(lookup_start >= lookup_stop);
      assert(lookup_stop >= (d_cp_length + d_fft_length - 1));

      // Array to store peak positions
      int peak_pos[d_fft_length];
      float d_phi[d_fft_length];

      // Calculate norm
#ifdef USE_VOLK
      if(is_unaligned())
        volk_32fc_magnitude_squared_32f_u(d_norm, in, d_cp_length);
      else
        volk_32fc_magnitude_squared_32f_a(d_norm, in, d_cp_length);
#else
      for (int i = lookup_start; i >= (lookup_stop - (d_cp_length + d_fft_length - 1)); i--)
        d_norm[i] = std::norm(in[i]);
#endif
      
      // Calculate gamma on each point
#ifdef USE_VOLK
      if (is_unaligned())
        volk_32fc_x2_multiply_conjugate_32fc_u(d_corr, &in[0], &in[d_fft_length], d_cp_length + d_fft_length);
      else
        volk_32fc_x2_multiply_conjugate_32fc_a(d_corr, &in[0], &in[d_fft_length], d_cp_length + d_fft_length);
#else
      for (int i = lookup_start; i >= (lookup_stop - d_cp_length - 1); i--)
        d_corr[i - d_fft_length] = in[i] * std::conj(in[i - d_fft_length]);
#endif

      // Calculate time delay and frequency correction
      // This looks like spagetti code but it is fast
      for (int i = lookup_start - 1; i >= lookup_stop; i--)
      {
        int k = i - lookup_stop;

        d_phi[k] = 0.0;
        d_gamma[k] = 0.0;

        // Moving sum for calculating gamma and phi
        for (int j = 0; j < d_cp_length; j++)
        {
          // Calculate gamma and store it
          d_gamma[k] += d_corr[i - j - d_fft_length];
          // Calculate phi and store it
          d_phi[k] += d_norm[i - j] + d_norm[i - j - d_fft_length];
        }
      }

      // Init lambda with gamma
#ifdef USE_VOLK
      if (is_unaligned())  
        volk_32fc_magnitude_32f_a(d_lambda, d_gamma, d_fft_length);
      else
        volk_32fc_magnitude_32f_u(d_lambda, d_gamma, d_fft_length);
#else
      for (int i = 0; i < (lookup_start - lookup_stop); i++)
        d_lambda[i] = std::abs(d_gamma[i]);
#endif

      // Calculate lambda
#ifdef USE_VOLK
      //
#else
      for (int i = 0; i < (lookup_start - lookup_stop); i++)
        d_lambda[i] -= (d_phi[i] * d_rho / 2.0);
#endif

#if 0
      // Print input
      for (int i = 0; i < (2 * d_fft_length + d_cp_length); i++)
        printf("in[%i]: %.10f\n", i, d_lambda[i]);

#endif
      // Print lambda
      for (int i = 0; i < (lookup_start - lookup_stop); i++)
        printf("lambda[%i]: %.10f\n", i, d_lambda[i]);

      // Find peaks of lambda
      int peak_length = peak_detect_process(&d_lambda[0], (lookup_start - lookup_stop), &peak_pos[0]);
      int peak = 0;
      // We found an end of symbol at peak_pos[0] + CP + FFT
      if (peak_length)
      {
        peak = peak_pos[0] + lookup_stop;
        printf("peak: %i, peak_pos[0]: %i\n",  peak, peak_pos[0]);

        *cp_pos = peak;

        // Calculate frequency correction
        float peak_epsilon = gr_fast_atan2f(d_gamma[peak_pos[0]]);
        double sensitivity = (double)(-1) / (double)d_fft_length;

        //printf("peak_epsilon: %.10f\n", peak_epsilon);
        //printf("d_phaseinc before fft: %.10f\n", d_phaseinc);

        // Store phases for derotating the signal
        // We always process CP len + FFT len
        for (int i = 0; i < (d_cp_length + d_fft_length); i++)
        {
          if (i == d_nextpos)
            d_phaseinc = d_nextphaseinc;

          // We are interested only in fft_length
          d_phase += d_phaseinc;

          while (d_phase > (float)M_PI)
            d_phase -= (float)(2.0 * M_PI);
          while (d_phase < (float)(-M_PI))
            d_phase += (float)(2.0 * M_PI);

          //printf("d_phase[%i]: %.10f, d_phaseinc: %.10f\n", d_index++, d_phase, d_phaseinc);

#ifdef USE_FIXEDARCTAN
          float oq, oi;
          gr_int32 angle = gr_fxpt::float_to_fixed (d_phase);
          gr_fxpt::sincos(angle, &oq, &oi);

          derot[i] = gr_complex(oi, oq);
#else
          derot[i] = gr_expj(d_phase);
#endif
        }
            
        d_nextphaseinc = sensitivity * peak_epsilon;
        d_nextpos = peak - (d_cp_length + d_fft_length);
              
        //printf("d_phaseinc after fft: %.10f\n", d_phaseinc);
        //printf("d_phase after fft: %f\n", d_phase);

        *to_consume = d_cp_length + d_fft_length;
        *to_out = 1;
      }
      else
      {
        printf("miss, phaseinc: %.10f, phase: %.10f, \n", d_phaseinc, d_phase);
#if 0
        // Print input
        for (int i = 0; i < (2 * d_fft_length + d_cp_length); i++)
          printf("in[%i]: %.10f\n", i, d_lambda[i]);

        // Print lambda
        for (int i = 0; i < d_fft_length; i++)
          printf("lambda[%i]: %.10f\n", i, d_lambda[i]);
#endif
        for (int i = 0; i < (d_cp_length + d_fft_length); i++)
        {
          d_phase += d_phaseinc;

          while (d_phase > (float)M_PI)
            d_phase -= (float)(2.0 * M_PI);
          while (d_phase < (float)(-M_PI))
            d_phase += (float)(2.0 * M_PI);
        }

        // We consume only fft_length
        *to_consume = d_cp_length + d_fft_length;
        *to_out = 0;
      }

      return (peak_length);
    }

    int
    ofdm_sym_acquisition_impl::cp_sync(const gr_complex * in, int * cp_pos, gr_complex * derot, int * to_consume, int * to_out)
    {
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
      d_index(0), d_phase(0.0), d_phaseinc(0.0), d_cp_found(0), d_count(0), d_nextphaseinc(0), d_nextpos(0), \
        d_sym_acq_count(0),d_sym_acq_timeout(100), d_initial_aquisition(0), \
        d_freq_correction_count(0), d_freq_correction_timeout(10)
    {
      set_relative_rate(1.0 / (double) (d_cp_length + d_fft_length));

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

      d_derot = new gr_complex[d_cp_length + d_fft_length];
      if (d_derot == NULL)
        std::cout << "cannot allocate memory" << std::endl;

      d_conj = new gr_complex[2 * d_fft_length + d_cp_length];
      if (d_conj == NULL)
        std::cout << "cannot allocate memory" << std::endl;

      d_norm = new float[2 * d_fft_length + d_cp_length];
      if (d_norm == NULL)
        std::cout << "cannot allocate memory" << std::endl;

      d_corr = new gr_complex[2 * d_fft_length + d_cp_length];
      if (d_corr == NULL)
        std::cout << "cannot allocate memory" << std::endl;

      //peak_detect_init(0.2, 0.25, 30, 0.0005);
      peak_detect_init(0.8, 0.9, 30, 0.9);
    }

    /*
     * Our virtual destructor.
     */
    ofdm_sym_acquisition_impl::~ofdm_sym_acquisition_impl()
    {
      delete [] d_gamma;
      delete [] d_lambda;
      delete [] d_derot;
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

        int cp_start = 0;

        // This is initial aquisition of symbol start
        while (!d_initial_aquisition && (++d_sym_acq_count < d_sym_acq_timeout)){
          d_initial_aquisition = ml_sync(in, 2 * d_fft_length + d_cp_length - 1, d_fft_length + d_cp_length - 1, \
              &d_cp_start, &d_derot[0], &d_to_consume, &d_to_out);
        }

        // This is fractional frequency correction (pre FFT)
        if (d_initial_aquisition)
        {
          cp_start = d_cp_start;

          d_cp_found = ml_sync(in, d_cp_start + 4, d_cp_start - 4, \
              &cp_start, &d_derot[0], &d_to_consume, &d_to_out);

          if (d_cp_found)
          {
            trigger[0] = 1;

            d_freq_correction_count = 0;

            // Derotate the signal and output
            int j = 0;
            for (int i = (d_cp_start - d_fft_length + 1); i <= d_cp_start; i++)
            {
              out[j] = d_derot[j] * in[i];
              j++;
            }
          }
          else
          {
            trigger[0] = 0;

            // If we have a number of consecutive misses then we restart aquisition
            if (++d_freq_correction_count > d_freq_correction_timeout)
            {
              d_initial_aquisition = 0;
              d_freq_correction_count = 0;

              printf("restart aquisition\n");
            }
          }
        }

        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each(d_to_consume);

        // Tell runtime system how many output items we produced.
        return (d_to_out);
    }
  } /* namespace dvbt */
} /* namespace gr */

