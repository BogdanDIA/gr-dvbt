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

namespace gr {
  namespace dvbt {

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
          gr_make_io_signature2(2, 2, sizeof (gr_complex) * blocks * fft_length, sizeof (char) * blocks * fft_length))
    {}

    /*
     * Our virtual destructor.
     */
    ofdm_sym_acquisition_impl::~ofdm_sym_acquisition_impl()
    {
    }

    void
    ofdm_sym_acquisition_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    ofdm_sym_acquisition_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const float *in = (const float *) input_items[0];
        float *out = (float *) output_items[0];

        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */

