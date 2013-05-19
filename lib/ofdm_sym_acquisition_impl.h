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

#ifndef INCLUDED_DVBT_OFDM_SYM_ACQUISITION_IMPL_H
#define INCLUDED_DVBT_OFDM_SYM_ACQUISITION_IMPL_H

#include <dvbt/ofdm_sym_acquisition.h>

namespace gr {
  namespace dvbt {

    class ofdm_sym_acquisition_impl : public ofdm_sym_acquisition
    {
    private:
      // Nothing to declare in this block.

    public:
      ofdm_sym_acquisition_impl(int blocks, int fft_length, int occupied_tones, int cp_length, float snr);
      ~ofdm_sym_acquisition_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_OFDM_SYM_ACQUISITION_IMPL_H */

