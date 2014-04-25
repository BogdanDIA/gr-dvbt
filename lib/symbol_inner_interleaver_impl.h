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

#ifndef INCLUDED_DVBT_SYMBOL_INNER_INTERLEAVER_IMPL_H
#define INCLUDED_DVBT_SYMBOL_INNER_INTERLEAVER_IMPL_H

#include <dvbt/symbol_inner_interleaver.h>
#include <dvbt/reference_signals.h>

namespace gr {
  namespace dvbt {
    /*!
     * \brief Symbol interleaver.
     * \ingroup dvbt
     * \param ninput length of input stream
     * \param noutput length of output stream
     * \param constellation constellation used
     * \param hierarchy hierarchy used
     */
    class symbol_inner_interleaver_impl : public symbol_inner_interleaver
    {
    private:
      const dvbt_config config;

      int d_symbols_per_frame;
      dvbt_transmission_mode_t d_transmission_mode;
      int d_nsize;
      int d_direction;
      int d_fft_length;
      int d_payload_length;

      int * d_h;
      const char * d_bit_perm;
      static const char d_bit_perm_2k[];
      static const char d_bit_perm_8k[];

      //Keeps the symbol index
      unsigned int d_symbol_index;

      void generate_H();
      int H(int q);
      int calculate_R(int i);

    public:
      symbol_inner_interleaver_impl(int nsize, \
        dvbt_transmission_mode_t transmission, int direction);
      ~symbol_inner_interleaver_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

        /*!
         * ETSI EN 300 744 Clause 4.3.4.2 \n
         * One block is 12groupsx126datawords=1512datawords \n
         *
         * Input format: \n
         * 000000I0I1 - QPSK \n
         * 0000I0I1I2I3 - 16QAM \n
         * 00I0I1I2I3I4I5 - 64QAM \n
         * 
         * Output format: \n
         * 000000Y0Y1 - QPSK \n
         * 0000Y0Y1Y2Y3 - 16QAM \n
         * 00Y0Y1Y2Y3Y4Y5 - 64QAM \n
         */
      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_SYMBOL_INNER_INTERLEAVER_IMPL_H */

