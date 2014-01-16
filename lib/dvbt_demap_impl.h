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

#ifndef INCLUDED_DVBT_DVBT_DEMAP_IMPL_H
#define INCLUDED_DVBT_DVBT_DEMAP_IMPL_H

#include <dvbt/dvbt_demap.h>
#include <dvbt/dvbt_config.h>

namespace gr {
  namespace dvbt {
    /*!
     * \brief DVBT demapper class.
     * \ingroup dvbt
     * \param nsize length of input stream \n
     * \param constellation constellation used \n
     * \param hierarchy hierarchy used \n
     * \param transmission transmission mode used \n
     * \param gain gian of complex output stream \n
     */
    class dvbt_demap_impl : public dvbt_demap
    {
    private:
      const dvbt_config config;

      int d_nsize;

      //Constellation size
      unsigned char d_constellation_size;
      //Transmission mode
      dvbt_transmission_mode_t d_transmission_mode;
      //Step on each axis of the constellation
      unsigned char d_step;
      //Keep Alpha internally
      unsigned char d_alpha;
      //Gain for the complex values
      float d_gain;

      gr_complex * d_constellation_points;

      void make_constellation_points(int size, int step, int alpha);
      int find_constellation_value(gr_complex val);
      int bin_to_gray(int val);

    public:
      dvbt_demap_impl(int nsize, dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, dvbt_transmission_mode_t transmission, float gain);
      ~dvbt_demap_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      /*!
       * ETSI EN 300 744 Clause 4.3.5. \n
       * Data input format: \n
       * complex(real(float), imag(float)) \n
       *
       * Data output format: \n
       * 000000Y0Y1 - QAM4 \n
       * 0000Y0Y1Y2Y3 - QAM16 \n
       * 00Y0Y1Y2Y3Y4Y5 - QAM64 \n
       */
      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_DVBT_DEMAP_IMPL_H */

