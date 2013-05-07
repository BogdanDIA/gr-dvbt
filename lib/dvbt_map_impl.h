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

#ifndef INCLUDED_DVBT_DVBT_MAP_IMPL_H
#define INCLUDED_DVBT_DVBT_MAP_IMPL_H

#include <dvbt/dvbt_map.h>
#include <dvbt/dvbt_config.h>

namespace gr {
  namespace dvbt {
    /*!
     * \brief DVBT mapper class.
     * \ingroup dvbt
     * \param nsize length of input stream \n
     * \param constellation constellation used \n
     * \param gain gian of complex output stream \n
     */

    class dvbt_map_impl : public dvbt_map
    {
    private:
      const dvbt_config config;

      int d_nsize;

      // Keeps transmission mode
      dvbt_transmission_mode_t d_transmission_mode;

      //Keep Alpha internally
      unsigned char d_alpha;
      //Quadrant of the current symbol
      unsigned char d_quadrant;
      //Number of bits per symbol
      unsigned char d_bits_per_symbol;
      //Number of points on x,y axis in a quadrant
      unsigned char d_qaxis_points;
      //Number of steps on x,y axis in a quadrant
      unsigned char d_qaxis_steps;
      //Gain for the complex values
      float d_gain;

      //Return gray representation from natural binary
      unsigned int bin_to_gray(unsigned int val);

    public:
      dvbt_map_impl(int nsize, dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, dvbt_transmission_mode_t transmission, float gain);
      ~dvbt_map_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

        /*!
	 * ETSI EN 300 744 Clause 4.3.5. \n
         * Data input format: \n
         * 000000Y0Y1 - QAM4 \n
         * 0000Y0Y1Y2Y3 - QAM16 \n
         * 00Y0Y1Y2Y3Y4Y5 - QAM64 \n
         *
         * Data output format: \n
         * complex(real(float), imag(float)) \n
         */
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_DVBT_MAP_IMPL_H */

