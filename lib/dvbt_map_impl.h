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

    class dvbt_map_impl : public dvbt_map
    {
    private:
      const dvbt_config config;

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

      //Return natural binary from grey representation
      unsigned int grey_to_bin(unsigned int);

    public:
      dvbt_map_impl(int nsize, dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, float gain);
      ~dvbt_map_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_DVBT_MAP_IMPL_H */

