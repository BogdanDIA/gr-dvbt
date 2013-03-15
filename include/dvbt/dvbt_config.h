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


#ifndef INCLUDED_DVBT_DVBT_CONFIG_H
#define INCLUDED_DVBT_DVBT_CONFIG_H

#include <dvbt/api.h>

namespace gr {
  namespace dvbt {

    enum dvbt_constellation_t {
      QPSK = 0,
      QAM16,
      QAM64,
      RES
    };

    enum dvbt_hierarchy_t {
      NH = 0,
      ALPHA1,
      ALPHA2,
      ALPHA4,
      HRES1,
      HRES2,
      HRES3,
      HRES4
    };

    enum dvbt_code_rate_t {
      C1_2 = 0,
      C2_3,
      C3_4,
      C5_6,
      C7_8,
      CRES1,
      CRES2,
      CRES3
    };

    enum dvbt_transmission_mode_t {
      T2k = 0,
      T8k = 1,
      TRES1 = 2,
      TRES2 = 3,
    };

    enum dvbt_guard_interval_t {
      G1_32 = 0,
      G1_16,
      G1_8,
      G1_4
    };

    /*!
     * \brief <+description+>
     *
     */
    class DVBT_API dvbt_config
    {

      public:

      unsigned char d_symbol_index;
      unsigned char d_frame_index;
      unsigned char d_superframe_index;

      //Transmission type parameters
      dvbt_transmission_mode_t d_transmission_mode;
      unsigned int d_Kmin;
      unsigned int d_Kmax;

      //Constelaltion parameters
      dvbt_constellation_t d_constellation;
      unsigned char d_m;

      //Inner Coding + puncturer parameters
      dvbt_code_rate_t d_code_rate_HP;
      dvbt_code_rate_t d_code_rate_LP;
      unsigned char d_cr_k;
      unsigned char d_cr_n;

      unsigned char d_include_cell_id;
      dvbt_hierarchy_t d_hierarchy;
      dvbt_guard_interval_t d_guard_interval;
      unsigned int d_cell_id;


      void set_frame_number(int fn);
      int get_frame_mumber();
      void set_constellation(dvbt_constellation_t constellation);
      dvbt_constellation_t get_constellation();
      void set_hierarchical(dvbt_hierarchy_t hierarchy);
      dvbt_hierarchy_t get_hierarchical();
      void set_code_rate_HP(dvbt_code_rate_t coderate);
      dvbt_code_rate_t get_code_rate_HP();
      void set_code_rate_LP(dvbt_code_rate_t coderate);
      dvbt_code_rate_t get_code_rate_LP();
      void set_transmission_mode(dvbt_transmission_mode_t transmission_mode);
      dvbt_transmission_mode_t get_transmission_mode();

      dvbt_config(dvbt_constellation_t constellation = gr::dvbt::QAM16, \
          dvbt_hierarchy_t hierarchy = gr::dvbt::NH, dvbt_code_rate_t code_rate_HP = gr::dvbt::C1_2, \
          dvbt_code_rate_t code_rate_LP = gr::dvbt::C1_2, dvbt_guard_interval_t guard_interval = gr::dvbt::G1_32, \
          dvbt_transmission_mode_t transmission_mode = gr::dvbt::T2k, int include_cell_id = 0, int cell_id = 0);
      ~dvbt_config();
    }; 
  } // namespace dvbt
} // namespace gr

    typedef gr::dvbt::dvbt_constellation_t dvbt_constellation_t;
    typedef gr::dvbt::dvbt_hierarchy_t dvbt_hierarchy_t;
    typedef gr::dvbt::dvbt_code_rate_t dvbt_code_rate_t;
    typedef gr::dvbt::dvbt_guard_interval_t dvbt_guard_interval_t;
    typedef gr::dvbt::dvbt_transmission_mode_t dvbt_transmission_mode_t;


#endif /* INCLUDED_DVBT_DVBT_CONFIG_H */

