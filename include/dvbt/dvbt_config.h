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

/*! 
 * \file dvbt_config.h 
 * \brief DVBT config header.
 */

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
     * \brief DVBT Configuration class.
     * \ingroup dvbt
     * Keeps all data related to configuration. \n
     * \param constellation constelaltion used
     * \param hierarchy hierarchy used
     * \param code_rate_HP high priority stream code rate
     * \param code_rate_LP low priority stream code rate
     * \param guard_interval gurad interval used
     * \param transmission_mode transmission mode used
     * \param include_cell_id include or not Cell ID
     * \param cell_id value of the Cell ID
     */
    class DVBT_API dvbt_config
    {

      public:

      int d_symbols_per_frame;
      int d_frames_per_superframe;

      int d_symbol_index;
      int d_frame_index;
      int d_superframe_index;

      //Transmission type parameters
      dvbt_transmission_mode_t d_transmission_mode;
      int d_Kmin;
      int d_Kmax;
      int d_fft_length;
      int d_payload_length;
      int d_zeros_on_left;
      int d_cp_length;

      //Constelaltion parameters
      dvbt_constellation_t d_constellation;
      int d_constellation_size;
      int d_step;
      int d_m;
      float d_norm;

      //Inner Coding + puncturer parameters
      dvbt_code_rate_t d_code_rate_HP;
      dvbt_code_rate_t d_code_rate_LP;
      int d_cr_k;
      int d_cr_n;

      //Hierarchy information
      dvbt_hierarchy_t d_hierarchy;
      int d_alpha;

      dvbt_guard_interval_t d_guard_interval;
      int d_include_cell_id;
      int d_cell_id;

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

