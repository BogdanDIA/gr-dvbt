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

#ifndef INCLUDED_DVBT_CONFIG_H
#define INCLUDED_DVBT_CONFIG_H

namespace gr {
  namespace dvbt {

    typedef enum {
      QPSK = 0,
      QAM16,
      QAM64,
      RES
    } dvbt_constellation_t;

    typedef enum {
      NH = 0,
      ALPHA1,
      ALPHA2,
      ALPHA4,
      HRES1,
      HRES2,
      HRES3,
      HRES4
    } dvbt_hierarchy_t;

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
      T8k,
      TRES1,
      TRES2,
    };

    enum dvbt_guard_interval_t {
      G1_32 = 0,
      G1_16,
      G1_8,
      G1_4
    };

    class dvbt_config {

      public:
      
      unsigned int Kmin;
      unsigned int Kmax;
      unsigned int length_with_cellid;
      unsigned int length_without_cellid;
      unsigned char symbol_index;
      unsigned char superframe_index;
      unsigned char frame_index;
      unsigned char cell_id_included;
      dvbt_constellation_t constellation;
      dvbt_hierarchy_t hierarchy;
      dvbt_code_rate_t code_rate_HP;
      dvbt_code_rate_t code_rate_LP;
      dvbt_guard_interval_t guard_interval;
      dvbt_transmission_mode_t transmission_mode;
      unsigned int cell_id;


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
    };
  }
}

typedef gr::dvbt::dvbt_constellation_t dvbt_constellation_t;

#if 0
typedef dvbt_config::dvbt_constellation_t dvbt_constellation_t;
typedef dvbt_config::dvbt_hierarchy_t dvbt_hierarchy_t;
typedef dvbt_config::dvbt_code_rate_t dvbt_code_rate_t;
typedef dvbt_config::dvbt_guard_interval_t dvbt_guard_interval_t;
typedef dvbt_config::dvbt_transmission_mode_t dvbt_transmission_mode_t;
#endif

#endif /* INCLUDED_DVBT_CONFIG_H */

