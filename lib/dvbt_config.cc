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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gr_io_signature.h>
#include <dvbt/dvbt_config.h>
#include <iostream>
#include <stdio.h>

namespace gr {
  namespace dvbt {

    void 
    dvbt_config::set_frame_number(int fn)
    {
      d_frame_index = fn;
    }
    int 
    dvbt_config::get_frame_mumber()
    {
      return (d_frame_index);
    }
    void
    dvbt_config::set_constellation(dvbt_constellation_t constellation)
    {
      d_constellation = constellation; 
    }
    dvbt_constellation_t 
    dvbt_config::get_constellation() 
    {
      return (d_constellation);
    }
    void
    dvbt_config::set_hierarchical(dvbt_hierarchy_t hierarchy)
    {
      d_hierarchy = hierarchy;
    }
    dvbt_hierarchy_t
    dvbt_config::get_hierarchical()
    {
      d_hierarchy;
    }
    void
    dvbt_config::set_code_rate_HP(dvbt_code_rate_t code_rate)
    {
      d_code_rate_HP = code_rate;
    }
    void
    dvbt_config::set_code_rate_LP(dvbt_code_rate_t code_rate)
    {
      d_code_rate_LP = code_rate;
    }
    dvbt_code_rate_t
    dvbt_config::get_code_rate_HP()
    {
      return d_code_rate_HP;
    }
    dvbt_code_rate_t
    dvbt_config::get_code_rate_LP()
    {
      return d_code_rate_LP;
    }
    void
    dvbt_config::set_transmission_mode(dvbt_transmission_mode_t transmission_mode)
    {
      d_transmission_mode = transmission_mode;
    }
    dvbt_transmission_mode_t
    dvbt_config::get_transmission_mode()
    {
      return d_transmission_mode;
    }

    dvbt_config::dvbt_config(dvbt_constellation_t constellation, \
	dvbt_hierarchy_t hierarchy, dvbt_code_rate_t code_rate_HP, \
	dvbt_code_rate_t code_rate_LP, dvbt_guard_interval_t guard_interval, \
	dvbt_transmission_mode_t transmission_mode, int include_cell_id, int cell_id) :
	    d_constellation(constellation), d_hierarchy(hierarchy), d_code_rate_HP(code_rate_HP),
	    d_code_rate_LP(code_rate_LP), d_guard_interval(guard_interval), d_transmission_mode(transmission_mode),
	    d_include_cell_id(include_cell_id), d_cell_id(cell_id)
    {
      d_symbols_per_frame = 68;
      d_frames_per_superframe = 4;

      switch (d_transmission_mode)
      {
	case gr::dvbt::T2k:
	  d_Kmin = 0; d_Kmax = 1704;
	  d_fft_length = 2048;
	  d_payload_length = 1512;
	  break;
	case gr::dvbt::T8k:
	  d_Kmin = 0; d_Kmax = 6816;
	  d_fft_length = 8192;
	  d_payload_length = 6048;
	  break;
	default:
	  d_Kmin = 0; d_Kmax = 1704;
	  d_fft_length = 2048;
	  d_payload_length = 1512;
	  break;
      }
      d_zeros_on_left = int(ceil((d_fft_length - (d_Kmax - d_Kmin + 1)) / 2.0));
      d_zeros_on_right = d_fft_length - d_zeros_on_left - (d_Kmax - d_Kmin + 1);

      switch (d_constellation)
      {
	case gr::dvbt::QPSK:
	  d_constellation_size = 2;
	  d_step = 2;
	  d_m = 2;
	  break;
	case gr::dvbt::QAM16:
	  d_constellation_size = 16;
	  d_step = 2;
	  d_m = 4;
	  break;
	case gr::dvbt::QAM64:
	  d_constellation_size = 64;
	  d_step = 2;
	  d_m = 6;
	  break;
	default:
	  d_constellation_size = 16;
	  d_step = 2;
	  d_m = 4;
	  break;
      }

      switch (d_code_rate_HP)
      {
        case gr::dvbt::C1_2:
          d_cr_k = 1; d_cr_n = 2;
          break;
        case gr::dvbt::C2_3:
          d_cr_k = 2; d_cr_n = 3;
          break;
        case gr::dvbt::C3_4:
          d_cr_k = 3; d_cr_n = 4;
          break;
        case gr::dvbt::C5_6:
          d_cr_k = 5; d_cr_k = 6;
        break;
        case gr::dvbt::C7_8:
          d_cr_k = 7; d_cr_n = 8;
          break;
        default:
          d_cr_k = 1; d_cr_n = 2;
          break;
      }

      switch (d_code_rate_LP)
      {
        case gr::dvbt::C1_2:
          d_cr_k = 1; d_cr_n = 2;
          break;
        case gr::dvbt::C2_3:
          d_cr_k = 2; d_cr_n = 3;
          break;
        case gr::dvbt::C3_4:
          d_cr_k = 3; d_cr_n = 4;
          break;
        case gr::dvbt::C5_6:
          d_cr_k = 5; d_cr_k = 6;
        break;
        case gr::dvbt::C7_8:
          d_cr_k = 7; d_cr_n = 8;
          break;
        default:
          d_cr_k = 1; d_cr_n = 2;
          break;
      }

      switch (guard_interval)
      {
        case gr::dvbt::G1_32:
          d_cp_length = d_fft_length / 32;
          break;
        case gr::dvbt::G1_16:
          d_cp_length = d_fft_length / 16;
          break;
        case gr::dvbt::G1_8:
          d_cp_length = d_fft_length / 8;
          break;
        case gr::dvbt::G1_4:
          d_cp_length = d_fft_length / 4;
          break;
        default:
          d_cp_length = d_fft_length / 32;
          break;
      }

      switch (d_hierarchy)
      {
        case (gr::dvbt::NH):
          d_alpha = 1; break;
        case (gr::dvbt::ALPHA1):
          d_alpha = 1; break;
        case (gr::dvbt::ALPHA2):
          d_alpha = 2; break;
        case (gr::dvbt::ALPHA4):
          d_alpha = 4;break;
        default:
          d_alpha = 1; break;
      }

      // ETSI EN 400 744 Clause 4.4
      // Normalization factor
      switch (d_m)
      {
	case 2:
	  d_norm = 1.0 / sqrt(2);
	  break;
	case 16:
	  if (d_alpha == 1) d_norm = 1.0 / sqrt(10);
	  if (d_alpha == 2) d_norm = 1.0 / sqrt(20);
	  if (d_alpha == 4) d_norm = 1.0 / sqrt(52);
	  break;
	case 64:
	  if (d_alpha == 1) d_norm = 1.0 / sqrt(42);
	  if (d_alpha == 2) d_norm = 1.0 / sqrt(60);
	  if (d_alpha == 4) d_norm = 1.0 / sqrt(108);
	  break;
	default:
	  if (d_alpha == 1) d_norm = 1.0 / sqrt(10);
	  if (d_alpha == 2) d_norm = 1.0 / sqrt(20);
	  if (d_alpha == 4) d_norm = 1.0 / sqrt(52);
	  break;
      }
    }

    dvbt_config::~dvbt_config()
    {
    }

  } /* namespace dvbt */
} /* namespace gr */

