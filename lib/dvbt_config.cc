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

    dvbt_config::dvbt_config(int ninput, int noutput, dvbt_constellation_t constellation, \
	dvbt_hierarchy_t hierarchy, dvbt_code_rate_t code_rate_HP, \
	dvbt_code_rate_t code_rate_LP, dvbt_guard_interval_t guard_interval, \
	dvbt_transmission_mode_t transmission_mode, int include_cell_id, int cell_id) :
	    d_ninput(ninput), d_noutput(noutput),
	    d_constellation(constellation), d_hierarchy(hierarchy), d_code_rate_HP(code_rate_HP),
	    d_code_rate_LP(code_rate_LP), d_guard_interval(guard_interval), d_transmission_mode(transmission_mode),
	    d_include_cell_id(include_cell_id), d_cell_id(cell_id)
    {
      switch (d_transmission_mode)
      {
	      case gr::dvbt::T2k:
		      d_Kmin = 0; d_Kmax = 1704;
		      break;
	      case gr::dvbt::T8k:
		      d_Kmin = 0; d_Kmax = 6816;
		      break;
	      default:
		      d_Kmin = 0; d_Kmax = 1704;
		      break;
      }

      switch (d_constellation)
      {
	      case gr::dvbt::QPSK:
		      d_m = 2;
		      break;
	      case gr::dvbt::QAM16:
		      d_m = 4;
		      break;
	      case gr::dvbt::QAM64:
		      d_m = 6;
		      break;
	      default:
		      d_m = 4;
		      break;
      }

      switch (d_code_rate_HP)
      {
	      case gr::dvbt::C1_2:
		d_k = 1; d_n = 2;
		break;
	      case gr::dvbt::C2_3:
      		d_k = 2; d_n = 3;
      		break;
	      case gr::dvbt::C3_4:
	        d_k = 3; d_n = 4;
		break;
	      case gr::dvbt::C5_6:
		d_k = 5; d_k = 6;
		break;
	      case gr::dvbt::C7_8:
		d_k = 7; d_n = 8;
		break;
	      default:
		d_k = 1; d_n = 2;
		break;
      }

      switch (d_code_rate_LP)
      {
	      case gr::dvbt::C1_2:
		d_k = 1; d_n = 2;
		break;
	      case gr::dvbt::C2_3:
      		d_k = 2; d_n = 3;
      		break;
	      case gr::dvbt::C3_4:
	        d_k = 3; d_n = 4;
		break;
	      case gr::dvbt::C5_6:
		d_k = 5; d_k = 6;
		break;
	      case gr::dvbt::C7_8:
		d_k = 7; d_n = 8;
		break;
	      default:
		d_k = 1; d_n = 2;
		break;
      }
    }

    dvbt_config::~dvbt_config()
    {
    }

  } /* namespace dvbt */
} /* namespace gr */

