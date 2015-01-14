/* -*- c++ -*- */
/* 
 * Copyright 2013,2014,2015 <Bogdan Diaconescu, yo3iiu@yo3iiu.ro>.
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

#ifndef INCLUDED_DVBT_DEMOD_REFERENCE_SIGNALS_IMPL_H
#define INCLUDED_DVBT_DEMOD_REFERENCE_SIGNALS_IMPL_H

#include <dvbt/demod_reference_signals.h>
#include <dvbt/reference_signals.h>
#include "reference_signals_impl.h"

namespace gr {
  namespace dvbt {

    class demod_reference_signals_impl : public demod_reference_signals
    {
      // configuration object for this class
      const dvbt_config config;

    private:
      // Pilot Generator object
      pilot_gen d_pg;

      //In and Out data length
      int d_ninput;
      int d_noutput;

      int d_init;
      int d_fi_start;

      int is_sync_start(int nitems);

    public:
      demod_reference_signals_impl(int itemsize, int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, \
        dvbt_code_rate_t code_rate_HP, dvbt_code_rate_t code_rate_LP, \
        dvbt_guard_interval_t guard_interval, \
        dvbt_transmission_mode_t transmission_mode = gr::dvbt::T2k, int include_cell_id = 0, int cell_id = 0);
      ~demod_reference_signals_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_DEMOD_REFERENCE_SIGNALS_IMPL_H */

