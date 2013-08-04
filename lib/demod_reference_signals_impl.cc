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
#include "demod_reference_signals_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    demod_reference_signals::sptr
    demod_reference_signals::make(int itemsize, int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, \
        dvbt_code_rate_t code_rate_HP, dvbt_code_rate_t code_rate_LP, \
        dvbt_guard_interval_t guard_interval, dvbt_transmission_mode_t transmission_mode, \
        int include_cell_id, int cell_id)
    {
      return gnuradio::get_initial_sptr (new demod_reference_signals_impl(itemsize, ninput, \
            noutput, constellation, hierarchy, code_rate_HP, code_rate_LP, \
            guard_interval, transmission_mode, include_cell_id, cell_id));
    }

    /*
     * The private constructor
     */
    demod_reference_signals_impl::demod_reference_signals_impl(int itemsize, int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, dvbt_code_rate_t code_rate_HP,\
          dvbt_code_rate_t code_rate_LP, dvbt_guard_interval_t guard_interval,\
          dvbt_transmission_mode_t transmission_mode, int include_cell_id, int cell_id)
      : gr_block("demod_reference_signals",
          gr_make_io_signature2(2, 2, itemsize * ninput, ninput),
		      gr_make_io_signature2(2, 2, itemsize * noutput, noutput)), // TODO trigger is char based stream!
		      config(constellation, hierarchy, code_rate_HP, code_rate_LP, \
            guard_interval, transmission_mode, include_cell_id, cell_id),
          d_ninput(ninput), d_noutput(noutput),
		      d_pg(config)
    {
      //
    }

    /*
     * Our virtual destructor.
     */
    demod_reference_signals_impl::~demod_reference_signals_impl()
    {
    }

    void
    demod_reference_signals_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      int ninputs = ninput_items_required.size();

      for (int i = 0; i < ninputs; i++)
        ninput_items_required[i] =  2 * noutput_items;
    }

    int
    demod_reference_signals_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      const unsigned char *trigger_in = (const unsigned char *) input_items[1];
      gr_complex *out = (gr_complex *) output_items[0];
      unsigned char *trigger_out = (unsigned char *) output_items[1];

      int to_out = 0;

      for (int i = 0; i < noutput_items; i++)
        to_out += d_pg.parse_input(&in[i * d_ninput], &trigger_in[i * d_ninput], &out[i * d_noutput], &trigger_out[i * d_noutput]);

      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return to_out;
    }

  } /* namespace dvbt */
} /* namespace gr */

