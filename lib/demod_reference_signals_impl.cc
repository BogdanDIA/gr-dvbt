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

#include <gnuradio/io_signature.h>
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

    int
    demod_reference_signals_impl::is_sync_start(int nitems)
    {
      std::vector<tag_t> tags;
      const uint64_t nread = this->nitems_read(0); //number of items read on port 0
      this->get_tags_in_range(tags, 0, nread, nread + nitems, pmt::string_to_symbol("sync_start"));

      return tags.size() ? 1 : 0;
    }

    /*
     * The private constructor
     */
    demod_reference_signals_impl::demod_reference_signals_impl(int itemsize, int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, dvbt_code_rate_t code_rate_HP,\
          dvbt_code_rate_t code_rate_LP, dvbt_guard_interval_t guard_interval,\
          dvbt_transmission_mode_t transmission_mode, int include_cell_id, int cell_id)
      : block("demod_reference_signals",
          io_signature::make(1, 1, itemsize * ninput),
          io_signature::make(1, 1, itemsize * noutput)),
          config(constellation, hierarchy, code_rate_HP, code_rate_LP, \
            guard_interval, transmission_mode, include_cell_id, cell_id),
          d_ninput(ninput), d_noutput(noutput),
          d_pg(config)
    {
      //
      d_skip = 1;
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
      gr_complex *out = (gr_complex *) output_items[0];

      int symbol_index;
      int to_out = 0;

      for (int i = 0; i < noutput_items; i++)
        to_out += d_pg.parse_input(&in[i * d_ninput], &out[i * d_noutput], &symbol_index);

      /*
       * Wait for a sync_start tag from upstream that signals when to start.
       * Skip OFDM symbols to start always with quad symbol (0, 4, 8, etc)
       */
      if (is_sync_start(noutput_items))
        d_skip = 1;

      if (d_skip)
      {
        if ((symbol_index % 4) != 0)
        {
          consume_each(1);
          return (0);
        }
        else
          d_skip = 0;
      }

      // Send a tag for each OFDM symbol informing about
      // symbol index.
      const uint64_t offset = this->nitems_written(0);
      pmt::pmt_t key = pmt::string_to_symbol("symbol_index");
      pmt::pmt_t value = pmt::from_long(symbol_index);
      this->add_item_tag(0, offset, key, value);

      // Consume from input stream
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return to_out;
    }

  } /* namespace dvbt */
} /* namespace gr */

