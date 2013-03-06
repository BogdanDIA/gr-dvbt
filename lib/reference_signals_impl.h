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

#ifndef INCLUDED_DVBT_REFERENCE_SIGNALS_IMPL_H
#define INCLUDED_DVBT_REFERENCE_SIGNALS_IMPL_H

#include <dvbt/reference_signals.h>
#include <dvbt/dvbt_config.h>
#include <vector>

// This should eventually go into a const file
    const int SYMBOLS_PER_FRAME = 68;
    const int FRAMES_PER_SUPERFRAME = 4;
    const int SCATTERED_PILOT_SIZE_2k = 142;
    const int CONTINUAL_PILOT_SIZE_2k = 45;
    const int TPS_PILOT_SIZE_2k = 17;

namespace gr {
  namespace dvbt {

class pilot_gen {
  public:
    // this should be first in order to be initialized first
    const dvbt_config &config;

    static const int d_symbols_per_frame;
    static const int d_frames_per_superframe;

    // scattered pilot carriers info
    static const int d_spilot_carriers_size_2k;

    // continual pilot carriers info
    static const int d_cpilot_carriers_size_2k;
    static const int d_cpilot_carriers_2k[];

    // TPS carriers info
    static const int d_tps_carriers_size_2k;
    static const int d_tps_carriers_2k[];
    gr_complex * d_tps_carriers_val_2k;

    // Keeps TPS data
    unsigned char * d_tps_data;


    int d_spilot_index;
    int d_cpilot_index;
    int d_tpilot_index;
    int d_symbol_index;
    int d_frame_index;
    int d_superframe_index;

    // PRPS generator data buffer
    char * d_wk;
    // Generate PRBS
    void generate_prbs();

    // TPS private methods
    void set_tps_bits(int start, int stop, unsigned int data);
    
  public:

    pilot_gen(const dvbt_config &config);
    ~pilot_gen();

    void set_symbol_index(int);
    int get_symbol_index();
    void set_tps_data();
    void get_tps_data();


    void reset_pilot_generator();
    // Scattered pilot generator methods
    int get_current_spilot() const;
    gr_complex get_spilot_value(int spilot);
    void advance_spilot();
    // Continual pilot denerator methods
    int get_current_cpilot() const;
    gr_complex get_cpilot_value(int cpilot);
    void advance_cpilot();
    // TPS generator methods
    int get_current_tpilot() const;
    gr_complex get_tpilot_value(int tpilot);
    void advance_tpilot();
    // TPS data
    void format_tps_data();
    // Encode TPS data
    void generate_bch_code();

    // Update out
    void update_output(const gr_complex *in, gr_complex *out);
};


    class reference_signals_impl : public reference_signals
    {
      // configuration object for this class
      const dvbt_config config;

    private:
	// Pilot Generator object
    	pilot_gen d_pg;

    public:

      
      reference_signals_impl(int itemsize, int ninput, int noutput, \
	dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, \
	dvbt_code_rate_t code_rate_HP, dvbt_code_rate_t code_rate_LP, \
	dvbt_guard_interval_t guard_interval, \
	dvbt_transmission_mode_t transmission_mode = gr::dvbt::T2k, int include_cell_id = 0, int cell_id = 0);
      ~reference_signals_impl();
      
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_REFERENCE_SIGNALS_IMPL_H */

