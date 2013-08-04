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

#ifndef INCLUDED_DVBT_REFERENCE_SIGNALS_IMPL_H
#define INCLUDED_DVBT_REFERENCE_SIGNALS_IMPL_H

#include <dvbt/reference_signals.h>
#include <dvbt/dvbt_config.h>
#include <vector>
#include <deque>

    // This should eventually go into a const file
    const int SYMBOLS_PER_FRAME = 68;
    const int FRAMES_PER_SUPERFRAME = 4;

    const int SCATTERED_PILOT_SIZE_2k = 142; // TODO - fix this - what should be here
    const int CONTINUAL_PILOT_SIZE_2k = 45;
    const int TPS_PILOT_SIZE_2k = 17;

    const int SCATTERED_PILOT_SIZE_8k = 567; //TODO
    const int CONTINUAL_PILOT_SIZE_8k = 177;
    const int TPS_PILOT_SIZE_8k = 68;

namespace gr {
  namespace dvbt {
    /*!
     * \brief Pilot generator class.
     * Adds DVBT pilots to a signal.
     * \param c config object to keep all config data
     */

class pilot_gen {
  private:
    // this should be first in order to be initialized first
    const dvbt_config &config;

    // Debug only
    int d_index;

    int d_Kmin;
    int d_Kmax;
    int d_fft_length;
    int d_payload_length;
    int d_zeros_on_left;
    int d_zeros_on_right;
    int d_cp_length;

    static const int d_symbols_per_frame;
    static const int d_frames_per_superframe;

    // 2k mode
    // scattered pilot carriers info
    static const int d_spilot_carriers_size_2k;

    // continual pilot carriers info
    static const int d_cpilot_carriers_size_2k;
    static const int d_cpilot_carriers_2k[];

    // TPS carriers info
    static const int d_tps_carriers_size_2k;
    static const int d_tps_carriers_2k[];

    //8k mode
    // scattered pilot carriers info
    static const int d_spilot_carriers_size_8k;

    // continual pilot carriers info
    static const int d_cpilot_carriers_size_8k;
    static const int d_cpilot_carriers_8k[];

    // TPS carriers info
    static const int d_tps_carriers_size_8k;
    static const int d_tps_carriers_8k[];

    // TPS sync data
    static const int d_tps_sync_size;
    static const int d_tps_sync_even[];
    static const int d_tps_sync_odd[];

    // Variables to keep data for 2k, 8k, 4k
    int d_spilot_carriers_size;
    gr_complex * d_spilot_carriers_val;
    gr_complex * d_channel_gain;

    int d_cpilot_carriers_size;
    const int * d_cpilot_carriers;
    float * d_known_phase_diff;
    float * d_cpilot_phase_diff;
    int d_freq_offset_max;
    int d_freq_offset;
    float d_carrier_freq_correction;
    float d_sampling_freq_correction;

    int d_tps_carriers_size;
    const int * d_tps_carriers;
    gr_complex * d_tps_carriers_val;
 
    // Keeps TPS data
    unsigned char * d_tps_data;
    // Keep TPS carriers values from previous symbol
    gr_complex * d_prev_tps_symbol;
    // Keep TPS carriers values from current symbol
    gr_complex * d_tps_symbol;
    // Keeps the rcv TPS data, is a FIFO
    std::deque<char> d_rcv_tps_data;
    // Keeps the TPS sync sequence
    std::deque<char> d_tps_sync_evenv;
    std::deque<char> d_tps_sync_oddv;

    // Keeps channel estimation carriers
    // we use both continual and scattered carriers
    int * d_chanestim_carriers;

    // Keeps paload carriers
    int * d_payload_carriers;

    // Indexes for all carriers
    int d_spilot_index;
    int d_cpilot_index;
    int d_tpilot_index;
    int d_payload_index;
    int d_chanestim_index;
    int d_symbol_index;
    int d_prev_mod_symbol_index;
    int d_mod_symbol_index;
    int d_symbol_index_known;
    int d_equalizer_ready;
    int d_frame_index;
    int d_superframe_index;

    // PRPS generator data buffer
    char * d_wk;
    // Generate PRBS
    void generate_prbs();

    // TPS private methods
    void set_tps_bits(int start, int stop, unsigned int data);
     
    void set_symbol_index(int index);
    int get_symbol_index();
    void set_tps_data();
    void get_tps_data();

    void reset_pilot_generator();

    // Scattered pilot generator methods
    int get_current_spilot(int spilot) const;
    gr_complex get_spilot_value(int spilot);
    void set_spilot_value(int spilot, gr_complex val);
    void advance_spilot(int sindex);
    // Methods used to quick iterate through all spilots
    int get_first_spilot();
    int get_last_spilot() const;
    int get_next_spilot();
    // Scattered pilot data processing method
    int process_spilot_data(const gr_complex * in);

    // Channel estimation methods
    void set_channel_gain(int spilot, gr_complex val);

    // Continual pilot generator methods
    int get_current_cpilot() const;
    gr_complex get_cpilot_value(int cpilot);
    void advance_cpilot();
    // Continual pilot data processing methods
    void process_cpilot_data(const gr_complex * in);
    gr_complex frequency_correction();
    void compute_oneshot_csft(const gr_complex * in);

    // TPS generator methods
    int get_current_tpilot() const;
    gr_complex get_tpilot_value(int tpilot);
    void advance_tpilot();
    // TPS data
    void format_tps_data();
    // Encode TPS data
    void generate_bch_code();
    // Verify parity on TPS data
    int verify_bch_code(std::deque<char> data);
    // TPS data processing metods
    int process_tps_data(const gr_complex * in, const int diff_symbo_index);

    // Channel estimation methods
    void set_chanestim_carrier(int k);

    // Payload data processing methods
    int get_current_payload();
    void advance_chanestim();
    void set_payload_carrier(int k);
    void advance_payload();
    void process_payload_data(const gr_complex *in, gr_complex *out);

    int d_trigger_index;

  public:

    pilot_gen(const dvbt_config &config);
    ~pilot_gen();

    /*!
     * ETSI EN 300 744 Clause 4.5. \n
     * Update a set of carriers with the pilot signals. \n
     */
    void update_output(const gr_complex *in, gr_complex *out);

    /*!
     * TODO
     * ETSI EN 300 744 Clause 4.5. \n
     * Extract data from a set of carriers using pilot signals. \n
     * This is doing frequency correcton, equalization. \n
     */
    int parse_input(const gr_complex *in, const unsigned char *trigger_in, gr_complex *out, unsigned char *trigger_out);

    };


    /*!
     * \brief Reference signals class. \n
     * \ingroup dvbt
     * \param itemsize size of an in/out item
     * \param ninput input stream length
     * \param noutput output stream length
     * \param constellation constelaltion used
     * \param hierarchy hierarchy used
     * \param code_rate_HP high priority stream code rate
     * \param code_rate_LP low priority stream code rate
     * \param guard_interval gurad interval used
     * \param transmission_mode transmission mode used
     * \param include_cell_id include or not Cell ID
     * \param cell_id value of the Cell ID
     */
    class reference_signals_impl : public reference_signals
    {
      // configuration object for this class
      const dvbt_config config;

    private:
	// Pilot Generator object
    	pilot_gen d_pg;

	//In and Out data length
	int d_ninput;
	int d_noutput;

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

