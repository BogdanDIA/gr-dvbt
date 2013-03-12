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
#include "reference_signals_impl.h"
#include <complex>
#include <stdio.h>

namespace gr {
  namespace dvbt {

   //Number of symbols in a frame
    const int pilot_gen::d_symbols_per_frame = SYMBOLS_PER_FRAME;
    //Number of frames in a superframe
    const int pilot_gen::d_frames_per_superframe = FRAMES_PER_SUPERFRAME;

    // Scattered pilots # of carriers
    const int pilot_gen::d_spilot_carriers_size_2k = SCATTERED_PILOT_SIZE_2k;
    // Continual pilots # of carriers and positions
    const int pilot_gen::d_cpilot_carriers_size_2k = CONTINUAL_PILOT_SIZE_2k;
    const int pilot_gen::d_cpilot_carriers_2k[pilot_gen::d_cpilot_carriers_size_2k] = {
       0,   48,   54,   87,  141,  156,  192, \
     201,  255,  279,  282,  333,  432,  450, \
     483,  525,  531,  618,  636,  714,  759, \
     765,  780,  804,  873,  888,  918,  939, \
     942,  969,  984, 1050, 1101, 1107, 1110, \
    1137, 1140, 1146, 1206, 1269, 1323, 1377, \
    1491, 1683, 1704
    };
    // TPS pilots # of carriers and positions
    const int pilot_gen::d_tps_carriers_size_2k = TPS_PILOT_SIZE_2k;
    const int pilot_gen::d_tps_carriers_2k[pilot_gen::d_tps_carriers_size_2k] = {
      34,   50,  209,  346,  413, \
     569,  595,  688,  790,  901, \
    1073, 1219, 1262, 1286, 1469, \
    1594, 1687
    };

    /*
     * Constructor of class
     */
    pilot_gen::pilot_gen(const dvbt_config &c) : config(c),
                                        d_spilot_index(0),
                                        d_cpilot_index(0),
                                        d_tpilot_index(0),
                                        d_symbol_index(0),
                                        d_frame_index(0),
                                        d_superframe_index(0)
    {
      //allocate PRBS buffer
      d_wk = new char[config.d_Kmax - config.d_Kmin + 1];
      if (d_wk == NULL)
      {
        std::cout << "error allocating d_wk" << std::endl;
        return;
      }

      // allocate buffer for first tps symbol constellation
      d_tps_carriers_val_2k = new gr_complex[d_tps_carriers_size_2k];
      if (d_tps_carriers_val_2k == NULL)
      {
        std::cout << "error allocating d_tps_carriers_val_2k" << std::endl;
        return;
      }

      // allocate tps data buffer
      d_tps_data = new unsigned char[d_symbols_per_frame];
      if (d_tps_data == NULL)
      {
        std::cout << "error allocating d_tps_data" << std::endl;
        return;
      }
      
      // Generate wk sequence
      generate_prbs();
      // Reset the pilot generator
      reset_pilot_generator();
      // Format TPS data with current values
      format_tps_data();
    }

    /*
     * Destructor of class
     */
    pilot_gen::~pilot_gen()
    {
      delete d_wk;
      delete d_tps_carriers_val_2k;
      delete d_tps_data;
    }

    /*
     * Generate PRBS sequence
     * X^11 + X^2 + 1
     * en 300 744 - section 4.5.2
     */
    void
    pilot_gen::generate_prbs()
    {  
      // init PRBS register with 1s
      unsigned int reg_prbs = (1 << 11) - 1;

      for (int k = 0; k < (config.d_Kmax - config.d_Kmin + 1); k++)
      {
        d_wk[k] = (char)(reg_prbs & 0x01);
        int new_bit = ((reg_prbs >> 2) ^ (reg_prbs >> 0)) & 0x01;
        reg_prbs = (reg_prbs >> 1) | (new_bit << 10);
      }
    }

    /*
     * Generate shortened BCH(67, 53) codes from TPS data
     * Extend the code with 60 bits and use BCH(127, 113)
     */
    void
    pilot_gen::generate_bch_code()
    {
      //TODO
      //DO other way: if (feedback == 1) reg = reg ^ polymomial
      //else nothing

      //(n, k) = (127, 113) = (60+67, 60+53)
      unsigned int reg_bch = 0;
      unsigned char data_in[113];

      //fill in 60 zeros
      memset(&data_in[0], 0, 60);
      //fill in TPS data - start bit not included
      memcpy(&data_in[60], &d_tps_data[1], 53);

      //X^14+X^9+X^8+X^6+X^5+X^4+X^2+X+1
      for (int i = 0; i < 113; i++)
      {
        int feedback = 0x1 & (data_in[i] ^ reg_bch);
        reg_bch  = reg_bch >> 1;
        reg_bch |= feedback << 13;
        reg_bch = reg_bch \
            ^ (feedback << 12) ^ (feedback << 11) \
            ^ (feedback << 9) ^ (feedback << 8) \
            ^ (feedback << 7) ^ (feedback << 5) \
            ^ (feedback << 4);
      }

      for (int i = 0; i < 14; i++)
        d_tps_data[i + 54] = 0x1 & (reg_bch >> i); 
    }

    void
    pilot_gen::set_symbol_index(int sindex)
    {
      d_symbol_index = sindex;
    }

    int
    pilot_gen::get_symbol_index()
    {
      return d_symbol_index;
    }

    void
    pilot_gen::set_tps_data()
    {
    }

    void
    pilot_gen::get_tps_data()
    {
    }

    /*
     * Reset pilot generator
     */
    void
    pilot_gen::reset_pilot_generator()
    {
      d_spilot_index = 0; d_cpilot_index = 0; d_tpilot_index = 0;
      d_symbol_index = 0; d_frame_index = 0; d_superframe_index = 0;
    }

    /*
     * Init scattered pilot generator
     */
    int
    pilot_gen::get_current_spilot() const
    {
      //TODO - can be optimised for same symbol_index
      return (config.d_Kmin + 3 * (d_symbol_index % 4) + 12 * d_spilot_index);
    }

    gr_complex
    pilot_gen::get_spilot_value(int spilot)
    {
      // TODO - can be calculated at the beginning
      return gr_complex(4 * 2 * (0.5 - d_wk[spilot]) / 3, 0);
    }

    void
    pilot_gen::advance_spilot()
    {
      d_spilot_index = (++d_spilot_index) % d_spilot_carriers_size_2k;
    }

    /*
     * Init continual pilot generator
     */
    int 
    pilot_gen::get_current_cpilot() const
    {
      return d_cpilot_carriers_2k[d_cpilot_index];
    }

    gr_complex
    pilot_gen::get_cpilot_value(int cpilot)
    {
      //TODO - can be calculated at the beginning
      return gr_complex((float)(4 * 2 * (0.5 - d_wk[cpilot])) / 3, 0);
    }

    void
    pilot_gen::advance_cpilot()
    {
      d_cpilot_index = (++d_cpilot_index) % d_cpilot_carriers_size_2k;
    }

    /*
     * Init tps sequence, return values for firs position
     * If first symbol then init tps DBPSK data
     */
    int
    pilot_gen::get_current_tpilot() const
    {
        return d_tps_carriers_2k[d_tpilot_index];
    }

    gr_complex
    pilot_gen::get_tpilot_value(int tpilot)
    {
      //TODO - it can be calculated at the beginnning
      if (d_symbol_index == 0)
        d_tps_carriers_val_2k[d_tpilot_index] = gr_complex(2 * (0.5 - d_wk[tpilot]), 0);
      else
        if (d_tps_data[d_symbol_index] == 1)
            d_tps_carriers_val_2k[d_tpilot_index] = gr_complex(-d_tps_carriers_val_2k[d_tpilot_index].real(), 0);

      return d_tps_carriers_val_2k[d_tpilot_index];
    }

    void
    pilot_gen::advance_tpilot()
    {
      d_tpilot_index = (++d_tpilot_index) % (d_tps_carriers_size_2k);
    }

    /*
     * Set a number of bits to a specified value
     */
    void
    pilot_gen::set_tps_bits(int start, int stop, unsigned int data)
    {
        for (int i = start; i >= stop; i--)
        {
          d_tps_data[i] = data & 0x1;
          data = data >> 1;
        }
    }

    /*
     * Format data that will be sent with TPS signals
     * en 300 744 - section 4.6.2
     * s0 Initialization
     * s1-s16 Synchronization word
     * s17-s22 Length Indicator
     * s23-s24 Frame Number
     * S25-s26 Constellation
     * s27, s28, s29 Hierarchy information
     * s30, s31, s32 Code rate, HP stream
     * s33, s34, s35 Code rate, LP stream
     * s36, s37 Guard interval
     * s38, s39 Transmission mode
     * s40, s47 Cell identifier
     * s48-s53 All set to "0"
     * s54-s67 Error protection (BCH code)
     */
    void
    pilot_gen::format_tps_data()
    {
      //Clause 4.6.3
      set_tps_bits(0, 0, d_wk[0]);
      //Clause 4.6.2.2
      if (d_frame_index % 2)
        set_tps_bits(16, 1, 0xca11);
      else
        set_tps_bits(16, 1, 0x35ee);
      //Clause 4.6.2.3
      if (config.d_include_cell_id)
        set_tps_bits(22, 17, 0x1f);
      else
        set_tps_bits(22, 17, 0x17);
      //Clause 4.6.2.4
      set_tps_bits(24, 23, d_frame_index);
      //Clause 4.6.2.5
      set_tps_bits(26, 25, config.d_constellation);
      //Clause 4.6.2.6
      set_tps_bits(29, 27, config.d_hierarchy);
      //Clause 4.6.2.7
      set_tps_bits(32, 30, config.d_code_rate_HP);
      set_tps_bits(35, 33, config.d_code_rate_LP);
      //Clause 4.6.2.8
      set_tps_bits(37, 36, config.d_guard_interval);
      //Clause 4.6.2.9
      set_tps_bits(39, 38, config.d_transmission_mode);
      //Clause 4.6.2.10
      set_tps_bits(47, 40, config.d_cell_id);
      //This bits are set to zero
      set_tps_bits(53, 48, 0);
      //Clause 4.6.2.11
      generate_bch_code();
}

    void
    pilot_gen::update_output(const gr_complex *in, gr_complex *out)
    {
      int is_payload = 1;
      int payload_count = 0;

      //move to the next symbol
      //re-genereate TPS data
      format_tps_data();

      //reset indexes
      payload_count = 0;
      d_spilot_index = 0; d_cpilot_index = 0; d_tpilot_index = 0;

      //process one block - one symbol
      for (int k = 0; k < (config.d_Kmax - config.d_Kmin + 1); k++)
      {
          is_payload = 1;
          if (k == get_current_spilot())
          {
            out[k] = get_spilot_value(k);
            advance_spilot();
            is_payload = 0;
          }

          if (k == get_current_cpilot())
          {
            out[k] = get_cpilot_value(k);
            advance_cpilot();
            is_payload = 0;
          }

          if (k == get_current_tpilot())
          {
            out[k] = get_tpilot_value(k);
            advance_tpilot();
            is_payload = 0;
          } 

          if (is_payload == 1)
            out[k] = in[payload_count++];
      }
      // update indexes
      if (++d_symbol_index == d_symbols_per_frame)
      {
        d_symbol_index = 0;
        if (++d_frame_index == d_frames_per_superframe)
        {
          d_frame_index == 0;
          d_superframe_index++;
        }
      }
    }

    reference_signals::sptr
    reference_signals::make(int itemsize, int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, \
        dvbt_code_rate_t code_rate_HP, dvbt_code_rate_t code_rate_LP, \
        dvbt_guard_interval_t guard_interval, dvbt_transmission_mode_t transmission_mode, \
        int include_cell_id, int cell_id)
    {
      return gnuradio::get_initial_sptr (new reference_signals_impl(itemsize, ninput, \
            noutput, constellation, hierarchy, code_rate_HP, code_rate_LP, \
            guard_interval, transmission_mode, include_cell_id, cell_id));
    }

    /*
     * The private constructor
     */
    reference_signals_impl::reference_signals_impl(int itemsize, int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, dvbt_code_rate_t code_rate_HP,\
          dvbt_code_rate_t code_rate_LP, dvbt_guard_interval_t guard_interval,\
          dvbt_transmission_mode_t transmission_mode, int include_cell_id, int cell_id)
      : gr_block("reference_signals",
		      gr_make_io_signature(1, 1, itemsize * ninput),
		      gr_make_io_signature(1, 1, itemsize * noutput)),
		      config(ninput, noutput, constellation, hierarchy, code_rate_HP, code_rate_LP, \
              guard_interval, transmission_mode, include_cell_id, cell_id),
		      d_pg(config)
    {
      //
      //
    }

    /*
     * Our virtual destructor.
     */
    reference_signals_impl::~reference_signals_impl()
    {
    }

    
    void
    reference_signals_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
        ninput_items_required[0] = noutput_items;
    }

    int
    reference_signals_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const gr_complex *in = (const gr_complex *) input_items[0];
        gr_complex *out = (gr_complex *) output_items[0];

        {
          for (int i = 0; i < noutput_items; i++)
            d_pg.update_output(&in[i * config.d_ninput], &out[i * config.d_noutput]);
        }

        // Do <+signal processing+>
        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */


