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

#ifndef INCLUDED_DVBT_BIT_INNER_INTERLEAVER_IMPL_H
#define INCLUDED_DVBT_BIT_INNER_INTERLEAVER_IMPL_H

#include <dvbt/bit_inner_interleaver.h>
#include <dvbt/dvbt_config.h>

namespace gr {
  namespace dvbt {
    /*!
     * \brief Bit Inner interleaver
     * \ingroup dvbt
     * \param ninput length of input stream \n
     * \param noutput length of output stream \n
     * \param constellation constelaltion used \n
     * \param hierarchy hierarchy used \n
     */
    class bit_inner_interleaver_impl : public bit_inner_interleaver
    {
    private:
      const dvbt_config config;

      int d_nsize;
      dvbt_hierarchy_t d_hierarchy;

      // constellation
      int d_v;
      // Bit interleaver block size
      static const int d_bsize;

      //Table to keep interleaved indices
      unsigned char * d_perm;

      // permutation function
      int H(int e, int w);

    public:
      bit_inner_interleaver_impl(int nsize, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy, dvbt_transmission_mode_t transmission);
      ~bit_inner_interleaver_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

	/*!
	* ETSI EN 300 744 Clause 4.3.4.1 \n
	* Data Input format Non-Hierarchical: \n
	* 000000X0X1 - QPSK \n
	* 0000X0X1X2X3 - QAM16 \n
	* 00X0X1X2X3X4X5 - QAM64 \n
	*
	* Data Input format Hierarchical: \n
	* 0000000X0 - H-QPSK \n
	* 0000000X1 - L-QPSK \n
	* 000000X0X1 - H-QAM16 \n
	* 000000X0X1 - L-QAM16 \n
	* 000000X0X1 - H-QAM64 \n
	* 0000X0X1X2X3 - L-QAM64 \n
	*
	* Data Output format: \n
	* 000000B0B1 - QPSK \n
	* 0000B0B1B2B3 - QAM16 \n
	* 00B0B1B2B3B4B5 - QAM64 \n
	* bit interleaver block size is 126 \n
	*/
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_BIT_INNER_INTERLEAVER_IMPL_H */

