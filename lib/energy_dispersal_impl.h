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

#ifndef INCLUDED_DVBT_ENERGY_DISPERSAL_IMPL_H
#define INCLUDED_DVBT_ENERGY_DISPERSAL_IMPL_H

#include <dvbt/energy_dispersal.h>

namespace gr {
  namespace dvbt {
    /*! 
     * \brief DVBT Energy dispersal class. \n 
     * \ingroup dvbt
     * Randomizes MPEG-2 packets using a PRBS generator. \n
     * Each block has 188 bytes. \n
     * \param nsize	number of blocks \n
     */
    class energy_dispersal_impl : public energy_dispersal
    {
    private:
      static const int d_nblocks;
      static const int d_bsize;
      static const int d_sync;
      static const int d_nsync;

      int d_reg;

      void init_prbs();	
      int clock_prbs(int clocks);

    public:
      energy_dispersal_impl(int nsize);
      ~energy_dispersal_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      /*!
       * ETSI EN 300 744 - Clause 4.3.1. \n
       * Input - MPEG-2 transport packets (including sync - 0x47). \n
       * Output - Randomized MPEG-2 transport packets. \n
       * We assume the first byte is a sync.\n
       * First sync in a row of 8 packets is reversed - 0xB8. \n
       * Block size is 188bytes. \n
       */
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_ENERGY_DISPERSAL_IMPL_H */

