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

#ifndef INCLUDED_DVBT_CONVOLUTIONAL_INTERLEAVER_IMPL_H
#define INCLUDED_DVBT_CONVOLUTIONAL_INTERLEAVER_IMPL_H

#include <dvbt/convolutional_interleaver.h>
#include <vector>
#include <deque>

namespace gr {
  namespace dvbt {
    /*!
     * \brief Convolutional interleaver
     * \ingroup dvbt
     * \param blocks number of blocks to process \n
     * \param I size of a block \n
     * \param M depth length for each element in shift registers \n
     */

    class convolutional_interleaver_impl : public convolutional_interleaver
    {
    private:
      int d_blocks;
      int d_I;
      int d_M;
      std::vector< std::deque<unsigned char> * > d_shift;

    public:
      convolutional_interleaver_impl(int nsize, int I, int M);
      ~convolutional_interleaver_impl();

      /*!
       * ETSI EN 300 744 Clause 4.3.1. \n
       * Forney (Ramsey type III) convolutional interleaver. \n
       * Data input: Blocks of I bytes size. \n
       * Data output: Stream of 1 byte elements. \n
       */
      int work(int noutput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_CONVOLUTIONAL_INTERLEAVER_IMPL_H */

