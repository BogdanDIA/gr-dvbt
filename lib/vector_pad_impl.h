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

#ifndef INCLUDED_DVBT_VECTOR_PAD_IMPL_H
#define INCLUDED_DVBT_VECTOR_PAD_IMPL_H

#include <dvbt/vector_pad.h>

namespace gr {
  namespace dvbt {
    /*!
     * \brief Pad a vector for IFFT. \n
     * \ingroup dvbt
     * \param itemsize size of in/out item. \n
     * \param ninput length of input stream. \n
     * \param noutput length of output stream. \n
     */
    class vector_pad_impl : public vector_pad
    {
    private:
      int d_itemsize;
      int d_ninput;
      int d_noutput;

      int d_prefix_len;
      int d_ninput_len;
      int d_suffix_len;

    public:
      vector_pad_impl(int itemsize, int ninput, int noutput);
      ~vector_pad_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      /*!
       * This version adds a prefix and a suffix to the input data. \n
       * Padding is with zero value.
       */
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_VECTOR_PAD_IMPL_H */

