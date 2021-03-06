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


#ifndef INCLUDED_DVBT_REED_SOLOMON_ENC_H
#define INCLUDED_DVBT_REED_SOLOMON_ENC_H

#include <dvbt/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace dvbt {

    /*!
     * \brief Reed Solomon encoder
     * \ingroup dvbt
     */
    class DVBT_API reed_solomon_enc : virtual public block
    {
    public:
       typedef boost::shared_ptr<reed_solomon_enc> sptr;

       /*!
        * \brief Return a shared_ptr to a new instance of dvbt::reed_solomon.
        *
        * To avoid accidental use of raw pointers, dvbt::reed_solomon's
        * constructor is in a private implementation
        * class. dvbt::reed_solomon::make is the public interface for
        * creating new instances.
        */
       static sptr make(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_REED_SOLOMON_ENC_H */

