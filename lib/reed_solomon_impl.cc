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
#include "reed_solomon_impl.h"
#include <stdio.h>

namespace gr {
  namespace dvbt {

    void
    reed_solomon_impl::gf_init(int p, int m, int gfpoly)
    {
      d_p = p; d_m = m;

      //maximum number of elements in the GF(p^m)
      int q = powl(p, m);

      d_gf_exp = new unsigned char[q];
      if (d_gf_exp == NULL)
      {
        std::cout << "Cannot allocate memory" << std::endl;
        return;
      }

      d_gf_log = new unsigned char[q];
      if (d_gf_log == NULL)
      {
        std::cout << "Cannot allocate memory" << std::endl;
        delete [] d_gf_exp;
        return;
      }

      int reg_rs = 1;

      d_gf_exp[q - 1] = 0;
      d_gf_log[0] = q - 1;

      for (int i = 0; i < (q - 1); i++)
      {
        d_gf_exp[i] = reg_rs;
        d_gf_log[reg_rs] = i;

        //This is equvalent with raise to power
        reg_rs = reg_rs << 1;

        if (reg_rs & (1 << m))
          reg_rs = reg_rs ^ gfpoly;

        reg_rs = reg_rs & ((1 << m) - 1);
      }
    }

    void
    reed_solomon_impl::gf_uninit()
    {
      delete [] d_gf_exp;
      delete [] d_gf_log;
    }

    int
    reed_solomon_impl::gf_exp(int a)
    {
      return d_gf_exp[a % ((1 << d_m) - 1)];
    }

    int
    reed_solomon_impl::gf_add(int a, int b)
    {
      return (a ^ b);
    }

    int
    reed_solomon_impl::gf_mul(int a, int b)
    {
      if (a == 0 || b == 0)
        return 0;
      else
        return gf_exp(d_gf_log[a] + d_gf_log[b]);
    }

    void
    reed_solomon_impl::rs_init(int lambda, int n, int k, int t)
    {
      d_n = n; d_k = k; d_t = t;

      d_l = new unsigned char[2 * d_t];
      if (d_l == NULL)
        return;

      d_g = new unsigned char[2 * d_t + 1];
      if (d_g == NULL)
      {
        delete [] d_l;
        return;
      }

      //Generate roots of lambda
      d_l[0] = 1;

      for (int i = 1; i < (2 * t); i++)
        d_l[i] = gf_mul(d_l[i - 1], lambda);

      //Init Generator polynomial buffer
      for (int i = 0; i <= (2*t); i++)
        d_g[i] = 0;

      //Start with x+lambda^0
      d_g[0] = 1;

      //Create generator polynomial
      for (int i = 1; i <= (2 * t); i++)
      {
        for (int j = i; j > 0; j--)
        {
          if (d_g[j] != 0)
            d_g[j] = gf_add(d_g[j - 1], gf_mul(d_g[j], d_l[i - 1]));
          else
            d_g[j] = d_g[j - 1];
        }

        d_g[0] = gf_mul(d_g[0], d_l[i - 1]);
      }
    }

    void
    reed_solomon_impl::rs_uninit()
    {
      if (d_l)
        delete [] d_l;
      if (d_g)
        delete [] d_g;
    }

    int
    reed_solomon_impl::rs_encode(const unsigned char *data_in, unsigned char *parity)
    {
      memset(parity, 0, 2 * d_t);

      for (int i = 0; i < d_k; i++)
      {
        int feedback = gf_add(data_in[i], parity[0]); 

        if (feedback != 0)
        {
          for (int j = 1; j < (2 * d_t); j++)
          {
            if (d_g[2 * d_t - j] != 0)
              parity[j] = gf_add(parity[j], gf_mul(feedback, d_g[2 * d_t - j]));
          }
        }

        //Shift the register
        memmove(&parity[0], &parity[1], (2 * d_t) - 1);

        if (feedback != 0)
          parity[2 * d_t - 1] = gf_mul(feedback, d_g[0]);
        else
          parity[2 * d_t - 1] = 0;
      }
    }

    reed_solomon::sptr
    reed_solomon::make(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks)
    {
      return gnuradio::get_initial_sptr (new reed_solomon_impl(p, m, gfpoly, n, k, t, s, blocks));
    }

    /*
     * The private constructor
     */
    reed_solomon_impl::reed_solomon_impl(int p, int m, int gfpoly, int n, int k, int t, int s, int blocks)
      : gr_block("reed_solomon",
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * blocks * (k - s)),
		      gr_make_io_signature(1, 1, sizeof(unsigned char) * blocks * (n - s))),
      d_p(p), d_m(m), d_gfpoly(gfpoly), d_n(n), d_k(k), d_t(t), d_s(s), d_blocks(blocks)
    {
      gf_init(d_p, d_m, d_gfpoly);
      rs_init(d_p, d_n, d_k, d_t);

      //Allocate buffer for RS input
      d_in = new unsigned char[d_k];
      if (d_in == NULL)
      {
        std::cout << "Cannot allocate memory" << std::endl;
      }
      //For shortened code, first s bytes are zero
      memset(&d_in[0], 0, d_s);
    }

    /*
     * Our virtual destructor.
     */
    reed_solomon_impl::~reed_solomon_impl()
    {
      if (d_in)
        delete [] d_in;

      rs_uninit();
      gf_uninit();
    }

    void
    reed_solomon_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
        ninput_items_required[0] = noutput_items;
    }

    int
    reed_solomon_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items[0];
        unsigned char *out = (unsigned char *) output_items[0];
 
        int in_bsize = d_k - d_s;
        int out_bsize = d_n - d_s;

        unsigned char parity[2 * d_t];

        int in_count = 0;
        int out_count = 0;

        for (int i = 0; i < (d_blocks * noutput_items); i++)
        {
          memcpy(&d_in[d_s], &in[i * in_bsize], in_bsize);

          //TODO - zero copy?
          rs_encode(d_in, parity);

          memcpy(&out[i * out_bsize], &in[i * in_bsize], in_bsize);
          memcpy(&out[i * out_bsize + in_bsize], parity, 2 * d_t);
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

