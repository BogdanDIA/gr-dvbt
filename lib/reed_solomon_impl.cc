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

#define min(a,b) ((a) < (b)) ? (a) : (b)

#define DEBUG

#ifdef DEBUG
#define PRINTF printf
#else
#define PRINTF
#endif

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
      return d_gf_exp[a % d_n];
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

    int
    reed_solomon_impl::gf_div(int a, int b)
    {
      if (a == 0 || b == 0)
        return (0);

      return (gf_exp(d_n + d_gf_log[a] - d_gf_log[b]));
    }

    int
    reed_solomon_impl::gf_pow(int a, int power)
    {
      if (a == 0)
        return (0);

      return gf_exp(d_n + d_gf_log[a] + power);
    }

    int
    reed_solomon_impl::gf_lpow(int power)
    {
      return d_l[power % d_n];
    }

    void
    reed_solomon_impl::rs_init(int lambda, int n, int k, int t)
    {
      d_lambda= lambda, d_n = n; d_k = k; d_t = t;
      // 2t = n - k, dmin = 2t + 1 = n -k + 1

      d_l = new unsigned char[d_n + 1]; // 2* d_t
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

      for (int i = 1; i <= d_n; i++) // TODO - 2 *d_t
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

      // Init syndrome array 
      d_syn = new unsigned char[2 * d_t];
      if (d_syn == NULL)
      {
        std::cout << "Cannot allocate memory" << std::endl;
        delete [] d_l;
        delete [] d_g;
        return;
      }
    }

    void
    reed_solomon_impl::rs_uninit()
    {
      if (d_l)
        delete [] d_l;
      if (d_g)
        delete [] d_g;
      if (d_syn)
        delete [] d_syn;
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

      return (0);
    }

    int
    reed_solomon_impl::rs_decode(unsigned char *data, unsigned char *eras, const int no_eras)
    {
      unsigned char *tau = new unsigned char[2 * d_t + 1];
      unsigned char *sigma = new unsigned char[2 * d_t + 1];
      unsigned char *b = new unsigned char [2 * d_t + 1];
      unsigned char *T = new unsigned char[2 * d_t + 1];
      unsigned char *reg = new unsigned char[2 * d_t + 1];
      unsigned char *root = new unsigned char[2 * d_t + 1];
      unsigned char *loc = new unsigned char[2 * d_t + 1];
      unsigned char *omega = new unsigned char[2 * d_t];

      // Compute erasure locator polynomial
      memset(sigma, 0, 2 * d_t + 1);
      sigma[0] = 1;

      if (no_eras > 0)
      {
        // In this case we know the locations of errors
        // Init sigma to be the erasure locator polynomial
        sigma[1] = gf_exp(d_n-1-eras[0]);

        for (int i = 1; i < no_eras; i++)
        {
          int u = d_n-1-eras[i];

          for (int j = i+1; j > 0; j--) 
            sigma[j] = gf_add(sigma[j], gf_pow(sigma[j - 1], u));
        }

        for (int i = 0; i <= no_eras; i++)
          PRINTF("sigma eras[%i]: %i\n", i, sigma[i]);
      }

      // Calculate syndrome
      int syn_error = 0;

      for (int j = 0; j < 2 * d_t; j++)
        d_syn[j] = data[0];

      for (int j = 1; j < d_n; j++)
      {
        for (int i = 0; i < 2 * d_t; i++)
          d_syn[i] = gf_add(data[j], gf_pow(d_syn[i], i));

        syn_error |= d_syn[j];
      }

      for (int i = 0; i < 2 * d_t; i++)
        PRINTF("S[%i]: %i\n", i, d_syn[i]);

      if (!syn_error)
      {
        // The syndrome is a codeword
        // Return data unmodified
        PRINTF("data is codeword\n");
        return (0);
      }

      // Use Modified (errors+erasures) BMA. Algorithm of Berlekamp-Massey
      // S(i)=r(lambda^i)=e(lambda^i)

      int r = no_eras;
      int el = no_eras;
      unsigned char d_discr;

      memcpy(b, sigma, 2 * d_t + 1);

      while (++r <= 2 * d_t)
      {
        int d_discr = 0;

        for (int i = 0; i < r; i++)
          d_discr = gf_add(d_discr, gf_mul(sigma[i], d_syn[r - i - 1]));

        PRINTF("r: %i, discr: %i\n", r, d_discr);

        if (d_discr == 0)
        {
          // b(x) = x * b(x)
          memmove(&b[1], b, 2 * d_t);
          b[0] = 0;
        }
        else
        {
          T[0] = sigma[0];

          // T(x) = sigma(x) + d*x*b(x)
          for (int i = 0; i < 2 * d_t; i++)
            T[i + 1] = gf_add(sigma[i + 1], gf_mul(d_discr, b[i]));

          if (2 * el <= r + no_eras - 1)
          {
            el = r + no_eras - el;

            // b(i) = sigma(i) / discr
            for (int i = 0; i <= 2 * d_t; i++)
              b[i] = gf_div(sigma[i], d_discr);
          }
          else
          {
            // b(x) = x*b(x)
            memmove(&b[1], b, 2 * d_t);
            b[0] = 0;
          }
          memcpy(sigma, T, 2 * d_t + 1);
        }
      }

      // Compute degree(sigma)
      int deg_sigma = 0;

      for (int i = 0; i < 2 * d_t + 1; i++)
      {
        PRINTF("sigma[%i]: %i\n", i, sigma[i]);
        if (sigma[i] != 0)
          deg_sigma = i;
      }

      PRINTF("deg_sigma: %i\n", deg_sigma);

      // Find the roots of sigma(x) by Chien search
      // Test sum(1)=1+sigma(1)*(lambda^1)+...+sigma(nu)*lambda(^nu)
      // Test sum(2)=1+sigma(1)*(lambda^2)+...+sigma(nu)*lambda(^nu*2)
      // ...
      // Test sum(l)=1+sigma(1)*(lambda^l)+...+sigma(nu)*lambda(^nu*l)
      // in order to see if lambda^(-1) is a root
      // where nu is degree(sigma)

      int no_roots = 0;

      memcpy(&reg[1], &sigma[1], 2 * d_t);

      for (int i = 1; i <= d_n; i++)
      {
        int q = 1;

        for (int j = deg_sigma; j > 0; j--)
        {
          reg[j] = gf_pow(reg[j], j);
          q = gf_add(q, reg[j]);
        }

        if (q != 0)
          continue;

        // We are here when we found roots of the sigma(x)
        // Keep roots in index form
        root[no_roots] = i;
        loc[no_roots] = i - 1;

        PRINTF("root[%i]: %i\n", no_roots, root[no_roots]);
        PRINTF("loc[%i]: %i\n", no_roots, loc[no_roots]);

        if (++no_roots == deg_sigma)
          break;
      }

      if (no_roots != deg_sigma)
      {
        // Uncorrectable error detected
        if (eras)
        {
          for (int i = 0; i < no_roots; i++)
            eras[i] = loc[i];
        }

        return (-1);
      }

      // Compute erros+erasures evaluator polynomial
      // omega(x)=sigma(x)S(x) mod (x ^ 2 * t)
      int deg_omega = 0;

      for (int i = 0; i < 2 * d_t; i++)
      {
        int tmp = 0;
        int j = (deg_sigma < i) ? deg_sigma : i;

        for(;j >= 0; j--)
          tmp = gf_add(tmp, gf_mul(d_syn[i - j], sigma[j]));

        if(tmp != 0)
          deg_omega = i;

        omega[i] = tmp;
      }
      omega[2 * d_t] = 0;

      PRINTF("deg_omega: %i\n", deg_omega);

      for (int i = 0; i < 2 * d_t; i++)
        PRINTF("omega[%i]: %i\n", i, omega[i]);

      // Compute error values using Forney formula (poly form)
      // e(j(l))) = (lambda(j(l)) ^ 2) * omega(lambda ^ (-j(l))) / sigma_pr(lambda ^ (-j(l)))
      // where sigma_pr is the formal derivative of sigma

      for (int j = no_roots - 1; j >= 0; j--)
      {
        int num1 = 0;

        // roots[] are in index form
        for (int i = deg_omega; i >= 0; i--)
          num1 = gf_add(num1, gf_pow(omega[i], i * root[j]));

        // root[] is in index form
        int num2 = gf_exp(root[j] * (-1) + d_n);

        int den = 0;

        // sigma[i+1] for i even is the formal derivative lambda_pr of sigma[i]
        int deg_max = min(deg_sigma, 2 * d_t - 1);

        for (int i = 1; i <= deg_max; i += 2)
          den = gf_add(den, gf_exp(d_gf_log[sigma[i]] + (i - 1) * root[j]));

        PRINTF("num1: %i, num2: %i\n", num1, num2);
        PRINTF("den: %i: %i\n", j, den);

        if (den == 0)
        {
          PRINTF("RS: denominator is null\n");
          if (eras)
          {
            for (int i = 0; i < no_roots; i++)
              eras[i] = loc[i];
          }
          return (-1);
        }

        int err = gf_div(gf_mul(num1, num2), den);

        data[loc[j]] = gf_add(data[loc[j]], err);

        PRINTF("data[%i]: %i\n", loc[j], data[loc[j]]);
      }

      return(no_roots);
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

      /************************************************/
      PRINTF("RS begin\n");

      unsigned char * datak = new unsigned char[d_k];
      unsigned char * datan = new unsigned char[d_n];
      unsigned char * datap = new unsigned char[d_n - d_k];

      // Init datak
      for (int i = 0; i < d_k; i++)
      {
        datak[i] = i;
        PRINTF("encode data_in[%i]: %i\n", i, datak[i]);
      }

      rs_encode(datak, datap);

      for (int i = 0; i < (d_n - d_k); i++)
        PRINTF("parity[%i]: %i\n", i, datap[i]);

      // Format data for decoder
      memcpy(&datan[0], &datak[0], d_k);
      memcpy(&datan[d_k], &datap[0], d_n - d_k);


      // Distort the data
      
      unsigned char no_eras = 6;
      unsigned char eras[no_eras];

      for (int i = 0; i < no_eras; i++)
      {
        datan[2 * i] = 0;
        eras[i] = 2 * i;
      }

      int err = no_eras;

      for (int i = 0; i < 5; i++)
        datan[3 * (err++)] = 0;

      rs_decode(datan, &eras[0], no_eras);
      //rs_decode(datan, NULL, 0);

      for (int i = 0; i < d_k; i++)
        PRINTF("decode data_out[%i]: %i\n", i, datan[i]);
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
          //TODO - zero copy?
          memcpy(&d_in[d_s], &in[i * in_bsize], in_bsize);

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

