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
#include "viterbi_decoder_impl.h"

static const float INF = 1.0e9;

namespace gr {
  namespace dvbt {

  template <class T>
  void viterbi_algorithm(int I, int S, int O,
              const std::vector<int> &NS,
              const std::vector<int> &OS,
              const std::vector< std::vector<int> > &PS,
              const std::vector< std::vector<int> > &PI,
              int K,
              int S0,int SK,
              const float *in, T *out)//,
              //std::vector<int> &trace)
  {
    std::vector<int> trace(S*K);
    std::vector<float> alpha(S*2);
    int alphai;
    float norm,mm,minm;
    int minmi;
    int st;

#if 0
    printf("I: %i\n", I);
    printf("S: %i\n", S);
    printf("O: %i\n", O);

    for (int i = 0; i < NS.size(); i++)
      printf("NS[%i]: %i\n", i, NS[i]);

    for (int i = 0; i < OS.size(); i++)
      printf("OS[%i]: %i\n", i, OS[i]);


    for (int i = 0; i < PS.size(); i++)
      for (int j = 0; j < PS[i].size(); j++)
        printf("PS[%i][%i]: %i\n", i, j, (PS[i])[j]);

    for (int i = 0; i < PI.size(); i++)
      for (int j = 0; j < PI[i].size(); j++)
        printf("PI[%i][%i]: %i\n", i, j, (PI[i])[j]);

    printf("K: %i\n", K);
    printf("S0: %i\n", S0);
    printf("SK: %i\n", SK);
#endif

    if(S0<0) { // initial state not specified
        for(int i=0;i<S;i++) alpha[0*S+i]=0;
    }
    else {
        for(int i=0;i<S;i++) alpha[0*S+i]=INF;
        alpha[0*S+S0]=0.0;
    }

    alphai=0;
    for(int k=0;k<K;k++) {
        norm=INF;
        for(int j=0;j<S;j++) { // for each next state do ACS
            minm=INF;
            minmi=0;
            for(unsigned int i=0;i<PS[j].size();i++) {
                //int i0 = j*I+i;
                //if((mm=alpha[alphai*S+PS[j][i]]+in[k*O+OS[PS[j][i]*I+PI[j][i]]])<minm)
                    //minm=mm,minmi=i;

          // Calculate Hamming distance between received codeword
          // and output symbol
          int metric = OS[PS[j][i]*I+PI[j][i]] ^ (int)in[k];
          int distance = 0;
          while (metric)
          {
      if (metric & 1)
        distance++;

      metric >>= 1;
          }

          if((mm=(int)alpha[alphai*S+PS[j][i]]+distance)<minm)
                    minm=mm,minmi=i;

#if 0
          printf("in[%i]: %f\n", k*O+OS[PS[j][i]*I+PI[j][i]], in[k*O+OS[PS[j][i]*I+PI[j][i]]]);
          printf("PS[%i][%i]: %i, PI[%i][%i]: %i\n", j, i, PS[j][i], j, i, PI[j][i]);
          printf("k: %i, S: %i, j: %i\n", k, S, j);
          printf("alpha[%i]: %f, mm: %f\n", alphai*S+PS[j][i], alpha[alphai*S+PS[j][i]], mm);
#endif
            }
            trace[k*S+j]=minmi;
            alpha[((alphai+1)%2)*S+j]=minm;
#if 0
      printf("******* minm: %f, minmi: %i\n", minm, minmi);
      printf("******* trace[%i]: %i\n", k*S+j, trace[k*S+j]);
      printf("******* alpha[%i]: %f\n", ((alphai+1)%2)*S+j, minm);
#endif
            if(minm<norm) norm=minm;
        }
        for(int j=0;j<S;j++)
            alpha[((alphai+1)%2)*S+j]-=norm; // normalize total metrics so they do not explode
        alphai=(alphai+1)%2;
    }

    if(SK<0) { // final state not specified
        minm=INF;
        minmi=0;
        for(int i=0;i<S;i++)
            if((mm=alpha[alphai*S+i])<minm) minm=mm,minmi=i;
        st=minmi;
    }
    else {
        st=SK;
    }

    for(int k=K-1;k>=0;k--) { // traceback
        int i0=trace[k*S+st];
        out[k]= (T) PI[st][i0];
#if 0
        printf("st: %i\n", st);
        printf("trace[%i]: %i\n", k*S+st, trace[k*S+st]);
        printf("out[%i]: %i, PI[%i][%i]: %i\n", k, out[k], st, i0, out[k]);
#endif
        st=PS[st][i0];
    }

  }


    viterbi_decoder::sptr
    viterbi_decoder::make(const fsm &FSM, int K, int S0, int SK)
    {
      return gnuradio::get_initial_sptr (new viterbi_decoder_impl(FSM, K, S0, SK));
    }

    /*
     * The private constructor
     */
    viterbi_decoder_impl::viterbi_decoder_impl(const fsm &FSM, int K, int S0, int SK)
      : gr_block("viterbi_decoder",
		      gr_make_io_signature(1, 1, sizeof (float)),
		      gr_make_io_signature(1, 1, sizeof (unsigned char))),
      d_FSM (FSM),
      d_K (K),
      d_S0 (S0),
      d_SK (SK)
    {
      set_relative_rate (1.0 / ((double) d_FSM.O()));
      set_output_multiple (d_K);
    }

    /*
     * Our virtual destructor.
     */
    viterbi_decoder_impl::~viterbi_decoder_impl()
    {
    }

    void
    viterbi_decoder_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
       assert (noutput_items % d_K == 0);
       int input_required =  d_FSM.O() * noutput_items ;
       unsigned ninputs = ninput_items_required.size();
       for (unsigned int i = 0; i < ninputs; i++) {
         ninput_items_required[i] = input_required;
       }
    }

    int
    viterbi_decoder_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
        assert (input_items.size() == output_items.size());
        int nstreams = input_items.size();
        assert (noutput_items % d_K == 0);
        int nblocks = noutput_items / d_K;

        for (int m=0;m<nstreams;m++) {
          const float *in = (const float *) input_items[m];
          unsigned char *out = (unsigned char *) output_items[m];
          for (int n=0;n<nblocks;n++) {
            viterbi_algorithm(d_FSM.I(),d_FSM.S(),d_FSM.O(),d_FSM.NS(),d_FSM.OS(),d_FSM.PS(),d_FSM.PI(),d_K,d_S0,d_SK,&(in[n*d_K]),&(out[n*d_K]));
          }
        }

        // Tell runtime system how many input items we consumed on
        // each input stream.
        consume_each (noutput_items);

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }

  } /* namespace dvbt */
} /* namespace gr */

