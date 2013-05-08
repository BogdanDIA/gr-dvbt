/* -*- c++ -*- */

#define DVBT_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "dvbt_swig_doc.i"

%{
#include "dvbt/test.h"
#include "dvbt/vector_pad.h"
#include "dvbt/reference_signals.h"
#include "dvbt/dvbt_config.h"
#include "dvbt/dvbt_map.h"
#include "dvbt/bit_inner_interleaver.h"
#include "dvbt/symbol_inner_interleaver.h"
#include "dvbt/inner_coder.h"
#include "dvbt/reed_solomon.h"
#include "dvbt/energy_dispersal.h"
#include "dvbt/convolutional_interleaver.h"
#include "dvbt/test2.h"
#include "dvbt/demod_reference_signals.h"
#include "dvbt/dvbt_demap.h"
#include "dvbt/bit_inner_deinterleaver.h"
#include "dvbt/convolutional_deinterleaver.h"
%}

%include "dvbt/test.h"
%include "dvbt/vector_pad.h"
%include "dvbt/reference_signals.h"
%include "dvbt/dvbt_config.h"
%include "dvbt/dvbt_map.h"
%include "dvbt/bit_inner_interleaver.h"
%include "dvbt/symbol_inner_interleaver.h"
%include "dvbt/inner_coder.h"
%include "dvbt/reed_solomon.h"
%include "dvbt/energy_dispersal.h"
%include "dvbt/convolutional_interleaver.h"
%include "dvbt/test2.h"
%include "dvbt/demod_reference_signals.h"
%include "dvbt/dvbt_demap.h"
%include "dvbt/bit_inner_deinterleaver.h"
%include "dvbt/convolutional_deinterleaver.h"
GR_SWIG_BLOCK_MAGIC2(dvbt, test);
GR_SWIG_BLOCK_MAGIC2(dvbt, vector_pad);
GR_SWIG_BLOCK_MAGIC2(dvbt, reference_signals);
GR_SWIG_BLOCK_MAGIC2(dvbt, dvbt_map);
GR_SWIG_BLOCK_MAGIC2(dvbt, bit_inner_interleaver);
GR_SWIG_BLOCK_MAGIC2(dvbt, symbol_inner_interleaver);
GR_SWIG_BLOCK_MAGIC2(dvbt, inner_coder);
GR_SWIG_BLOCK_MAGIC2(dvbt, reed_solomon);
GR_SWIG_BLOCK_MAGIC2(dvbt, energy_dispersal);
GR_SWIG_BLOCK_MAGIC2(dvbt, convolutional_interleaver);
GR_SWIG_BLOCK_MAGIC2(dvbt, test2);
GR_SWIG_BLOCK_MAGIC2(dvbt, demod_reference_signals);
GR_SWIG_BLOCK_MAGIC2(dvbt, dvbt_demap);
GR_SWIG_BLOCK_MAGIC2(dvbt, bit_inner_deinterleaver);
GR_SWIG_BLOCK_MAGIC2(dvbt, convolutional_deinterleaver);
