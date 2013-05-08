/*
 * Copyright 2012 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

/*
 * This class gathers together all the test cases for the gr-filter
 * directory into a single test suite.  As you create new test cases,
 * add them here.
 */

#include "qa_dvbt.h"

#include "qa_test.h"
#include "qa_vector_pad.h"
#include "qa_reference_signals.h"
#include "qa_dvbt_config.h"
#include "qa_dvbt_map.h"
#include "qa_bit_inner_interleaver.h"
#include "qa_symbol_inner_interleaver.h"
#include "qa_inner_coder.h"
#include "qa_reed_solomon.h"
#include "qa_energy_dispersal.h"
#include "qa_convolutional_interleaver.h"
#include "qa_test2.h"
#include "qa_demod_reference_signals.h"
#include "qa_bit_inner_deinterleaver.h"
#include "qa_convolutional_deinterleaver.h"

CppUnit::TestSuite *
qa_dvbt::suite()
{
  CppUnit::TestSuite *s = new CppUnit::TestSuite("dvbt");

  s->addTest(gr::dvbt::qa_test::suite());
  s->addTest(gr::dvbt::qa_vector_pad::suite());
  s->addTest(gr::dvbt::qa_reference_signals::suite());
  s->addTest(gr::dvbt::qa_dvbt_config::suite());
  s->addTest(gr::dvbt::qa_dvbt_map::suite());
  s->addTest(gr::dvbt::qa_bit_inner_interleaver::suite());
  s->addTest(gr::dvbt::qa_symbol_inner_interleaver::suite());
  s->addTest(gr::dvbt::qa_inner_coder::suite());
  s->addTest(gr::dvbt::qa_reed_solomon::suite());
  s->addTest(gr::dvbt::qa_energy_dispersal::suite());
  s->addTest(gr::dvbt::qa_convolutional_interleaver::suite());
  s->addTest(gr::dvbt::qa_test2::suite());
  s->addTest(gr::dvbt::qa_demod_reference_signals::suite());
  s->addTest(gr::dvbt::qa_bit_inner_deinterleaver::suite());
  s->addTest(gr::dvbt::qa_convolutional_deinterleaver::suite());

  return s;
}
