# Copyright 2011,2012 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

########################################################################
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
add_library(gnuradio-dvbt SHARED
    test_impl.cc
    vector_pad_impl.cc
    reference_signals_impl.cc
    dvbt_config.cc
    dvbt_config.cc
    dvbt_map_impl.cc
    bit_inner_interleaver_impl.cc
    symbol_inner_interleaver_impl.cc
    inner_coder_impl.cc
    reed_solomon_enc_impl.cc
    energy_dispersal_impl.cc
    convolutional_interleaver_impl.cc
    test2_impl.cc
    demod_reference_signals_impl.cc
    dvbt_demap_impl.cc
    bit_inner_deinterleaver_impl.cc
    convolutional_deinterleaver_impl.cc
    energy_descramble_impl.cc
    reed_solomon.cc
    reed_solomon_dec_impl.cc
    ofdm_sym_acquisition_impl.cc
    fsm.cc
    base.cc
    viterbi_decoder_impl.cc )
target_link_libraries(gnuradio-dvbt ${Boost_LIBRARIES} ${GRUEL_LIBRARIES} ${GNURADIO_CORE_LIBRARIES})
set_target_properties(gnuradio-dvbt PROPERTIES DEFINE_SYMBOL "gnuradio_dvbt_EXPORTS")

########################################################################
# Install built library files
########################################################################
install(TARGETS gnuradio-dvbt
    LIBRARY DESTINATION lib${LIB_SUFFIX} # .so/.dylib file
    ARCHIVE DESTINATION lib${LIB_SUFFIX} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)

########################################################################
# Build and register unit test
########################################################################
include(GrTest)

include_directories(${CPPUNIT_INCLUDE_DIRS})
list(APPEND test_dvbt_sources
  ${CMAKE_CURRENT_SOURCE_DIR}/test_dvbt.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_dvbt.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_viterbi_decoder.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_ofdm_sym_acquisition.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_reed_solomon_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_reed_solomon_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_energy_descramble.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_convolutional_deinterleaver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_bit_inner_deinterleaver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_demod_reference_signals.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_test2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_convolutional_interleaver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_energy_dispersal.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_reed_solomon.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_inner_coder.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_symbol_inner_interleaver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_bit_inner_interleaver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_dvbt_map.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_dvbt_config.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_reference_signals.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_vector_pad.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/qa_test.cc
  )

add_executable(test-dvbt ${test_dvbt_sources})

target_link_libraries(
  test-dvbt
  ${GNURADIO_CORE_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CPPUNIT_LIBRARIES}
  gnuradio-dvbt
)

GR_ADD_TEST(test_dvbt test-dvbt)
