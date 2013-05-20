# Copyright 2011 Free Software Foundation, Inc.
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
# Include python install macros
########################################################################
include(GrPython)
if(NOT PYTHONINTERP_FOUND)
    return()
endif()

########################################################################
# Install python sources
########################################################################
GR_PYTHON_INSTALL(
    FILES
    __init__.py
    ofdm_receiver.py
    #ofdm.py
    ofdm_sync_pn.py
    ofdm_sync_pnac.py
    ofdm_sync_ml.py
    ofdm_sync_fixed.py
    DESTINATION ${GR_PYTHON_DIR}/dvbt
    COMPONENT "dvbt_python"
)

########################################################################
# Handle the unit tests
########################################################################
include(GrTest)

set(GR_TEST_TARGET_DEPS gnuradio-dvbt)
set(GR_TEST_PYTHON_DIRS ${CMAKE_BINARY_DIR}/swig)
GR_ADD_TEST(qa_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_test.py)
GR_ADD_TEST(qa_vector_pad ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_vector_pad.py)
GR_ADD_TEST(qa_reference_signals ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_reference_signals.py)
GR_ADD_TEST(qa_dvbt_map ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_dvbt_map.py)
GR_ADD_TEST(qa_bit_inner_interleaver ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_bit_inner_interleaver.py)
GR_ADD_TEST(qa_symbol_inner_interleaver ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_symbol_inner_interleaver.py)
GR_ADD_TEST(qa_inner_coder ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_inner_coder.py)
GR_ADD_TEST(qa_reed_solomon ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_reed_solomon.py)
GR_ADD_TEST(qa_reed_solomon ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_reed_solomon_enc.py)
GR_ADD_TEST(qa_energy_dispersal ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_energy_dispersal.py)
GR_ADD_TEST(qa_convolutional_interleaver ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_convolutional_interleaver.py)
GR_ADD_TEST(qa_test2 ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_test2.py)
GR_ADD_TEST(qa_demod_reference_signals ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_demod_reference_signals.py)
GR_ADD_TEST(qa_dvbt_demap ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_dvbt_demap.py)
GR_ADD_TEST(qa_bit_inner_deinterleaver ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_bit_inner_deinterleaver.py)
GR_ADD_TEST(qa_convolutional_deinterleaver ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_convolutional_deinterleaver.py)
GR_ADD_TEST(qa_energy_descramble ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_energy_descramble.py)
GR_ADD_TEST(qa_reed_solomon_dec ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_reed_solomon_dec.py)
GR_ADD_TEST(qa_ofdm_sym_acquisition ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/qa_ofdm_sym_acquisition.py)
