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
install(FILES
    dvbt_test.xml
    dvbt_vector_pad.xml
    dvbt_reference_signals.xml
    dvbt_dvbt_map.xml
    dvbt_bit_inner_interleaver.xml
    dvbt_symbol_inner_interleaver.xml
    dvbt_inner_coder.xml
    dvbt_reed_solomon_enc.xml
    dvbt_energy_dispersal.xml
    dvbt_convolutional_interleaver.xml
    dvbt_ofdm_demod.xml
    dvbt_test2.xml
    dvbt_demod_reference_signals.xml
    dvbt_dvbt_demap.xml
    dvbt_bit_inner_deinterlever.xml
    dvbt_convolutional_deinterleaver.xml
    dvbt_energy_descramble.xml
    dvbt_reed_solomon_dec.xml
    dvbt_ofdm_sym_acquisition.xml DESTINATION share/gnuradio/grc/blocks
)
