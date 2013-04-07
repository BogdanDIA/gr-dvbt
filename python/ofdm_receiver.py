#!/usr/bin/env python
#
# Copyright 2006, 2007, 2008 Free Software Foundation, Inc.
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
# 

from gnuradio import blks2
import math
from numpy import fft
from gnuradio import gr

from gnuradio.digital import digital_swig
from dvbt import ofdm_sync_pn
from dvbt import ofdm_sync_fixed
from dvbt import ofdm_sync_pnac
from dvbt import ofdm_sync_ml

class ofdm_receiver(gr.hier_block2):
    """
    Performs receiver synchronization on OFDM symbols.

    The receiver performs channel filtering as well as symbol, frequency, and phase synchronization.
    The synchronization routines are available in three flavors: preamble correlator (Schmidl and Cox),
    modifid preamble correlator with autocorrelation (not yet working), and cyclic prefix correlator
    (Van de Beeks).
    """

    def __init__(self, fft_length, cp_length, occupied_tones, snr, logging=False):
        """
	Hierarchical block for receiving OFDM symbols.

	The input is the complex modulated signal at baseband.
        Synchronized packets are sent back to the demodulator.

        @param fft_length: total number of subcarriers
        @type  fft_length: int
        @param cp_length: length of cyclic prefix as specified in subcarriers (<= fft_length)
        @type  cp_length: int
        @param occupied_tones: number of subcarriers used for data
        @type  occupied_tones: int
        @param snr: estimated signal to noise ratio used to guide cyclic prefix synchronizer
        @type  snr: float
        @param logging: turn file logging on or off
        @type  logging: bool
	"""

	gr.hier_block2.__init__(self, "ofdm_receiver",
				gr.io_signature(1, 1, gr.sizeof_gr_complex), # Input signature
                                gr.io_signature2(2, 2, gr.sizeof_gr_complex*fft_length, gr.sizeof_char*fft_length)) # Output signature
        
        bw = (float(occupied_tones) / float(fft_length)) / 2.0
        tb = bw*0.08
        chan_coeffs = gr.firdes.low_pass (1.0,                     # gain
                                          1.0,                     # sampling rate
                                          bw+tb,                   # midpoint of trans. band
                                          tb,                      # width of trans. band
                                          gr.firdes.WIN_HAMMING)   # filter type
        #self.chan_filt = gr.fft_filter_ccc(1, chan_coeffs)
        self.chan_filt = gr.multiply_cc(1)
        
        win = [1 for i in range(fft_length)]

        ks0time = fft_length*(0,)

        SYNC = "ml"
        if SYNC == "ml":
            #TODO -1.0/...
            nco_sensitivity = -1.0/fft_length   # correct for fine frequency
            self.ofdm_sync = ofdm_sync_ml.ofdm_sync_ml(fft_length,
                                          cp_length,
                                          snr,
                                          ks0time,
                                          logging)
        elif SYNC == "pn":
            nco_sensitivity = -2.0/fft_length   # correct for fine frequency
            self.ofdm_sync = ofdm_sync_pn.ofdm_sync_pn(fft_length,
                                          cp_length,
                                          logging)
        elif SYNC == "pnac":
            nco_sensitivity = -2.0/fft_length   # correct for fine frequency
            self.ofdm_sync = ofdm_sync_pnac.ofdm_sync_pnac(fft_length,
                                            cp_length,
                                            ks0time,
                                            logging)
        # for testing only; do not user over the air
        # remove filter and filter delay for this
        elif SYNC == "fixed":
            self.chan_filt = gr.multiply_const_cc(1.0) 
            nsymbols = 18      # enter the number of symbols per packet
            freq_offset = 0.0  # if you use a frequency offset, enter it here
            nco_sensitivity = -2.0/fft_length   # correct for fine frequency
            self.ofdm_sync = ofdm_sync_fixed.ofdm_sync_fixed(fft_length,
                                             cp_length,
                                             nsymbols,
                                             freq_offset,
                                             logging)
        # Set up blocks

        self.nco = gr.frequency_modulator_fc(nco_sensitivity)         # generate a signal proportional to frequency error of sync block
        self.sigmix = gr.multiply_cc()
        self.sampler = digital_swig.ofdm_sampler(fft_length, fft_length+cp_length)
        self.fft_demod = gr.fft_vcc(fft_length, True, win, True)

        self.connect(self, self.chan_filt)                            # filter the input channel
        self.connect(self.chan_filt, self.ofdm_sync)                  # into the synchronization alg.
        self.connect((self.ofdm_sync,0), self.nco, (self.sigmix,1))   # use sync freq. offset output to derotate input signal
        self.connect(self.chan_filt, (self.sigmix,0))                 # signal to be derotated

        self.connect(self.sigmix, (self.sampler,0))                   # sample off timing signal detected in sync alg
        #self.connect(self.chan_filt, (self.sampler,0))                   # sample off timing signal detected in sync alg

        self.connect((self.ofdm_sync,1), (self.sampler,1))            # timing signal to sample at
        self.connect((self.sampler,0), self.fft_demod)                # send derotated sampled signal to FFT
        self.connect(self.fft_demod, (self,0))                        # frequency domain signal sent to output 0
        self.connect((self.sampler,1), (self,1))                      # timing sent to output 1
        print "setup OK"

        logging = 0

        if logging:
            self.connect(self.chan_filt, gr.file_sink(gr.sizeof_gr_complex, "ofdm_receiver-chan_filt_c.dat"))
            self.connect((self.ofdm_sync, 0), gr.file_sink(gr.sizeof_float, "ofdm_receiver-sync_0.dat"))
            self.connect((self.ofdm_sync, 1), gr.file_sink(gr.sizeof_char, "ofdm_receiver-sync_1.dat"))
            self.connect(self.nco, gr.file_sink(gr.sizeof_gr_complex, "ofdm_receiver-nco_c.dat"))
            self.connect(self.sigmix, gr.file_sink(gr.sizeof_gr_complex, "ofdm_receiver-sigmix_c.dat"))
            self.connect((self.sampler, 0), gr.file_sink(gr.sizeof_gr_complex*fft_length, "ofdm_receiver-sampler_0.dat"))
            self.connect((self.sampler, 1), gr.file_sink(gr.sizeof_char*fft_length, "ofdm_receiver-sampler_1.dat"))
            self.connect(self.fft_demod, gr.file_sink(gr.sizeof_gr_complex*fft_length, "ofdm_receiver-fft_out_c.dat"))
