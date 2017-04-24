### PROJECT: DVB-T implementation using gnuradio

**_Late note: As of 2015, I donated the gr-dvbt project to gnuradio. It is now integrated in the mainiline of gnuradio._**  

**This is a opensource implementation of DVB-T encoder/decoder according to ETSI 300 744**  
Works with gnuradio 3.7 and more details can be seen here:  
http://yo3iiu.ro/blog/?p=1191  
http://yo3iiu.ro/blog/?p=1220  
http://yo3iiu.ro/blog/?p=1244  


### HowTo:
**Run TX**   
The dvbt TX implementation now supports 2k/8k mode, QPSK, QAM16, QAM64 constelaltions and rates 1/2, 2/3, 3/4, 5/6, 7/8.  
The simplest way to run DVB-T encoding is to use the gnuradio-companion flowgraphs apps/dvbt_tx_demo.grc. It will start with a MPEG-2 TS file and will eventually generate the 10Msps baseband samples.  
Open dvbt_tx_demo.grc and run it for generating testBB.bin from test.ts. This specific flowgraph has the parameters set as: 2k OFDM, FEC code 1/2, Modulation 16-QAM, Guard Interval 1/32  

**Run RX** 
The dvbt RX implementation now supports all constellation types QPSK, QAM16, QAM64 and rates: 1/2, 2/3, 3/4, 5/6, 7/8.  
All combinations of rates and constellations are tested and working for 2k OFDM.  
To run DVB-T decoding just run apps/dvbt_rx_demo.grc. It will take the baseband samples and turn them into MPEG-2 TS file.  
testBB.bin(BB)->test_out.ts(MPEG-2 TS)  

Then the output can be played with any video player that supports MPEG-2 TS:  

mplayer test_out.ts  

**Notes**  
The baseband samples can be sent to USRP N210 on TX and received with another USRP N210 on RX. The decoder/encoder consumes a lot of processing power and therefore the realtime functionality will depend on the available computing power. On my computer - i7 Sandybridge 2600K - it is possible to send and receive in realtime the stream.  

The encoder/decoder currently works with gnuradio 3.7.2.x and 3.6 (if taken from Gnuradio_v_3_6 tag) but it is not dependent heavily on it. There is also a list of todo tasks in TODO.txt.  

The dvb-t module requires SSE2 SIMD instructions available in the processor. If the SSE2 is not available then illegal instruction will be seen at runtime.  

Demo grc files for encoding/decoding:  
- for 2k: QAM16 and rate 1/2: dvbt_tx_demo.grc and dvbt_rx_demo.grc.  
- for 2k: QAM64 and rate 7/8: dvbt_tx_demo_2k_QAM64_rate78.grc and dvbt_rx_demo_2k_QAM64_rate78.grc.  
- for 8k, QAM16 and rate 1/2: dvbt_tx_demo_8k.grc and dvbt_rx_demo_8k.grc.  
- for 8k, QAM64 and rate 7/8: dvbt_tx_demo_8k_QAM64_rate78.grc and dvbt_rx_demo_8k_QAM64_rate78.grc.  
- for 8k, QPSK and rate 7/8: dvbt_tx_demo_8k_QPSK_rate78.grc and dvbt_rx_demo_8k_QPSK_rate78.grc.  


**Build and run**  
git clone https://github.com/BogdanDIA/gr-dvbt.git  
cd gr-dvbt  
mkdir build  
cd build  
cmake ../  
make && sudo make install

if gr-dvbt is installed for the first time:

sudo ldconfig

