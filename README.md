###PROJECT: DVB-T implementation using gnuradio
**This is a opensource implementation of DVB-T encoder/decoder according to ETSI 300 744**  
More details can be seen here:  
http://yo3iiu.ro/blog/?p=1191  
http://yo3iiu.ro/blog/?p=1220  
http://yo3iiu.ro/blog/?p=1244  


###HowTo:
**Run TX**   
The simplest way to run DVB-T encoding is to use the gnuradio-companion flowgraphs apps/dvbt_tx_demo.grc. It will start with a MPEG-2 TS file and will eventually generate the 10Msps baseband samples.  
test.ts(MPEG-2 TS)->test.bin(BB)  

The parameters used are: 2k OFDM, FEC code 1/2, Modulation 16-QAM, Guard Interval 1/32

**Run RX** 
To run DVB-T decoding just run apps/dvbt_rx_demo.grc. It will take the baseband samples and turn them into MPEG-2 TS file.  
test.bin(BB)->test_out.ts(MPEG-2 TS)  

Then the output can be played with any video player that supports MPEG-2 TS:  

mplayer test_out.ts  

**Notes**  
The baseband samples can be sent to USRP N210 on TX and received with another USRP N210 on RX. The decoder/encoder consumes a lot of processing power and therefore the realtime functionality will depend on the available computing power. On my computer - Sandybridge 2600K it is possible to send and receive in realtime the stream.  



