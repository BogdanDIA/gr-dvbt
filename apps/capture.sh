#!/bin/bash

#./benchmark_rx.py -f 570e6 -m qam16 --fft-length=2048 --occupied-tones=1512 --cp-length=64 --snr=50 -W 10e6

sudo uhd_rx_cfile -a addr=192.168.10.2 -f 1280022000 --gain=19  --samp-rate=10e6 ./captureBB_3.bin

#sudo uhd_rx_cfile -f 570e6 --gain=0  --samp-rate=10e6 ./captureBB.bin

#./benchmark_rx.py -m qam16 --fft-length=2048 --occupied-tones=1512 --cp-length=64 --snr=40 -W 10e6 --from-file=./captureBB_benchmark.bin
