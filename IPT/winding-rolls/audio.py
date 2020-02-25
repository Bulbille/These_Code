#!/usr/bin/python
# -*- coding: utf8

import scipy.io.wavfile as wav
import numpy as np
import matplotlib.pyplot as plt

# Lecture du fichier audio
a = wav.read('./sample.wav')
rate = a[0]
signal = a[1]
tot = signal.size
signal = signal[:int(tot/10)]

# Transformée de fourier
transfo = np.fft.rfft(signal)
# Plage des fréquences
freq = np.fft.rfftfreq(signal.size, d=1/rate)

#Plot
real = plt.subplot(211)
fourier = plt.subplot(212)

real.plot(np.arange(signal.size)/rate,signal)
fourier.plot(freq,transfo)

real.set_xlabel('Temps (s)')
real.set_xlim([0,signal.size/rate])

fourier.set_xlabel('Frequence (Hz)')
fourier.set_xlim([0,500])

plt.show()
