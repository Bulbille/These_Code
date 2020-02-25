#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import scipy.io.wavfile as wav
import numpy as np
import matplotlib.pyplot as plt
import os #pour les chemins de fichiers

# Lecture du fichier audio
a = wav.read('./sample.wav')
rate = a[0]
signal = a[1]

# Transformée de fourier
transfo = np.fft.fft(signal)
# Plage des fréquences
freq = np.fft.fftfreq(signal.size, d=1/rate)

#Plot
real = plt.subplot(211)
fourier = plt.subplot(212)

real.plot(np.arange(signal.size),signal)
fourier.plot(freq,transfo)

plt.show()
