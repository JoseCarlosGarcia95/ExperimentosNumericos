#!/usr/bin/python
import matplotlib.pyplot as plt

mallados             = [5, 25, 50, 70, 100, 1000, 10000, 2*10**4]
tiempos_mathematica  = [0.008, 0.608, 9.66, 39.508, 181.252, 400, 1000, 2000]
tiempos_cnoparalell  = [0.001, 0.001, 0.001, 0.001, 0.001, 0.005, 1.178, 8.961]
tiempos_cparalell    = [0.001, 0.001, 0.001, 0.001, 0.001, 0.005, 0.624, 7.754]
tiempos_gpu          = [0.440, 0.440, 0.440, 0.440, 0.484, 0.484, 0.997, 2.594]

#plt.plot(mallados, tiempos_mathematica, label="Mathematica")
plt.plot(mallados, tiempos_cnoparalell, label="C No paralelo")
plt.plot(mallados, tiempos_cparalell, label="C Paralelo")
plt.plot(mallados, tiempos_gpu, label="GPU")

plt.legend()
plt.show()
