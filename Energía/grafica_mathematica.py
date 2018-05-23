#!/usr/bin/python
import matplotlib.pyplot as plt

mallados             = [5, 25, 50, 70, 100, 1000, 10000, 2*10**4]
tiempos_mathematica  = [0.008*186, 0.608*186, 9.66*186, 39.508*186, 181.252*186, 400*186, 1000*186, 2000*186]
tiempos_cnoparalell  = [0.001 * 100, 0.001 * 100, 0.001 *100, 0.001 *100, 0.001*100, 0.005*100, 1.178*100, 8.961*100]
tiempos_cparalell    = [230 * 0.001, 230 * 0.001, 230 * 0.001, 230 * 0.001, 230 * 0.001, 230 * 0.005, 0.624 * 230 , 7.754 * 230]
tiempos_gpu          = [130 * 0.440, 130 * 0.440, 130* 0.440, 130 * 0.440, 130 * 0.484,130 * 0.484, 0.997 * 130, 2.594 * 130]


plt.plot(mallados, tiempos_mathematica, label="Mathematica")
plt.plot(mallados, tiempos_cnoparalell, label="C No paralelo")
plt.plot(mallados, tiempos_cparalell, label="C Paralelo")
plt.plot(mallados, tiempos_gpu, label="GPU")

plt.legend()
plt.show()
