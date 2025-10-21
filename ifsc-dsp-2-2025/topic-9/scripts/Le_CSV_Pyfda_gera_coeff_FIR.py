import pandas as pd
import numpy as np

def fir_header(fname_out, val):

    Ns= len(val)
    if ((Ns%2) != 0):
        val = np.append(val,0)
        Ns= len(val)
    val_fixo = np.round(val * 2**15)

    f = open(fname_out, 'wt')
    f.write('//define a FIR CMSIS-DSP coefficient array\n\n')
    f.write('#ifndef INCLUDE_COEFFS_FIR_H_\n')
    f.write('#define INCLUDE_COEFFS_FIR_H_\n\n')
    f.write('#include <stdint.h>\n\n')
    f.write('#include "filter.h"\n\n')
    f.write('#ifndef NUM_TAPS\n')
    f.write('#define NUM_TAPS %d\n' % Ns)
    f.write('#endif\n\n')
    f.write('/*************************************/\n');
    f.write('/*     FIR Filter Coefficients       */\n');
    f.write('const float32_t firCoeffs32[%d] = { ' % Ns)
    for k in range(Ns):
        if (k < Ns - 1):
            f.write(' %+-13e, ' % val[k] )
        else:
            f.write(' %+-13e ' % val[k] )

    f.write('};\n\n')
    f.write('const q15_t firCoeffsQ15[%d] = { ' % Ns)
    for k in range(Ns):
        if (k < Ns - 1):
            f.write(' %d, ' % val_fixo[k] )
        else:
            f.write(' %d ' % val_fixo[k] )

    f.write('};\n')

    f.write('FilterTypeDef filterType=FIR_FLOAT32;\n');

    f.write('/***********************************/\n\n')
    f.write('#endif /* INCLUDE_COEFFS_FIR_H_ */')
    f.close()




# Lê um arquivo .CSV e o transforma em um numpy array
# gravar table orientation = Vertical no pyfda
dados = pd.read_csv('coeffs_FIR.csv', header = None)
valores = dados.values
valores = valores[:,0]


print(valores)


from scipy import signal
import matplotlib.pyplot as plt


w, h = signal.freqz(valores,[1])

fs_kit = 48000
fig, ax1 = plt.subplots()
ax1.set_title('Digital filter frequency response')
ax1.plot((w/np.pi)*(fs_kit/2), 20 * np.log10(abs(h)), 'b')
ax1.set_ylabel('Amplitude [dB]', color='b')
ax1.set_xlabel('Frequency [rad/sample]')
ax1.grid()
#ax1.set_ylim([-80, 5])
plt.show()

fir_header('coeffs_FIR.h',valores)
