import numpy as np
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt

def iir_sos_header(fname_out, sos):
    Ns = sos.shape[0]
    with open(fname_out, 'wt') as f:
        f.write('//define a IIR SOS CMSIS-DSP coefficient array\n\n')
        f.write('#ifndef INCLUDE_COEFFS_IIR_H_\n')
        f.write('#define INCLUDE_COEFFS_IIR_H_\n\n')
        f.write('#include <stdint.h>\n\n')
        f.write('#ifndef NUM_STAGES\n')
        f.write('#define NUM_STAGES %d\n' % Ns)
        f.write('#endif\n\n')
        f.write('/*********************************************************/\n')
        f.write('/*                     IIR SOS Filter Coefficients       */\n')
        f.write('const float32_t iirCoeffsF32[%d] = { // b0,b1,b2,a1,a2,... by stage\n' % (5*Ns))
        for k in range(Ns):
            b0, b1, b2, a0, a1, a2 = sos[k]
            # CMSIS expects denom = 1 + a1 z^-1 + a2 z^-2 (NO sign flip here)
            line = f'    {b0:+.12e}, {b1:+.12e}, {b2:+.12e},\n' \
                   f'    {a1:+.12e}, {a2:+.12e}'
            f.write(line + (',\n' if k < Ns-1 else '\n'))
        f.write('};\n')
        f.write('/*********************************************************/\n\n')
        f.write('#endif /* INCLUDE_COEFFS_IIR_H_ */\n')

# ---- Load your b/a CSV (vertical orientation from pyfdax/pyFDA) ----
dados = pd.read_csv('coeffs_IIR.csv', header=None).values
# dados is 9x2: first col are b_k, second col are a_k
b = dados[:, 0]          # b0..b8
a = dados[:, 1]          # a0..a8  (a0 should be 1.0)
if not np.isclose(a[0], 1.0):
    # normalize to make a0 = 1 just in case
    b = b / a[0]
    a = a / a[0]

# ---- Convert full polynomials to SOS ----
sos = signal.tf2sos(b, a)  # shape (n_sections, 6) with rows [b0 b1 b2 a0 a1 a2]

# ---- Check and plot frequency response ----
w, h = signal.sosfreqz(sos, worN=2048)
fs = 48000
plt.figure()
plt.title('Digital filter frequency response')
plt.plot((w/np.pi)*(fs/2), 20*np.log10(np.maximum(np.abs(h), 1e-12)))
plt.ylabel('Amplitude [dB]')
plt.xlabel('Frequency [Hz]')
plt.grid(True)
plt.ylim([-50, 5])
plt.show()

# ---- Emit CMSIS header in {b0,b1,b2,a1,a2} packing ----
iir_sos_header('coeffs_IIR.h', sos.astype(np.float32))
