# funções para filtros analógicos e digitais

import numpy as np
from scipy import signal

def u_buttap(N,Omegac):
    # Unnormalized Butterworth Analog Lowpass Filter Prototype
    # --------------------------------------------------------
    # [b,a] = u_buttap(N,Omegac);
    # b = numerator polynomial coefficients of Ha(s)
    # a = denominator polynomial coefficients of Ha(s)
    # N = Order of the Butterworth Filter
    # Omegac = Cutoff frequency in radians/sec
    [z,p,k] = signal.buttap(N)
    p = p*Omegac
    k = k*Omegac**N
    B = np.real(np.poly(z))
    b0 = k
    b = k*B
    if type(b)!=np.ndarray: b=np.array([b])
    a = np.real(np.poly(p))
    if type(a)!=np.ndarray: a=np.array([a])
    return [b,a]

def u_chb1ap(N,Rp,Omegac):
    # Unnormalized Chebyshev-1 Analog Lowpass Filter Prototype
    # --------------------------------------------------------
    # [b,a] = u_chb1ap(N,Rp,Omegac)
    # b = numerator polynomial coefficients of Ha(s)
    # a = denominator polynomial coefficients of Ha(s)
    # N = Order of the Elliptic Filter
    # Rp = Passband Ripple in dB; Rp > 0
    # Omegac = Cutoff frequency in radians/sec
    [z,p,k] = signal.cheb1ap(N,Rp)
    a = np.real(np.poly(p))
    aNn = a[int(N)]
    p = p*Omegac
    a = np.real(np.poly(p))
    aNu = a[int(N)]  
    k = k*aNu/aNn
    b0 = k
    B = np.real(np.poly(z))
    b = k*B
    if type(b)!=np.ndarray: b=np.array([b])
    a = np.real(np.poly(p))
    if type(a)!=np.ndarray: a=np.array([a])
    return [b,a]

def afd_butt(Wp,Ws,Rp,As):
    # Analog Lowpass Filter Design: Butterworth
    #  -----------------------------------------
    #  [b,a] = afd_butt(Wp,Ws,Rp,As);
    #   b = Numerator coefficients of Ha(s)
    #   a = Denominator coefficients of Ha(s)
    #  Wp = Passband edge frequency in rad/sec; Wp > 0
    #  Ws = Stopband edge frequency in rad/sec; Ws > Wp > 0
    # Rp = Passband ripple in +dB; (Rp > 0)
    # As = Stopband attenuation in +dB; (As > 0)
    #
    if Wp <= 0:
        print("error: Passband edge must be larger than 0")
    if Ws <= Wp:
        print("error: Stopband edge must be larger than Passband edge")
    if (Rp <= 0) or (As < 0):
        print("error: PB ripple and/or SB attenuation ust be larger than 0")
    N = np.ceil((np.log10((10**(Rp/10)-1)/(10**(As/10)-1)))/(2*np.log10(Wp/Ws)))
    print("\n*** Butterworth Filter Order =",int(N))
    OmegaC = Wp/((10**(Rp/10)-1)**(1/(2*N)))
    [b,a] = u_buttap(N,OmegaC)
    return b,a

def afd_chb1(Wp,Ws,Rp,As):
    # Analog Lowpass Filter Design: Chebyshev-1
    #  -----------------------------------------
    #  [b,a] = afd_chb1(Wp,Ws,Rp,As);
    #   b = Numerator coefficients of Ha(s)
    #   a = Denominator coefficients of Ha(s)
    #  Wp = Passband edge frequency in rad/sec; Wp > 0
    #  Ws = Stopband edge frequency in rad/sec; Ws > Wp > 0
    # Rp = Passband ripple in +dB; (Rp > 0)
    # As = Stopband attenuation in +dB; (As > 0)
    #
    if Wp <= 0:
        print("error: Passband edge must be larger than 0")
    if Ws <= Wp:
        print("error: Stopband edge must be larger than Passband edge")
    if (Rp <= 0) or (As < 0):
        print("error: PB ripple and/or SB attenuation ust be larger than 0")
    ep = np.sqrt(10**(Rp/10)-1)
    A = 10**(As/20)
    OmegaC = Wp
    OmegaR = Ws/Wp
    g = np.sqrt(A*A-1)/ep
    N = np.ceil(np.log10(g+np.sqrt(g*g-1))/np.log10(OmegaR+np.sqrt(OmegaR*OmegaR-1)))
    print("\n*** Chebyshev-1 Filter Order =",int(N))
    [b,a]=u_chb1ap(N,Rp,OmegaC)
    return b,a


def freqs_m(b,a,wmax):
    # Computation of s-domain frequency response: Modified version
    # ------------------------------------------------------------
    # [db,mag,pha,w] = freqs_m(b,a,wmax);
    #   db = Relative magnitude in db over [0 to wmax]
    #  mag = Absolute magnitude over [0 to wmax]
    #  pha = Phase response in radians over [0 to wmax]
    #    w = array of 500 frequency samples between [0 to wmax]
    #    b = Numerator polynomial coefficents of Ha(s)
    #    a = Denominator polynomial coefficents of Ha(s)
    # wmax = Maximum frequency in rad/sec over which response is desired
    #
    #w = np.arange(501)*wmax/500
    w,H = signal.freqs(b,a,worN=501)
    mag = np.abs(H)
    db = 20*np.log10((mag+np.spacing(1))/np.max(mag))
    pha = np.angle(H)
    return db,mag,pha,w

def imp_invr(c,d,T):
    # Impulse Invariance Transformation from Analog to Digital Filter
    # ---------------------------------------------------------------
    # [b,a] = imp_invr(c,d,T)
    # b = Numerator polynomial in z^(-1) of the digital filter
    # a = Denominator polynomial in z^(-1) of the digital filter
    # c = Numerator polynomial in s of the analog filter
    # d = Denominator polynomial in s of the analog filter
    # T = Sampling (transformation) parameter
    #
    [R,p,k] = signal.residue(c,d)
    p = np.exp(p*T)
    [b,a] = signal.invresz(R,p,k)
    b = np.real(b)
    a = np.real(a)
    return b,a

def freqz_m(b,a):
    # Modified version of freqz subroutine
    # ------------------------------------
    # [db,mag,pha,grd,w] = freqz_m(b,a)
    #  db = Relative magnitude in dB computed over 0 to pi radians
    # mag = absolute magnitude computed over 0 to pi radians 
    # pha = Phase response in radians over 0 to pi radians
    # grd = Group delay over 0 to pi radians
    #   w = 501 frequency samples between 0 to pi radians
    #   b = numerator polynomial of H(z)   (for FIR: b=h)
    #   a = denominator polynomial of H(z) (for FIR: a=[1])
    #
    [w, H] = signal.freqz(b,a,worN=2048)
    db = 20*np.log10((np.abs(H))/np.max(np.abs(H)))
    w,grd = signal.group_delay((b,a),w)
    mag = np.abs(H)
    pha = np.angle(H)
    return db,mag,pha,grd,w


def QCoeff(x,N):
    #  [y,L,B] = QCoeff(x,N)
    #  Coefficient Quantization using N=1+L+B bit Representation
    #  with Rounding operation
    #  y: quantized array (same dim as x)
    #  L: number of integer bits
    #  B: number of fractional bits
    #  x: a scalar, vector, or matrix
    #  N: total number of bits
    xm = np.abs(x)
    eps = np.spacing(1)
    L = np.max(np.maximum(0,np.fix(np.log2(xm[:]+eps)+1)))   # Integer bits
    if L > N:
        print("errmsg = *** N deve ser pelo menos ", int(L))
    B = N-L   # Fractional bits
    y = xm/(2**L)
    y = np.round(y*(2**N))   # Rounding to N bits
    y = np.sign(x)*y*(2**(-B))   # L+B+1 bit representation
    return y,L,B

