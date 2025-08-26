import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

def zplane(b, a, ax=None, filename=None):
    """Plot the complex z-plane given a transfer function.
    
    Parameters:
    - b: Coefficients of the numerator (zeros).
    - a: Coefficients of the denominator (poles).
    - ax: Matplotlib axis object where the plot will be drawn. If None, a new figure and axis are created.
    - filename: If provided, the plot will be saved to this file.
    
    Returns:
    - z: Zeros of the transfer function.
    - p: Poles of the transfer function.
    - k: Gain factor.
    """
    
    if ax is None:
        ax = plt.subplot(111)

    # Create the unit circle
    uc = patches.Circle((0, 0), radius=1, fill=False, color='black', ls='dashed')
    ax.add_patch(uc)

    # Normalize the coefficients if necessary
    if np.max(b) > 1:
        kn = np.max(b)
        b = b / float(kn)
    else:
        kn = 1

    if np.max(a) > 1:
        kd = np.max(a)
        a = a / float(kd)
    else:
        kd = 1

    # Get the poles and zeros
    p = np.roots(a)
    z = np.roots(b)
    k = kn / float(kd)

    # Plot the zeros and set marker properties    
    t1 = ax.plot(z.real, z.imag, 'go', ms=10)
    plt.setp(t1, markersize=10.0, markeredgewidth=1.0, markeredgecolor='k', markerfacecolor='g')

    # Plot the poles and set marker properties
    t2 = ax.plot(p.real, p.imag, 'rx', ms=10)
    plt.setp(t2, markersize=12.0, markeredgewidth=3.0, markeredgecolor='r', markerfacecolor='r')

    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Set the ticks
    r = 1.5
    ax.axis('scaled')
    ax.axis([-r, r, -r, r])
    ticks = [-1, -0.5, 0.5, 1]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    if filename is not None:
        plt.savefig(filename)

    return z, p, k
