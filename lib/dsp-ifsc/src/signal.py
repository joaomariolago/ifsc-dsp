import numpy
from typing import Tuple

def add(x1: numpy.ndarray, n1: numpy.ndarray, x2: numpy.ndarray, n2: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """
    Implements y(n) = x1(n) + x2(n)
    Args:
        x1: The first sequence over n1.
        n1: The first sequence over n1.
        x2: The second sequence over n2.
        n2: The second sequence over n2.

    Returns:
        y: The sum sequence over n, which includes n1 and n2.
        n: The duration of y(n).
    """

    n = numpy.arange(min(n1.min(0), n2.min(0)), max(n1.max(0), n2.max(0)) + 1)
    y1 = numpy.zeros(len(n))
    y1[numpy.logical_and((n >= n1.min(0)), (n <= n1.max(0)))] = x1
    y2 = numpy.zeros(len(n))
    y2[numpy.logical_and((n >= n2.min(0)), (n <= n2.max(0)))] = x2
    y = y1 + y2

    return y, n

def multiply(x1: numpy.ndarray, n1: numpy.ndarray, x2: numpy.ndarray, n2: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """
    Implements y(n) = x1(n) * x2(n)
    Args:
        x1: The first sequence over n1.
        n1: The first sequence over n1.
        x2: The second sequence over n2.
        n2: The second sequence over n2.

    Returns:
        y: The product sequence over n, which includes n1 and n2.
        n: The duration of y(n).
    """

    n = numpy.arange(min(n1.min(0), n2.min(0)), max(n1.max(0), n2.max(0)) + 1)
    y1 = numpy.zeros(len(n))
    y1[numpy.logical_and((n >= n1.min(0)), (n <= n1.max(0)))] = x1
    y2 = numpy.zeros(len(n))
    y2[numpy.logical_and((n >= n2.min(0)), (n <= n2.max(0)))] = x2
    y = y1 * y2

    return y, n

def shift(x: numpy.ndarray, n: numpy.ndarray, k: int) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """
    Implements y(n) = x(n - k)
    Args:
        x: The sequence over n.
        n: The sequence over n.
        k: The shift value.

    Returns:
        y: The shifted sequence over n.
        n: The shifted sequence over n.
    """

    y = x
    n = n + k

    return y, n

def fold(x: numpy.ndarray, n: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """
    Implements y(n) = x(-n)
    Args:
        x: The sequence over n.
        n: The sequence over n.

    Returns:
        y: The folded sequence over n.
        n: The folded sequence over n.
    """

    y = numpy.flip(x)
    n = -numpy.flip(n)

    return y, n

def convolution(x: numpy.ndarray, n: numpy.ndarray, h: numpy.ndarray, m: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """
    Implements y(n) = x(n) * h(n)
    Args:
        x: The first sequence over n.
        n: The first sequence over n.
        h: The second sequence over m.
        m: The second sequence over m.

    Returns:
        y: The convolution sequence over n.
        n: The convolution sequence over n.
    """

    y = numpy.convolve(x, h)
    n = numpy.arange(n.min(0) + m.min(0), n.max(0) + m.max(0) + 1)

    return y, n

def correlation(x: numpy.ndarray, n: numpy.ndarray, h: numpy.ndarray, m: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """
    Implements y(n) = x(n) * h(-n)
    Args:
        x: The first sequence over n.
        n: The first sequence over n.
        h: The second sequence over m.
        m: The second sequence over m.

    Returns:
        y: The correlation sequence over n.
        n: The correlation sequence over n.
    """

    y = numpy.correlate(x, h)
    n = numpy.arange(n.min(0) - m.max(0), n.max(0) - m.min(0) + 1)

    return y, n
