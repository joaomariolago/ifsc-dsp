from typing import Self, Optional
# External
import numpy
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
# Signals
import dsp_ifsc.sequence as sequence


class Signal:
    """
    The Signal class represents a signal in the discrete domain.

    Attributes:
        x: The sequence over n.
        n: The sequence over n.
    """

    def __init__(self, x: numpy.ndarray, n: numpy.ndarray):
        """
        Initializes the Signal class.
        Args:
            x: The sequence over n.
            n: The sequence over n.
        """

        if len(n) < len(x):
            raise ValueError("The n sequence must have the same length or more than x sequence")

        if len(n) > len(x):
            x = numpy.pad(x, (0, len(n) - len(x)))

        self.x = x
        self.n = n

    # Methods

    def _adjust_signals(self, other: Self) -> numpy.ndarray:
        """
        Adjusts the n sequence of two signals.
        Args:
            signal: The Signal class to be adjusted.

        Returns:
            n: The adjusted n sequence.
        """

        n = numpy.arange(min(self.n.min(0), other.n.min(0)), max(self.n.max(0), other.n.max(0)) + 1)
        y1 = numpy.zeros(len(n))
        y1[numpy.logical_and((n >= self.n.min(0)), (n <= self.n.max(0)))] = self.x
        y2 = numpy.zeros(len(n))
        y2[numpy.logical_and((n >= other.n.min(0)), (n <= other.n.max(0)))] = other.x

        return n, y1, y2


    def add(self, other: Self) -> 'Signal':
        """
        Implements y(n) = x1(n) + x2(n)
        Args:
            signal: The Signal class to be added.

        Returns:
            signal: The sum Signal class.
        """

        n, y1, y2 = self._adjust_signals(other)
        y = y1 + y2

        return Signal(y, n)


    def subtract(self, other: Self) -> 'Signal':
        """
        Implements y(n) = x1(n) - x2(n)
        Args:
            signal: The Signal class to be subtracted.

        Returns:
            signal: The subtracted Signal class.
        """

        n, y1, y2 = self._adjust_signals(other)
        y = y1 - y2

        return Signal(y, n)


    def multiply(self, other: Self) -> 'Signal':
        """
        Implements y(n) = x1(n) * x2(n)
        Args:
            signal: The Signal class to be multiplied.

        Returns:
            signal: The product Signal class.
        """

        n, y1, y2 = self._adjust_signals(other)
        y = y1 * y2

        return Signal(y, n)


    def divide(self, other: Self) -> 'Signal':
        """
        Implements y(n) = x1(n) / x2(n)
        Args:
            signal: The Signal class to be divided.

        Returns:
            signal: The division Signal class.
        """

        n, y1, y2 = self._adjust_signals(other)
        y = y1 / y2

        return Signal(y, n)


    def shift(self, k: int) -> 'Signal':
        """
        Implements y(n) = x(n - k)
        Args:
            k: The shift value.

        Returns:
            signal: The shifted Signal class.
        """

        y = self.x
        n = self.n + k

        return Signal(y, n)


    def fold(self) -> 'Signal':
        """
        Implements y(n) = x(-n)
        Returns:
            signal: The folded Signal class.
        """

        y = numpy.flip(self.x)
        n = -numpy.flip(self.n)

        return Signal(y, n)


    def negate(self) -> 'Signal':
        """
        Implements y(n) = -x(n)
        Returns:
            signal: The negated Signal class.
        """

        y = -self.x
        n = self.n

        return Signal(y, n)


    def convolution(self, other: Self) -> 'Signal':
        """
        Implements y(n) = x(n) * h(n)
        Args:
            signal: The Signal class to be convolved.

        Returns:
            signal: The convolved Signal class.
        """

        y = numpy.convolve(self.x, other.x)
        n = numpy.arange(self.n.min(0) + other.n.min(0), self.n.max(0) + other.n.max(0) + 1)

        return Signal(y, n)


    def correlation(self, other: Self) -> 'Signal':
        """
        Implements y(n) = x(n) * h(-n)
        Args:
            signal: The Signal class to be correlated.

        Returns:
            signal: The correlated Signal class.
        """

        y = numpy.correlate(self.x, other.x)
        n = numpy.arange(self.n.min(0) - other.n.max(0), self.n.max(0) - other.n.min(0) + 1)

        return Signal(y, n)


    def stem(self, plot: Optional[Axes] = None, auto_plot: Optional[bool] = True) -> Axes:
        """
        Plots the signal.
        """

        if plot is None:
            plot = plt.subplot()

        plot.stem(self.n, self.x)
        # Adjust grid and axes scales and intervals
        plot.grid(alpha = 0.3)
        x_ticks = numpy.round(numpy.linspace(self.n.min(), self.n.max(), 10), 2)
        y_ticks = numpy.round(numpy.linspace(self.x.min(), self.x.max() + self.x.max() / 10, 10), 2)
        plot.set_xticks(x_ticks)
        plot.set_yticks(y_ticks)
        # Details
        plot.set_xlabel('n')
        plot.set_ylabel('x(n)')
        plot.set_title('Signal')

        if auto_plot:
            plt.show()

        return plot


    # Overloads

    def __add__(self, other: Self | int | float) -> 'Signal':
        """
        Overloads the + operator.
        Args:
            signal: The Signal class to be added.

        Returns:
            signal: The sum Signal class.
        """

        if isinstance(other, (int, float)):
            return Signal(self.x + other, self.n)

        return self.add(other)


    __radd__ = __add__


    def __sub__(self, other: Self | int | float) -> 'Signal':
        """
        Overloads the - operator.
        Args:
            signal: The Signal class to be subtracted.

        Returns:
            signal: The subtracted Signal class.
        """

        if isinstance(other, (int, float)):
            return Signal(self.x - other, self.n)

        return self.subtract(other)


    def __rsub__(self, other: Self | int | float) -> 'Signal':
        """
        Overloads the - operator.
        Args:
            signal: The Signal class to be subtracted.

        Returns:
            signal: The subtracted Signal class.
        """

        return self.negate() + other


    def __mul__(self, other: Self | int | float) -> 'Signal':
        """
        Overloads the * operator.
        Args:
            signal: The Signal class to be multiplied.

        Returns:
            signal: The product Signal class.
        """

        if isinstance(other, (int, float)):
            return Signal(self.x * other, self.n)

        return self.multiply(other)


    __rmul__ = __mul__


    def __truediv__(self, other: Self | int | float) -> 'Signal':
        """
        Overloads the / operator.
        Args:
            signal: The Signal class to be divided.

        Returns:
            signal: The division Signal class.
        """

        if isinstance(other, (int, float)):
            return Signal(self.x / other, self.n)

        return self.divide(other)


    def __neg__(self) -> 'Signal':
        """
        Overloads the - operator.
        Returns:
            signal: The negative Signal class.
        """

        return self.negate()


    def __lshift__(self, k: int) -> 'Signal':
        """
        Overloads the << operator.
        Args:
            k: The shift value.

        Returns:
            signal: The shifted Signal class.
        """

        return self.shift(k)


    def __rshift__(self, k: int) -> 'Signal':
        """
        Overloads the >> operator.
        Args:
            k: The shift value.

        Returns:
            signal: The shifted Signal class.
        """

        return self.shift(-k)


    def __invert__(self) -> 'Signal':
        """
        Overloads the ~ operator.
        Returns:
            signal: The folded Signal class.
        """

        return self.fold()


    def __str__(self) -> str:
        """
        Overloads the str() function.
        Returns:
            str: The string representation of the Signal class.
        """

        return f"Signal(x={self.x}, n={self.n})"


    def __repr__(self) -> str:
        """
        Overloads the repr() function.
        Returns:
            str: The string representation of the Signal class.
        """

        return str(self)


    def __matmul__(self, other: Self) -> 'Signal':
        """
        Overloads the @ operator.
        Args:
            signal: The Signal class to be convolved.

        Returns:
            signal: The convolved Signal class.
        """

        return self.convolution(other)


    def __rmatmul__(self, other: Self) -> 'Signal':
        """
        Overloads the @ operator.
        Args:
            signal: The Signal class to be convolved.

        Returns:
            signal: The convolved Signal class.
        """

        return other.convolution(self)


    def __eq__(self, other: Self) -> bool:
        """
        Overloads the == operator.
        Args:
            signal: The Signal class to be compared.

        Returns:
            bool: True if the signals are equal, False otherwise.
        """

        return numpy.array_equal(self.x, other.x) and numpy.array_equal(self.n, other.n)


    def __ne__(self, other: Self) -> bool:
        """
        Overloads the != operator.
        Args:
            signal: The Signal class to be compared.

        Returns:
            bool: True if the signals are different, False otherwise.
        """

        return not (self == other)


    def __mod__(self, other: Self) -> 'Signal':
        """
        Overloads the % operator.
        Args:
            signal: The Signal class to be compared.

        Returns:
            signal: The Signal class.
        """

        return self.correlation(other)


    def __rmod__(self, other: Self) -> 'Signal':
        """
        Overloads the % operator.
        Args:
            signal: The Signal class to be compared.

        Returns:
            signal: The Signal class.
        """

        return other.correlation(self)


    # Static methods


    @staticmethod
    def from_scalar(scalar: float, n: numpy.ndarray) -> 'Signal':
        """
        Generates a scalar signal.
        Args:
            scalar: The scalar value.
            n: The sequence of n values.

        Returns:
            signal: The scalar Signal class.
        """

        x = scalar * numpy.ones(len(n))

        return Signal(x, n)


    @staticmethod
    def from_impulse(position: int, n: numpy.ndarray) -> 'Signal':
        """
        Generates an impulse signal.
        Args:
            position: The position of the impulse.
            n: The sequence of n values.

        Returns:
            signal: The impulse Signal class.
        """

        x, _n = sequence.impulse(position, n.min(), n.max())

        return Signal(x, _n)


    @staticmethod
    def from_step(position: int, n: numpy.ndarray):
        """
        Generates a step signal.
        Returns:
            signal: The step Signal class.
        """

        x, _n = sequence.step(position, n.min(), n.max())

        return Signal(x, _n)


    @staticmethod
    def from_ramp(slope: float, position: int, n: numpy.ndarray):
        """
        Generates a ramp signal.
        Returns:
            signal: The ramp Signal class.
        """

        x, _n = sequence.ramp(slope, position, n.min(), n.max())

        return Signal(x, _n)


    @staticmethod
    def from_exponential(amplitude: float, decay: float, position: int, n: numpy.ndarray):
        """
        Generates an exponential signal.
        Returns:
            signal: The exponential Signal class.
        """

        x, _n = sequence.exponential(amplitude, decay, position, n.min(), n.max())

        return Signal(x, _n)
