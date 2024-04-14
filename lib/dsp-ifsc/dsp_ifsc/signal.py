import numpy
# Signals
import sequence


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

        self.x = x
        self.n = n

    # Methods

    def add(self, other: 'Signal') -> 'Signal':
        """
        Implements y(n) = x1(n) + x2(n)
        Args:
            signal: The Signal class to be added.

        Returns:
            signal: The sum Signal class.
        """

        n = numpy.arange(min(self.n.min(0), other.n.min(0)), max(self.n.max(0), other.n.max(0)) + 1)
        y1 = numpy.zeros(len(n))
        y1[numpy.logical_and((n >= self.n.min(0)), (n <= self.n.max(0)))] = self.x
        y2 = numpy.zeros(len(n))
        y2[numpy.logical_and((n >= other.n.min(0)), (n <= other.n.max(0)))] = other.x
        y = y1 + y2

        return y, n

    def multiply(self, other: 'Signal') -> 'Signal':
        """
        Implements y(n) = x1(n) * x2(n)
        Args:
            signal: The Signal class to be multiplied.

        Returns:
            signal: The product Signal class.
        """

        n = numpy.arange(min(self.n.min(0), other.n.min(0)), max(self.n.max(0), other.n.max(0)) + 1)
        y1 = numpy.zeros(len(n))
        y1[numpy.logical_and((n >= self.n.min(0)), (n <= self.n.max(0)))] = self.x
        y2 = numpy.zeros(len(n))
        y2[numpy.logical_and((n >= other.n.min(0)), (n <= other.n.max(0)))] = other.x
        y = y1 * y2

        return y, n


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

        return y, n


    def fold(self) -> 'Signal':
        """
        Implements y(n) = x(-n)
        Returns:
            signal: The folded Signal class.
        """

        y = numpy.flip(self.x)
        n = -numpy.flip(self.n)

        return y, n


    def negate(self) -> 'Signal':
        """
        Implements y(n) = -x(n)
        Returns:
            signal: The negated Signal class.
        """

        y = -self.x
        n = self.n

        return y, n


    def convolution(self, other: 'Signal') -> 'Signal':
        """
        Implements y(n) = x(n) * h(n)
        Args:
            signal: The Signal class to be convolved.

        Returns:
            signal: The convolved Signal class.
        """

        y = numpy.convolve(self.x, other.x)
        n = numpy.arange(self.n.min(0) + other.n.min(0), self.n.max(0) + other.n.max(0) + 1)

        return y, n


    def correlation(self, other: 'Signal') -> 'Signal':
        """
        Implements y(n) = x(n) * h(-n)
        Args:
            signal: The Signal class to be correlated.

        Returns:
            signal: The correlated Signal class.
        """

        y = numpy.correlate(self.x, other.x)
        n = numpy.arange(self.n.min(0) - other.n.max(0), self.n.max(0) - other.n.min(0) + 1)

        return y, n

    # Overloads

    def __add__(self, other: 'Signal') -> 'Signal':
        """
        Overloads the + operator.
        Args:
            signal: The Signal class to be added.

        Returns:
            signal: The sum Signal class.
        """

        return self.add(other)


    def __mul__(self, other: 'Signal') -> 'Signal':
        """
        Overloads the * operator.
        Args:
            signal: The Signal class to be multiplied.

        Returns:
            signal: The product Signal class.
        """

        return self.multiply(other)


    def __sub__(self, other: 'Signal') -> 'Signal':
        """
        Overloads the - operator.
        Args:
            signal: The Signal class to be subtracted.

        Returns:
            signal: The subtracted Signal class.
        """

        return self.add(-other)


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


    # Static methods

    @staticmethod
    def impulse(position: int, n_start: int, n_end: int) -> 'Signal':
        """
        Generates an impulse signal.
        Returns:
            signal: The impulse Signal class.
        """

        x, n = sequence.impulse(position, n_start, n_end)

        return Signal(x, n)


    @staticmethod
    def from_step(position: int, n_start: int, n_end: int):
        """
        Generates a step signal.
        Returns:
            signal: The step Signal class.
        """

        x, n = sequence.step(position, n_start, n_end)

        return Signal(x, n)


    @staticmethod
    def from_ramp(slope: float, position: int, n_start: int, n_end: int):
        """
        Generates a ramp signal.
        Returns:
            signal: The ramp Signal class.
        """

        x, n = sequence.ramp(slope, position, n_start, n_end)

        return Signal(x, n)


    @staticmethod
    def from_exponential(amplitude: float, decay: float, position: int, n_start: int, n_end: int):
        """
        Generates an exponential signal.
        Returns:
            signal: The exponential Signal class.
        """

        x, n = sequence.exponential(amplitude, decay, position, n_start, n_end)

        return Signal(x, n)
