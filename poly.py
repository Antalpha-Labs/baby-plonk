from curve import Scalar
from enum import Enum
from numpy.polynomial import polynomial as P
import numpy as np

class Basis(Enum):
    LAGRANGE = 1
    MONOMIAL = 2


class Polynomial:
    values: list[Scalar]
    basis: Basis

    def __init__(self, values: list[Scalar], basis: Basis):
        assert all(isinstance(x, Scalar) for x in values)
        assert isinstance(basis, Basis)
        self.values = values
        self.basis = basis

    def __eq__(self, other):
        return (self.basis == other.basis) and (self.values == other.values)

    def __add__(self, other):
        if isinstance(other, Polynomial):
            assert self.basis == other.basis
            if (self.basis == Basis.LAGRANGE):
                assert len(self.values) == len(other.values)
                return Polynomial(
                    [x + y for x, y in zip(self.values, other.values)],
                    self.basis,
                )

            if (self.basis == Basis.MONOMIAL):
                res = P.polyadd(self.values, other.values)
                return Polynomial(
                    res,
                    self.basis,
                )
        else:
            assert isinstance(other, Scalar)
            if (self.basis == Basis.LAGRANGE):
                return Polynomial(
                    [x + other for x in self.values],
                    self.basis,
                )

            if (self.basis == Basis.MONOMIAL):
                res = P.polyadd(self.values, [other])
                return Polynomial(
                    res,
                    self.basis,
                )


    def __sub__(self, other):
        if isinstance(other, Polynomial):
            assert self.basis == other.basis
            if (self.basis == Basis.LAGRANGE):
                assert len(self.values) == len(other.values)
                return Polynomial(
                    [x - y for x, y in zip(self.values, other.values)],
                    self.basis,
                )

            if (self.basis == Basis.MONOMIAL):
                res = P.polysub(self.values, other.values)
                return Polynomial(
                    res,
                    self.basis,
                )
        else:
            assert isinstance(other, Scalar)
            if (self.basis == Basis.LAGRANGE):
                return Polynomial(
                    [x - other for x in self.values],
                    self.basis,
                )

            if (self.basis == Basis.MONOMIAL):
                res = P.polysub(self.values, [other])
                return Polynomial(
                    res,
                    self.basis,
                )

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            assert self.basis == other.basis
            if (self.basis == Basis.LAGRANGE):
                assert len(self.values) == len(other.values)
                res = [x * y for x, y in zip(self.values, other.values)]
            if (self.basis == Basis.MONOMIAL):
                c1 = self.values
                c2 = other.values
                res = P.polymul(c1,c2)

            return Polynomial(
                res,
                self.basis,
            )
        else:
            assert isinstance(other, Scalar)
            if (self.basis == Basis.LAGRANGE):
                return Polynomial(
                    [x * other for x in self.values],
                    self.basis,
                )

            if (self.basis == Basis.MONOMIAL):
                c1 = self.values
                c2 = [other]
                res = P.polymul(c1,c2)
                return Polynomial(
                    res,
                    self.basis,
                )

    # division without remainder
    def __truediv__(self, other):
        if isinstance(other, Polynomial):
            assert self.basis == other.basis
            if (self.basis == Basis.LAGRANGE):
                assert len(self.values) == len(other.values)
                return Polynomial(
                    [x / y for x, y in zip(self.values, other.values)],
                    self.basis,
                )
            if (self.basis == Basis.MONOMIAL):
                qx, rx = P.polydiv(self.values, other.values)
                # here we only consider the scenario of remainder is 0
                assert list(rx) == [0]

                return Polynomial(
                    qx,
                    self.basis,
                )
        else:
            assert isinstance(other, Scalar)
            if (self.basis == Basis.LAGRANGE):
                return Polynomial(
                    [x / other for x in self.values],
                    self.basis,
                )

            if (self.basis == Basis.MONOMIAL):
                c1 = self.values
                c2 = [other]
                quo, rx = P.polydiv(c1,c2)
                assert list(rx) == [0]
                return Polynomial(
                    quo,
                    self.basis,
                )

    def div_with_remainder(self, other):
        assert isinstance(other, Polynomial)
        assert self.basis == other.basis
        if (self.basis == Basis.LAGRANGE):
            assert len(self.values) == len(other.values)
            return Polynomial(
                [x / y for x, y in zip(self.values, other.values)],
                self.basis,
            )
        if (self.basis == Basis.MONOMIAL):
            qx, rx = P.polydiv(self.values, other.values)
            return Polynomial(
                qx,
                self.basis,
            ), Polynomial(
                rx,
                self.basis,
            )

    def shift(self, shift: int):
        assert self.basis == Basis.LAGRANGE
        assert shift < len(self.values)

        return Polynomial(
            self.values[shift:] + self.values[:shift],
            self.basis,
        )

    # Convenience method to do FFTs specifically over the subgroup over which
    # all of the proofs are operating
    def fft(self, inv=False):
        # Fast Fourier transform, used to convert between polynomial coefficients
        # and a list of evaluations at the roots of unity
        # See https://vitalik.ca/general/2019/05/12/fft.html
        def _fft(vals, modulus, roots_of_unity):
            if len(vals) == 1:
                return vals
            L = _fft(vals[::2], modulus, roots_of_unity[::2])
            R = _fft(vals[1::2], modulus, roots_of_unity[::2])
            o = [0] * len(vals)
            for i, (x, y) in enumerate(zip(L, R)):
                y_times_root = y * roots_of_unity[i]
                o[i] = (x + y_times_root) % modulus
                o[i + len(L)] = (x - y_times_root) % modulus
            return o

        roots = [x.n for x in Scalar.roots_of_unity(len(self.values))]
        o, nvals = Scalar.field_modulus, [x.n for x in self.values]
        if inv:
            assert self.basis == Basis.LAGRANGE
            # Inverse FFT
            invlen = Scalar(1) / len(self.values)
            reversed_roots = [roots[0]] + roots[1:][::-1]
            return Polynomial(
                [Scalar(x) * invlen for x in _fft(nvals, o, reversed_roots)],
                Basis.MONOMIAL,
            )
        else:
            assert self.basis == Basis.MONOMIAL
            # Regular FFT
            return Polynomial(
                [Scalar(x) for x in _fft(nvals, o, roots)], Basis.LAGRANGE
            )

    def ifft(self):
        return self.fft(True)

    # add two polynomial for all cases
    # this may be slower than the normal +
    def force_add(self, other):
        assert isinstance(other, Polynomial)
        if self.basis != other.basis:
            if self.basis == Basis.LAGRANGE:
                return self.ifft() + other
            else:
                return self + other.ifft()
        else:
            if self.basis == Basis.LAGRANGE and len(self.values) != len(other.values):
                return self.ifft() + other.ifft()
            else:
                return self + other

    # Given a polynomial expressed as a list of evaluations at roots of unity,
    # evaluate it at x directly, without using an FFT to covert to coeffs first
    # https://hackmd.io/@vbuterin/barycentric_evaluation
    def barycentric_eval(self, x: Scalar):
        assert self.basis == Basis.LAGRANGE

        order = len(self.values)
        roots_of_unity = Scalar.roots_of_unity(order)
        elem = x.n % Scalar.field_modulus
        if elem in roots_of_unity:
            return self.values[roots_of_unity.index(elem)]

        return (
            (Scalar(x) ** order - 1)
            / order
            * sum(
                [
                    value * root / (x - root)
                    for value, root in zip(self.values, roots_of_unity)
                ]
            )
        )

    # Evaluate at x directly for polynomial of MONOMIAL
    # This is inefficient, just for study usage
    def coeff_eval(self, x: Scalar):
        assert self.basis == Basis.MONOMIAL
        coeffs = self.values
        result = coeffs[0]
        x_pow = Scalar(1)
        for i in range(1, len(coeffs)):
            x_pow = x_pow * x
            result = result + coeffs[i] * x_pow
        return result

    def eval(self, x: Scalar):
        if self.basis == Basis.LAGRANGE:
            return self.barycentric_eval(x)
        else:
            return self.coeff_eval(x)

    def to_mononial(self):
        if self.basis == Basis.LAGRANGE:
            return self.ifft()
        else:
            return self

class PolyUtil:
    # f(X) = X - a
    def root_poly(self, x_val: Scalar) -> Polynomial:
        return Polynomial([-x_val, Scalar(1)], Basis.MONOMIAL)

    # f(X) = a
    def const_poly(self, x_val: Scalar) -> Polynomial:
        return Polynomial([x_val], Basis.MONOMIAL)

    # vanishing polynomial on multiplicative subgroup
    # z_H(X) = X^n - 1
    def vanishing_poly(self, n: int) -> Polynomial:
        return [Scalar(-1)] + [Scalar(0)] * (n - 1) + [Scalar(1)]

    # generate polynomial: X^n
    def x_exponent_poly(self, n: int) -> Polynomial:
        values = [Scalar(0)] * (n - 1) + [Scalar(1)]
        return Polynomial(values, Basis.MONOMIAL)

# construct MONOMIAL Polynomial with any X and Y values
# Note: do not use with FFT due to X is probably not multiplicative subgroup
class InterpolationPoly:
    n: int
    X: list[Scalar]
    Y: list[Scalar]
    def __init__(self, X: list[Scalar], Y: list[Scalar]):
        assert len(X) == len(Y), "Error: X should have the same length with Y"
        self.n = len(X)
        self.X = np.array(X)
        self.Y = np.array(Y)
        self.poly_util = PolyUtil()

    # z_H(X) = (X - self.X[0])(X - self.X[1])(X - self.X[2])...
    def vanishing_poly(self) -> Polynomial:
        v_poly = self.poly_util.const_poly(Scalar(1))
        for i in range(self.n):
            v_poly *= self.poly_util.root_poly(self.X[i])
        return v_poly

    # compute the derivative
    def vanishing_poly_diff(self) -> Polynomial:
        v_poly = self.vanishing_poly()
        v_diff_poly = self.poly_util.const_poly(Scalar(0))
        for i in range(self.n):
            v_diff_poly += v_poly / self.poly_util.root_poly(self.X[i])
        return v_diff_poly

    # Give i, return ith Lagrange polynomial L_i(X)
    # L_i(X) = z_H(X) / z_H'(a_i) / (X - a_i)
    def lagrange_poly(self, i: int) -> Polynomial:
        v_poly = self.vanishing_poly()
        v_diff_poly = self.vanishing_poly_diff()
        x_val = self.X[i]
        v_diff_poly_at_i = v_diff_poly.coeff_eval(x_val)
        x_root_poly = self.poly_util.root_poly(x_val)
        return v_poly / x_root_poly / v_diff_poly_at_i

    # f(X) = Σ(L_i(X) * y_i)
    def poly(self) -> Polynomial:
        poly = self.poly_util.const_poly(Scalar(0))
        for i in range(self.n):
            lagrange_poly = self.lagrange_poly(i)
            poly += lagrange_poly * self.Y[i]
        return poly