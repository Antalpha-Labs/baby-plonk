from utils import *
import py_ecc.bn128 as b
from curve import ec_lincomb, G1Point, G2Point
from compiler.program import CommonPreprocessedInput
from verifier import VerificationKey
from dataclasses import dataclass
from poly import Polynomial, Basis

@dataclass
class Setup(object):
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( G,    xG,  ...,  x^{d-1}G ), where G is a generator of G_1
    powers_of_x: list[G1Point]
    # [x]₂ = xH, where H is a generator of G_2
    X2: G2Point

    @classmethod
    # tau: a random number whatever you choose
    # powers: length of SRS(structured reference string)
    def generate_srs(cls, powers: int, tau: int):
        print("Start to generate structured reference string")

        # Initialize powers_of_x with 0 values
        powers_of_x = [0] * powers
        # powers_of_x[0] =  b.G1 * tau**0 = b.G1
        # powers_of_x[1] =  b.G1 * tau**1 = powers_of_x[0] * tau
        # powers_of_x[2] =  b.G1 * tau**2 = powers_of_x[1] * tau
        # ...
        # powers_of_x[i] =  b.G1 * tau**i = powers_of_x[i - 1] * tau
        # TODO: generate powers_of_x
        powers_of_x = [b.multiply(b.G1, tau ** i) for i in range(powers)]
        # reference: https://github.com/sec-bit/learning-zkp/blob/master/plonk-intro-cn/5-plonk-polycom.md

        print("Generated G1 side, X^1 point: {}".format(powers_of_x[1]))

        X2 = b.multiply(b.G2, tau)
        print("Generated G2 side, X^1 point: {}".format(X2))

        # verify point is on the curve
        assert b.is_on_curve(powers_of_x[1], b.b)
        assert b.is_on_curve(X2, b.b2)
        # check pairing
        assert b.pairing(b.G2, powers_of_x[1]) == b.pairing(X2, b.G1)

        print("Finished to generate structured reference string")

        return cls(powers_of_x, X2)

    # Encodes the KZG commitment that evaluates to the given values in the group
    def commit(self, values: Polynomial) -> G1Point:
        if (values.basis == Basis.LAGRANGE):
            # inverse FFT from Lagrange basis to monomial basis
            coeffs = values.ifft().values
        elif (values.basis == Basis.MONOMIAL):
            coeffs = values.values
        if len(coeffs) > len(self.powers_of_x):
            raise Exception("Not enough powers in setup")
        return ec_lincomb([(s, x) for s, x in zip(self.powers_of_x, coeffs)])

    # Generate the verification key for this program with the given setup
    def verification_key(self, pk: CommonPreprocessedInput) -> VerificationKey:
        return VerificationKey(
            pk.group_order,
            self.commit(pk.QM),
            self.commit(pk.QL),
            self.commit(pk.QR),
            self.commit(pk.QO),
            self.commit(pk.QC),
            self.commit(pk.S1),
            self.commit(pk.S2),
            self.commit(pk.S3),
            self.X2,
            Scalar.root_of_unity(pk.group_order),
        )
