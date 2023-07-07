from utils import Scalar
from curve import G1Point
from merlin import MerlinTranscript
from py_ecc.secp256k1.secp256k1 import bytes_to_int
from dataclasses import dataclass


@dataclass
class Message1:
    # [a(x)]₁ (commitment to left wire polynomial)
    a_1: G1Point
    # [b(x)]₁ (commitment to right wire polynomial)
    b_1: G1Point
    # [c(x)]₁ (commitment to output wire polynomial)
    c_1: G1Point


@dataclass
class Message2:
    # [z(x)]₁ (commitment to permutation polynomial)
    z_1: G1Point


@dataclass
class Message3:
    # [quot(x)]₁ (commitment to the quotient polynomial t(X))
    W_quot: G1Point


@dataclass
class Message4:
    # Evaluation of a(X) at evaluation challenge ζ
    a_eval: Scalar
    # Evaluation of b(X) at evaluation challenge ζ
    b_eval: Scalar
    # Evaluation of c(X) at evaluation challenge ζ
    c_eval: Scalar
    # Evaluation of QL(X) at evaluation challenge ζ
    ql_eval: Scalar
    # Evaluation of QR(X) at evaluation challenge ζ
    qr_eval: Scalar
    # Evaluation of QM(X) at evaluation challenge ζ
    qm_eval: Scalar
    # Evaluation of QO(X) at evaluation challenge ζ
    qo_eval: Scalar
    # Evaluation of QC(X) at evaluation challenge ζ
    qc_eval: Scalar
    # Evaluation of the first permutation polynomial S_σ1(X) at evaluation challenge ζ
    s1_eval: Scalar
    # Evaluation of the second permutation polynomial S_σ2(X) at evaluation challenge ζ
    s2_eval: Scalar
    # Evaluation of the second permutation polynomial S_σ3(X) at evaluation challenge ζ
    s3_eval: Scalar
    # Evaluation of the permutation polynomial z(X) at the evaluation challenge ζ
    z_eval: Scalar
    # Evaluation of the shifted permutation polynomial z(X) at the shifted evaluation challenge ζω
    z_shifted_eval: Scalar
    # Evaluation of Quotient Polynomial Quot(X) at evaluation challenge ζ
    quot_eval: Scalar


@dataclass
class Message5:
    # (commitment to the opening proof polynomial)
    W_a: G1Point
    W_a_quot: G1Point
    W_b: G1Point
    W_b_quot: G1Point
    W_c: G1Point
    W_c_quot: G1Point
    W_z: G1Point
    W_z_quot: G1Point
    W_zw: G1Point
    W_zw_quot: G1Point
    W_t: G1Point
    W_t_quot: G1Point


class Transcript(MerlinTranscript):
    def append(self, label: bytes, item: bytes) -> None:
        self.append_message(label, item)

    def append_scalar(self, label: bytes, item: Scalar):
        self.append_message(label, item.n.to_bytes(32, "big"))

    def append_point(self, label: bytes, item: G1Point):
        self.append_message(label, item[0].n.to_bytes(32, "big"))
        self.append_message(label, item[1].n.to_bytes(32, "big"))

    def get_and_append_challenge(self, label: bytes) -> Scalar:
        while True:
            challenge_bytes = self.challenge_bytes(label, 255)
            f = Scalar(bytes_to_int(challenge_bytes))
            if f != Scalar.zero():  # Enforce challenge != 0
                self.append(label, challenge_bytes)
                return f

    def round_1(self, message: Message1) -> tuple[Scalar, Scalar]:
        self.append_point(b"a_1", message.a_1)
        self.append_point(b"b_1", message.b_1)
        self.append_point(b"c_1", message.c_1)

        # The first two Fiat-Shamir challenges
        beta = self.get_and_append_challenge(b"beta")
        gamma = self.get_and_append_challenge(b"gamma")

        return beta, gamma

    def round_2(self, message: Message2) -> tuple[Scalar, Scalar]:
        self.append_point(b"z_1", message.z_1)

        alpha = self.get_and_append_challenge(b"alpha")
        # This value could be anything, it just needs to be unpredictable. Lets us
        # have evaluation forms at cosets to avoid zero evaluations, so we can
        # divide polys without the 0/0 issue
        fft_cofactor = self.get_and_append_challenge(b"fft_cofactor")

        return alpha, fft_cofactor

    def round_3(self, message: Message3) -> Scalar:
        self.append_point(b"W_quot", message.W_quot)

        zeta = self.get_and_append_challenge(b"zeta")
        return zeta

    def round_4(self, message: Message4) -> Scalar:
        self.append_scalar(b"a_eval", message.a_eval)
        self.append_scalar(b"b_eval", message.b_eval)
        self.append_scalar(b"c_eval", message.c_eval)
        self.append_scalar(b"s1_eval", message.s1_eval)
        self.append_scalar(b"s2_eval", message.s2_eval)
        self.append_scalar(b"s3_eval", message.s3_eval)
        self.append_scalar(b"z_shifted_eval", message.z_shifted_eval)

        v = self.get_and_append_challenge(b"v")
        return v

    def round_5(self, message: Message5) -> Scalar:
        self.append_point(b"W_a", message.W_a)
        self.append_point(b"W_a_quot", message.W_a_quot)
        self.append_point(b"W_b", message.W_b)
        self.append_point(b"W_b_quot", message.W_b_quot)
        self.append_point(b"W_c", message.W_c)
        self.append_point(b"W_c_quot", message.W_c_quot)
        self.append_point(b"W_z", message.W_z)
        self.append_point(b"W_z_quot", message.W_z_quot)
        self.append_point(b"W_zw", message.W_zw)
        self.append_point(b"W_zw_quot", message.W_zw_quot)
        self.append_point(b"W_t", message.W_t)
        self.append_point(b"W_t_quot", message.W_t_quot)

        u = self.get_and_append_challenge(b"u")
        return u
