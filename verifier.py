import py_ecc.bn128 as b
from utils import *
from dataclasses import dataclass
from curve import *
from transcript import Transcript
from poly import Polynomial, Basis


@dataclass
class VerificationKey:
    """Verification key"""

    group_order: int
    # [q_M(x)]₁ (commitment to multiplication selector polynomial)
    Qm: G1Point
    # [q_L(x)]₁ (commitment to left selector polynomial)
    Ql: G1Point
    # [q_R(x)]₁ (commitment to right selector polynomial)
    Qr: G1Point
    # [q_O(x)]₁ (commitment to output selector polynomial)
    Qo: G1Point
    # [q_C(x)]₁ (commitment to constants selector polynomial)
    Qc: G1Point
    # [S_σ1(x)]₁ (commitment to the first permutation polynomial S_σ1(X))
    S1: G1Point
    # [S_σ2(x)]₁ (commitment to the second permutation polynomial S_σ2(X))
    S2: G1Point
    # [S_σ3(x)]₁ (commitment to the third permutation polynomial S_σ3(X))
    S3: G1Point
    # [x]₂ = xH, where H is a generator of G_2
    X_2: G2Point
    # nth root of unity, where n is the program's group order.
    w: Scalar

    # More optimized version that tries hard to minimize pairings and
    # elliptic curve multiplications, but at the cost of being harder
    # to understand and mixing together a lot of the computations to
    # efficiently batch them
    def verify_proof(self, group_order: int, pf, public=[]) -> bool:
        # 4. Compute challenges
        beta, gamma, alpha, zeta, v, u = self.compute_challenges(pf)
        proof = pf.flatten()

        # 5. Compute zero polynomial evaluation Z_H(ζ) = ζ^n - 1
        ZH_ev = zeta**group_order - 1

        # 6. Compute Lagrange polynomial evaluation L_0(ζ)
        L0_ev = ZH_ev / (group_order * (zeta - 1))

        # 7. Compute public input polynomial evaluation PI(ζ).
        PI = Polynomial(
            [Scalar(-x) for x in public]
            + [Scalar(0) for _ in range(group_order - len(public))],
            Basis.LAGRANGE,
        )
        PI_ev = PI.barycentric_eval(zeta)

        # verify A with KZG10 commitment
        W_a = proof["W_a"]
        W_a_quot = proof["W_a_quot"]
        a_eval = proof["a_eval"]
        ec_comb_a = ec_lincomb(
                [
                    (W_a, 1),
                    (W_a_quot, zeta),
                    (b.G1, -a_eval),
                ]
            )
        assert b.pairing(self.X_2, W_a_quot) == b.pairing(b.G2, ec_comb_a)
        print("Done KZG10 commitment check for A polynomial")

        # # verify A with KZG10 commitment
        W_b = proof["W_b"]
        W_b_quot = proof["W_b_quot"]
        b_eval = proof["b_eval"]
        ec_comb_b = ec_lincomb(
                [
                    (W_b, 1),
                    (W_b_quot, zeta),
                    (b.G1, -b_eval),
                ]
            )

        assert b.pairing(self.X_2, W_b_quot) == b.pairing(b.G2, ec_comb_b)
        print("Done KZG10 commitment check for B polynomial")

        # # verify A with KZG10 commitment
        W_c = proof["W_c"]
        W_c_quot = proof["W_c_quot"]
        c_eval = proof["c_eval"]
        ec_comb_c = ec_lincomb(
                [
                    (W_c, 1),
                    (W_c_quot, zeta),
                    (b.G1, -c_eval),
                ]
            )

        assert b.pairing(self.X_2, W_c_quot) == b.pairing(b.G2, ec_comb_c)
        print("Done KZG10 commitment check for C polynomial")

        # # verify Z with KZG10 commitment
        W_z = proof["W_z"]
        W_z_quot = proof["W_z_quot"]
        z_eval = proof["z_eval"]
        ec_comb_z = ec_lincomb(
                [
                    (W_z, 1),
                    (W_z_quot, zeta),
                    (b.G1, -z_eval),
                ]
            )

        assert b.pairing(self.X_2, W_z_quot) == b.pairing(b.G2, ec_comb_z)
        print("Done KZG10 commitment check for Z polynomial")

        # verify ZW with KZG10 commitment
        W_zw = proof["W_zw"]
        W_zw_quot = proof["W_zw_quot"]
        zw_eval = proof["z_shifted_eval"]
        ec_comb_zw = ec_lincomb(
                [
                    (W_zw, 1),
                    (W_zw_quot, zeta),
                    (b.G1, -zw_eval),
                ]
            )

        assert b.pairing(self.X_2, W_zw_quot) == b.pairing(b.G2, ec_comb_zw)
        print("Done KZG10 commitment check for ZW polynomial")

        # verify T with KZG10 commitment
        W_t = proof["W_t"]
        W_t_quot = proof["W_t_quot"]
        t_eval = proof["quot_eval"]
        ec_comb_t = ec_lincomb(
                [
                    (W_t, 1),
                    (W_t_quot, zeta),
                    (b.G1, -t_eval),
                ]
            )

        assert b.pairing(self.X_2, W_t_quot) == b.pairing(b.G2, ec_comb_t)
        print("Done KZG10 commitment check for quotient polynomial")

        s1_eval = proof["s1_eval"]
        s2_eval = proof["s2_eval"]
        s3_eval = proof["s3_eval"]
        ql_eval = proof["ql_eval"]
        qr_eval = proof["qr_eval"]
        qm_eval = proof["qm_eval"]
        qo_eval = proof["qo_eval"]
        qc_eval = proof["qc_eval"]

        f_eval = (
            (a_eval + beta * zeta + gamma)
            * (b_eval + beta * zeta * 2 + gamma)
            * (c_eval + beta * zeta * 3 + gamma)
        )
        g_eval = (
            (a_eval + beta * s1_eval + gamma)
            * (b_eval + beta * s2_eval + gamma)
            * (c_eval + beta * s3_eval + gamma)
        )

        gate_constraints_eval = (
            ql_eval * a_eval
            + qr_eval * b_eval
            + qm_eval * a_eval * b_eval
            + qo_eval * c_eval
            + qc_eval
            + PI_ev
        )

        permutation_grand_product_eval = z_eval * f_eval - zw_eval * g_eval

        permutation_first_row_eval = L0_ev * (z_eval - 1)

        left = (
            gate_constraints_eval
            + alpha * permutation_grand_product_eval
            +  alpha ** 2 * permutation_first_row_eval
        )

        right = t_eval * ZH_ev

        assert left == right

        print("Done equation check for all constraints")
        return True

    # Compute challenges (should be same as those computed by prover)
    def compute_challenges(
        self, proof
    ) -> tuple[Scalar, Scalar, Scalar, Scalar, Scalar, Scalar]:
        transcript = Transcript(b"plonk")
        beta, gamma = transcript.round_1(proof.msg_1)
        alpha, _fft_cofactor = transcript.round_2(proof.msg_2)
        zeta = transcript.round_3(proof.msg_3)
        v = transcript.round_4(proof.msg_4)
        u = transcript.round_5(proof.msg_5)

        return beta, gamma, alpha, zeta, v, u
