import py_ecc.bn128 as b
from utils import *
from dataclasses import dataclass
from curve import *
from transcript import Transcript
from poly import Polynomial, Basis


@dataclass
class VerificationKey:
    # https://github.com/sec-bit/learning-zkp/blob/develop/plonk-intro-cn/plonk-constraints.md
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
        # Compute challenges
        beta, gamma, alpha, zeta = self.compute_challenges(pf)
        proof = pf.flatten()

        # Compute zero polynomial evaluation Z_H(ζ) = ζ^n - 1
        ZH_ev = zeta**group_order - 1

        # Compute Lagrange polynomial evaluation L_0(ζ)
        L0_ev = ZH_ev / (group_order * (zeta - 1))

        # Compute public input polynomial evaluation PI(ζ).
        PI = Polynomial(
            [Scalar(-x) for x in public]
            + [Scalar(0) for _ in range(group_order - len(public))],
            Basis.LAGRANGE,
        )
        PI_ev = PI.barycentric_eval(zeta)

        # Verify KZG10 commitment
        self.verify_commitment(proof, proof["W_a"], "W_a_quot", "a_eval", zeta)
        self.verify_commitment(proof, proof["W_b"], "W_b_quot", "b_eval", zeta)
        self.verify_commitment(proof, proof["W_c"], "W_c_quot", "c_eval", zeta)
        self.verify_commitment(proof, proof["W_z"], "W_z_quot", "z_eval", zeta)
        self.verify_commitment(proof, proof["W_zw"], "W_zw_quot", "zw_eval", zeta)
        self.verify_commitment(proof, proof["W_t"], "W_t_quot", "t_eval", zeta)
        self.verify_commitment(proof, self.Ql, "W_ql_quot", "ql_eval", zeta)
        self.verify_commitment(proof, self.Qr, "W_qr_quot", "qr_eval", zeta)
        self.verify_commitment(proof, self.Qm, "W_qm_quot", "qm_eval", zeta)
        self.verify_commitment(proof, self.Qo, "W_qo_quot", "qo_eval", zeta)
        self.verify_commitment(proof, self.Qc, "W_qc_quot", "qc_eval", zeta)
        self.verify_commitment(proof, self.S1, "W_s1_quot", "s1_eval", zeta)
        self.verify_commitment(proof, self.S2, "W_s2_quot", "s2_eval", zeta)
        self.verify_commitment(proof, self.S3, "W_s3_quot", "s3_eval", zeta)

        # Verify constraints and permutation
        a_eval = proof["a_eval"]
        b_eval = proof["b_eval"]
        c_eval = proof["c_eval"]
        ql_eval = proof["ql_eval"]
        qr_eval = proof["qr_eval"]
        qm_eval = proof["qm_eval"]
        qo_eval = proof["qo_eval"]
        qc_eval = proof["qc_eval"]
        s1_eval = proof["s1_eval"]
        s2_eval = proof["s2_eval"]
        s3_eval = proof["s3_eval"]
        z_eval = proof["z_eval"]
        zw_eval = proof["zw_eval"]
        t_eval = proof["t_eval"]

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
    ) -> tuple[Scalar, Scalar, Scalar, Scalar]:
        transcript = Transcript(b"plonk")
        beta, gamma = transcript.round_1(proof.msg_1)
        alpha = transcript.round_2(proof.msg_2)
        zeta = transcript.round_3(proof.msg_3)

        return beta, gamma, alpha, zeta

    def verify_commitment(self, proof, W, W_quot_key, eval_key, zeta):
        W_quot = proof[W_quot_key]
        eval = proof[eval_key]
        ec_comb = ec_lincomb(
            [
                (W, 1),
                (W_quot, zeta),
                (b.G1, -eval),
            ]
        )

        assert b.pairing(self.X_2, W_quot) == b.pairing(b.G2, ec_comb)
        print(f"Done KZG10 commitment check for {eval_key} polynomial")
