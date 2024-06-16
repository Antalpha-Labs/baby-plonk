from compiler.program import Program, CommonPreprocessedInput
from utils import *
from setup import *
from typing import Optional
from dataclasses import dataclass
from transcript import Transcript, Message1, Message2, Message3, Message4, Message5
from poly import Polynomial, Basis


@dataclass
class Proof:
    msg_1: Message1
    msg_2: Message2
    msg_3: Message3
    msg_4: Message4
    msg_5: Message5

    def flatten(self):
        proof = {}
        proof["a_1"] = self.msg_1.a_1
        proof["b_1"] = self.msg_1.b_1
        proof["c_1"] = self.msg_1.c_1
        proof["z_1"] = self.msg_2.z_1
        proof["W_t"] = self.msg_3.W_t
        proof["a_eval"] = self.msg_4.a_eval
        proof["b_eval"] = self.msg_4.b_eval
        proof["c_eval"] = self.msg_4.c_eval
        proof["ql_eval"] = self.msg_4.ql_eval
        proof["qr_eval"] = self.msg_4.qr_eval
        proof["qm_eval"] = self.msg_4.qm_eval
        proof["qo_eval"] = self.msg_4.qo_eval
        proof["qc_eval"] = self.msg_4.qc_eval
        proof["s1_eval"] = self.msg_4.s1_eval
        proof["s2_eval"] = self.msg_4.s2_eval
        proof["s3_eval"] = self.msg_4.s3_eval
        proof["z_eval"] = self.msg_4.z_eval
        proof["zw_eval"] = self.msg_4.zw_eval
        proof["t_eval"] = self.msg_4.t_eval
        proof["W_a"] = self.msg_5.W_a
        proof["W_a_quot"] = self.msg_5.W_a_quot
        proof["W_b"] = self.msg_5.W_b
        proof["W_b_quot"] = self.msg_5.W_b_quot
        proof["W_c"] = self.msg_5.W_c
        proof["W_c_quot"] = self.msg_5.W_c_quot
        proof["W_ql"] = self.msg_5.W_ql
        proof["W_ql_quot"] = self.msg_5.W_ql_quot
        proof["W_qr"] = self.msg_5.W_qr
        proof["W_qr_quot"] = self.msg_5.W_qr_quot
        proof["W_qm"] = self.msg_5.W_qm
        proof["W_qm_quot"] = self.msg_5.W_qm_quot
        proof["W_qo"] = self.msg_5.W_qo
        proof["W_qo_quot"] = self.msg_5.W_qo_quot
        proof["W_qc"] = self.msg_5.W_qc
        proof["W_qc_quot"] = self.msg_5.W_qc_quot
        proof["W_s1"] = self.msg_5.W_s1
        proof["W_s1_quot"] = self.msg_5.W_s1_quot
        proof["W_s2"] = self.msg_5.W_s2
        proof["W_s2_quot"] = self.msg_5.W_s2_quot
        proof["W_s3"] = self.msg_5.W_s3
        proof["W_s3_quot"] = self.msg_5.W_s3_quot
        proof["W_z"] = self.msg_5.W_z
        proof["W_z_quot"] = self.msg_5.W_z_quot
        proof["W_zw"] = self.msg_5.W_zw
        proof["W_zw_quot"] = self.msg_5.W_zw_quot
        proof["W_t"] = self.msg_5.W_t
        proof["W_t_quot"] = self.msg_5.W_t_quot
        return proof


@dataclass
class Prover:
    group_order: int
    setup: Setup
    program: Program
    pk: CommonPreprocessedInput

    def __init__(self, setup: Setup, program: Program):
        self.group_order = program.group_order
        self.setup = setup
        self.program = program
        self.pk = program.common_preprocessed_input()

    def prove(self, witness: dict[Optional[str], int]) -> Proof:
        # Initialise Fiat-Shamir transcript
        transcript = Transcript(b"plonk")

        # Collect fixed and public information
        # FIXME: Hash pk and PI into transcript
        public_vars = self.program.get_public_assignments()
        # Public input polynomial
        PI = Polynomial(
            [Scalar(-witness[v]) for v in public_vars]
            + [Scalar(0) for _ in range(self.group_order - len(public_vars))],
            Basis.LAGRANGE,
        )
        self.PI = PI

        # Round 1
        msg_1 = self.round_1(witness)
        self.beta, self.gamma = transcript.round_1(msg_1)

        # Round 2
        msg_2 = self.round_2()
        self.alpha = transcript.round_2(msg_2)

        # Round 3
        msg_3 = self.round_3()
        self.zeta = transcript.round_3(msg_3)

        # Round 4
        msg_4 = self.round_4()

        # Round 5
        msg_5 = self.round_5()

        return Proof(msg_1, msg_2, msg_3, msg_4, msg_5)

    def round_1(
        self,
        witness: dict[Optional[str], int],
    ) -> Message1:
        # https://github.com/sec-bit/learning-zkp/blob/master/plonk-intro-cn/1-plonk-arithmetization.md
        program = self.program
        setup = self.setup
        group_order = self.group_order

        if None not in witness:
            witness[None] = 0

        # Compute wire assignments
        A_values = [Scalar(0) for _ in range(group_order)]
        B_values = [Scalar(0) for _ in range(group_order)]
        C_values = [Scalar(0) for _ in range(group_order)]

        for i, gate_wires in enumerate(program.wires()):
            A_values[i] = Scalar(witness[gate_wires.L])
            B_values[i] = Scalar(witness[gate_wires.R])
            C_values[i] = Scalar(witness[gate_wires.O])

        self.A = Polynomial(A_values, Basis.LAGRANGE)
        self.B = Polynomial(B_values, Basis.LAGRANGE)
        self.C = Polynomial(C_values, Basis.LAGRANGE)

        a_1 = setup.commit(self.A)
        b_1 = setup.commit(self.B)
        c_1 = setup.commit(self.C)

        # Sanity check that witness fulfils gate constraints
        assert (
            self.A * self.pk.QL
            + self.B * self.pk.QR
            + self.A * self.B * self.pk.QM
            + self.C * self.pk.QO
            + self.PI
            + self.pk.QC
            == Polynomial([Scalar(0)] * group_order, Basis.LAGRANGE)
        )

        return Message1(a_1, b_1, c_1)

    def round_2(self) -> Message2:
        # https://github.com/sec-bit/learning-zkp/blob/master/plonk-intro-cn/3-plonk-permutation.md
        group_order = self.group_order
        setup = self.setup

        Z_values = [Scalar(1)]
        roots_of_unity = Scalar.roots_of_unity(group_order)
        for i in range(group_order):
            Z_values.append(
                Z_values[-1]
                * self.rlc(self.A.values[i], roots_of_unity[i])
                * self.rlc(self.B.values[i], 2 * roots_of_unity[i])
                * self.rlc(self.C.values[i], 3 * roots_of_unity[i])
                / self.rlc(self.A.values[i], self.pk.S1.values[i])
                / self.rlc(self.B.values[i], self.pk.S2.values[i])
                / self.rlc(self.C.values[i], self.pk.S3.values[i])
            )
        # The last value is 1
        assert Z_values.pop() == 1

        # Sanity-check that Z was computed correctly
        for i in range(group_order):
            assert (
                self.rlc(self.A.values[i], roots_of_unity[i])
                * self.rlc(self.B.values[i], 2 * roots_of_unity[i])
                * self.rlc(self.C.values[i], 3 * roots_of_unity[i])
            ) * Z_values[i] - (
                self.rlc(self.A.values[i], self.pk.S1.values[i])
                * self.rlc(self.B.values[i], self.pk.S2.values[i])
                * self.rlc(self.C.values[i], self.pk.S3.values[i])
            ) * Z_values[
                (i + 1) % group_order
            ] == 0

        Z = Polynomial(Z_values, Basis.LAGRANGE)
        z_1 = setup.commit(Z)
        print("Permutation accumulator polynomial successfully generated")

        self.Z = Z
        return Message2(z_1)

    def round_3(self) -> Message3:
        # https://github.com/sec-bit/learning-zkp/blob/master/plonk-intro-cn/4-plonk-constraints.md
        group_order = self.group_order
        setup = self.setup

        # Compute the quotient polynomial

        alpha = self.alpha

        roots_of_unity = Scalar.roots_of_unity(group_order)

        A_coeff, B_coeff, C_coeff, S1_coeff, S2_coeff, S3_coeff, Z_coeff, QL_coeff, QR_coeff, QM_coeff, QO_coeff, QC_coeff, PI_coeff = (
            x.ifft()
            for x in (
                self.A,
                self.B,
                self.C,
                self.pk.S1,
                self.pk.S2,
                self.pk.S3,
                self.Z,
                self.pk.QL,
                self.pk.QR,
                self.pk.QM,
                self.pk.QO,
                self.pk.QC,
                self.PI,
            )
        )

        L0_coeff = (
            Polynomial([Scalar(1)] + [Scalar(0)] * (group_order - 1), Basis.LAGRANGE)
        ).ifft()

        # x^8 - 1 coeffs are [-1, 0, 0, 0, 0, 0, 0, 0, 1]
        # which needs 9 points(n + 1) to determine the polynomial
        ZH_array = [Scalar(-1)] + [Scalar(0)] * (group_order - 1) + [Scalar(1)]
        ZH_coeff = Polynomial(ZH_array, Basis.MONOMIAL)

        gate_constraints_coeff = (
            A_coeff * QL_coeff
            + B_coeff * QR_coeff
            + A_coeff * B_coeff * QM_coeff
            + C_coeff * QO_coeff
            + PI_coeff
            + QC_coeff
        )

        normal_roots = Polynomial(
            roots_of_unity, Basis.LAGRANGE
        )

        roots_coeff = normal_roots.ifft()
        # z * w
        ZW = self.Z.shift(1)
        ZW_coeff = ZW.ifft()

        for i in range(group_order):
            assert (
                self.rlc(self.A.values[i], roots_of_unity[i])
                * self.rlc(self.B.values[i], 2 * roots_of_unity[i])
                * self.rlc(self.C.values[i], 3 * roots_of_unity[i])
            ) * self.Z.values[i] - (
                self.rlc(self.A.values[i], self.pk.S1.values[i])
                * self.rlc(self.B.values[i], self.pk.S2.values[i])
                * self.rlc(self.C.values[i], self.pk.S3.values[i])
            ) * ZW.values[
                i % group_order
            ] == 0

        permutation_grand_product_coeff = (
            (
                self.rlc(A_coeff, roots_coeff)
                * self.rlc(B_coeff, roots_coeff * Scalar(2))
                * self.rlc(C_coeff, roots_coeff * Scalar(3))
            )
            * Z_coeff
            - (
                self.rlc(A_coeff, S1_coeff)
                * self.rlc(B_coeff, S2_coeff)
                * self.rlc(C_coeff, S3_coeff)
            )
            * ZW_coeff
        )

        permutation_first_row_coeff = (Z_coeff - Scalar(1)) * L0_coeff

        all_constraints = (
            gate_constraints_coeff
            + permutation_grand_product_coeff * alpha
            + permutation_first_row_coeff * alpha**2
        )

        # quotient polynomial
        T_coeff = all_constraints / ZH_coeff

        print("Generated the quotient polynomial")

        W_t = setup.commit(T_coeff)

        self.A_coeff = A_coeff
        self.B_coeff = B_coeff
        self.C_coeff = C_coeff
        self.S1_coeff = S1_coeff
        self.S2_coeff = S2_coeff
        self.S3_coeff = S3_coeff
        self.Z_coeff = Z_coeff
        self.ZW_coeff = ZW_coeff
        self.QL_coeff = QL_coeff
        self.QR_coeff = QR_coeff
        self.QM_coeff = QM_coeff
        self.QO_coeff = QO_coeff
        self.QC_coeff = QC_coeff
        self.PI_coeff = PI_coeff
        self.T_coeff = T_coeff

        return Message3(W_t)

    def round_4(self) -> Message4:
        # https://github.com/sec-bit/learning-zkp/blob/master/plonk-intro-cn/4-plonk-constraints.md
        group_order = self.group_order
        zeta = self.zeta

        a_eval = self.A_coeff.coeff_eval(zeta)
        b_eval = self.B_coeff.coeff_eval(zeta)
        c_eval = self.C_coeff.coeff_eval(zeta)
        s1_eval = self.S1_coeff.coeff_eval(zeta)
        s2_eval = self.S2_coeff.coeff_eval(zeta)
        s3_eval = self.S3_coeff.coeff_eval(zeta)
        root_of_unity = Scalar.root_of_unity(group_order)
        z_eval = self.Z_coeff.coeff_eval(zeta)
        zw_eval = self.Z_coeff.coeff_eval(zeta * root_of_unity)
        ql_eval = self.QL_coeff.coeff_eval(zeta)
        qr_eval = self.QR_coeff.coeff_eval(zeta)
        qm_eval = self.QM_coeff.coeff_eval(zeta)
        qo_eval = self.QO_coeff.coeff_eval(zeta)
        qc_eval = self.QC_coeff.coeff_eval(zeta)
        t_eval = self.T_coeff.coeff_eval(zeta)

        self.a_eval = a_eval
        self.b_eval = b_eval
        self.c_eval = c_eval
        self.ql_eval = ql_eval
        self.qr_eval = qr_eval
        self.qm_eval = qm_eval
        self.qo_eval = qo_eval
        self.qc_eval = qc_eval
        self.s1_eval = s1_eval
        self.s2_eval = s2_eval
        self.s3_eval = s3_eval
        self.z_eval = z_eval
        self.zw_eval = zw_eval
        self.t_eval = t_eval

        return Message4(
            a_eval,
            b_eval,
            c_eval,
            ql_eval,
            qr_eval,
            qm_eval,
            qo_eval,
            qc_eval,
            s1_eval,
            s2_eval,
            s3_eval,
            z_eval,
            zw_eval,
            t_eval
        )

    def round_5(self) -> Message5:
        W_a, W_a_quot = self.generate_commitment(self.A_coeff, self.a_eval)
        W_b, W_b_quot = self.generate_commitment(self.B_coeff, self.b_eval)
        W_c, W_c_quot = self.generate_commitment(self.C_coeff, self.c_eval)
        W_ql, W_ql_quot = self.generate_commitment(self.QL_coeff, self.ql_eval)
        W_qr, W_qr_quot = self.generate_commitment(self.QR_coeff, self.qr_eval)
        W_qm, W_qm_quot = self.generate_commitment(self.QM_coeff, self.qm_eval)
        W_qo, W_qo_quot = self.generate_commitment(self.QO_coeff, self.qo_eval)
        W_qc, W_qc_quot = self.generate_commitment(self.QC_coeff, self.qc_eval)
        W_s1, W_s1_quot = self.generate_commitment(self.S1_coeff, self.s1_eval)
        W_s2, W_s2_quot = self.generate_commitment(self.S2_coeff, self.s2_eval)
        W_s3, W_s3_quot = self.generate_commitment(self.S3_coeff, self.s3_eval)
        W_z, W_z_quot = self.generate_commitment(self.Z_coeff, self.z_eval)
        W_zw, W_zw_quot = self.generate_commitment(self.ZW_coeff, self.zw_eval)
        W_t, W_t_quot = self.generate_commitment(self.T_coeff, self.t_eval)

        print("Generated final quotient witness polynomials")
        return Message5(
            W_a, W_a_quot,
            W_b, W_b_quot,
            W_c, W_c_quot,
            W_ql, W_ql_quot,
            W_qr, W_qr_quot,
            W_qm, W_qm_quot,
            W_qo, W_qo_quot,
            W_qc, W_qc_quot,
            W_s1, W_s1_quot,
            W_s2, W_s2_quot,
            W_s3, W_s3_quot,
            W_z, W_z_quot,
            W_zw, W_zw_quot,
            W_t, W_t_quot,
        )

    def rlc(self, term_1, term_2):
        return term_1 + term_2 * self.beta + self.gamma

    def generate_commitment(self, coeff: Polynomial, eval: Scalar):
        setup = self.setup
        zeta = self.zeta
        # Polynomial for (X - zeta)
        ZH_zeta_coeff = Polynomial([-zeta, Scalar(1)], Basis.MONOMIAL)
        quot_coeff = (coeff - eval) / ZH_zeta_coeff
        # witness for polynomial itself
        w = setup.commit(coeff)
        # witness for quotient polynomial
        w_quot = setup.commit(quot_coeff)
        return w, w_quot
