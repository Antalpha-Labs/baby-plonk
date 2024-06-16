from compiler.program import Program
from setup import Setup
from prover import Prover
from test.mini_poseidon import rc, mds, poseidon_hash
from utils import *
import random

# Generate a random integer between a specified range by user
tau = random.randint(0, 100)
print("Random number tau: ", tau)

def prover_test():
    print("Beginning prover test")
    # powers should be 2^n so that we can use roots of unity for FFT
    # and should be bigger than len(coeffs) of polynomial to do KZG commitment
    # the value here is: powers = 4 * group_order
    # which is bigger than the order of quotient polynomial
    group_order = 8
    powers = group_order * 4
    setup = Setup.generate_srs(powers, tau)

    program = Program(["e public", "c <== a * b", "e <== c * d"], group_order)
    assignments = {"a": 3, "b": 4, "c": 12, "d": 5, "e": 60}
    prover = Prover(setup, program)
    proof = prover.prove(assignments)
    print("Prover test success")
    return setup, proof, group_order

def verifier_test(setup, proof, group_order):
    print("Beginning verifier test")
    program = Program(["e public", "c <== a * b", "e <== c * d"], group_order)
    public = [60]
    vk = setup.verification_key(program.common_preprocessed_input())
    assert vk.verify_proof(group_order, proof, public)
    print("Verifier test success")


def factorization_test():
    print("Beginning test: prove you know small integers that multiply to 91")
    group_order = 16
    powers = group_order * 4
    setup = Setup.generate_srs(powers, tau)

    program = Program.from_str(
        """n public
        pb0 === pb0 * pb0
        pb1 === pb1 * pb1
        pb2 === pb2 * pb2
        pb3 === pb3 * pb3
        qb0 === qb0 * qb0
        qb1 === qb1 * qb1
        qb2 === qb2 * qb2
        qb3 === qb3 * qb3
        pb01 <== pb0 + 2 * pb1
        pb012 <== pb01 + 4 * pb2
        p <== pb012 + 8 * pb3
        qb01 <== qb0 + 2 * qb1
        qb012 <== qb01 + 4 * qb2
        q <== qb012 + 8 * qb3
        n <== p * q""",
        16,
    )
    public = [91]
    vk = setup.verification_key(program.common_preprocessed_input())
    print("Generated verification key")
    assignments = program.fill_variable_assignments(
        {
            "pb3": 1,
            "pb2": 1,
            "pb1": 0,
            "pb0": 1,
            "qb3": 0,
            "qb2": 1,
            "qb1": 1,
            "qb0": 1,
        }
    )
    prover = Prover(setup, program)
    proof = prover.prove(assignments)
    print("Generated proof")
    assert vk.verify_proof(group_order, proof, public)
    print("Factorization test success!")

def output_proof_lang() -> str:
    o = []
    o.append("L0 public")
    o.append("M0 public")
    o.append("M64 public")
    o.append("R0 <== 0")
    for i in range(64):
        for j, pos in enumerate(("L", "M", "R")):
            f = {"x": i, "r": rc[i][j], "p": pos}
            if i < 4 or i >= 60 or pos == "L":
                o.append("{p}adj{x} <== {p}{x} + {r}".format(**f))
                o.append("{p}sq{x} <== {p}adj{x} * {p}adj{x}".format(**f))
                o.append("{p}qd{x} <== {p}sq{x} * {p}sq{x}".format(**f))
                o.append("{p}qn{x} <== {p}qd{x} * {p}adj{x}".format(**f))
            else:
                o.append("{p}qn{x} <== {p}{x} + {r}".format(**f))
        for j, pos in enumerate(("L", "M", "R")):
            f = {"x": i, "p": pos, "m": mds[j]}
            o.append("{p}suma{x} <== Lqn{x} * {m}".format(**f))
            f = {"x": i, "p": pos, "m": mds[j + 1]}
            o.append("{p}sumb{x} <== {p}suma{x} + Mqn{x} * {m}".format(**f))
            f = {"x": i, "xp1": i + 1, "p": pos, "m": mds[j + 2]}
            o.append("{p}{xp1} <== {p}sumb{x} + Rqn{x} * {m}".format(**f))
    return "\n".join(o)


def poseidon_test():
    # PLONK-prove the correctness of a Poseidon execution. Note that this is
    # a very suboptimal way to do it: an optimized implementation would use
    # a custom PLONK gate to do a round in a single gate
    expected_value = poseidon_hash(1, 2)
    group_order = 1024
    powers = group_order * 4
    setup = Setup.generate_srs(powers, tau)
    # Generate code for proof
    program = Program.from_str(output_proof_lang(), group_order)
    print("Generated code for Poseidon test")
    assignments = program.fill_variable_assignments({"L0": 1, "M0": 2})
    vk = setup.verification_key(program.common_preprocessed_input())
    print("Generated verification key")
    prover = Prover(setup, program)
    proof = prover.prove(assignments)
    print("Generated proof")
    assert vk.verify_proof(group_order, proof, [1, 2, expected_value])
    print("Verified proof!")

if __name__ == "__main__":
    setup, proof, group_order = prover_test()
    verifier_test(setup, proof, group_order)
    # comment out them if you need to test them
    # factorization_test()
    # poseidon_test()
