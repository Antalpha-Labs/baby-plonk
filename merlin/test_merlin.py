from merlin.merlin_transcript import MerlinTranscript
from merlin.strobe import Strobe128


def test_strobe_conformance():
    s: Strobe128 = Strobe128.new(b"Conformance Test Protocol")

    msg = int(99).to_bytes(1, "big") * 1024

    s.meta_ad(b"ms", False)
    s.meta_ad(b"g", True)
    s.ad(msg, False)

    s.meta_ad(b"prf", False)
    prf = s.prf(32, False)
    print("PRF: ", prf.hex())
    assert (
        prf.hex() == "b48e645ca17c667fd5206ba57a6a228d72d8e1903814d3f17f622996d7cfefb0"
    )

    s.meta_ad(b"key", False)
    s.key(prf, False)

    s.meta_ad(b"prf", False)
    prf = s.prf(32, False)

    print("PRF: ", prf.hex())
    assert (
        prf.hex() == "07e45cce8078cee259e3e375bb85d75610e2d1e1201c5f645045a194edd49ff8"
    )


def test_merlin():
    transcript = MerlinTranscript(b"test protocol")
    transcript.append_message(b"some label", b"some data")
    challenge = transcript.challenge_bytes(b"challenge", 32)
    print("challenge", challenge.hex())
    assert (
        challenge.hex()
        == "d5a21972d0d5fe320c0d263fac7fffb8145aa640af6e9bca177c03c7efcf0615"
    )
