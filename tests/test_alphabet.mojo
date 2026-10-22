from std.testing import assert_equal, TestSuite

from plegma import N_CODE, CodeBuffer, encode_base


def test_encode_base() raises:
    var text = String("ACGTUacgtuNn-")
    var want: List[UInt8] = [0, 1, 2, 3, 3, 0, 1, 2, 3, 3, N_CODE, N_CODE, N_CODE]
    var bytes = text.as_bytes()
    for i in range(len(bytes)):
        assert_equal(encode_base(bytes[i]), want[i])


def test_code_buffer_matches_scalar() raises:
    # Longer than one SIMD block, with a scalar tail, so both paths of encode_into run.
    var text = String()
    for _ in range(5):
        text += "ACGTNacgtnRYKM*.TTGCA"
    var bytes = text.as_bytes()
    var buf = CodeBuffer.encode(bytes, padding=8)
    assert_equal(len(buf), len(bytes))
    for i in range(len(bytes)):
        assert_equal(buf[i], encode_base(bytes[i]))
    for i in range(len(bytes), len(bytes) + 8):
        assert_equal(buf[i], N_CODE)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
