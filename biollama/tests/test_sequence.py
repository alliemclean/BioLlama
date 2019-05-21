import sys
import biollama.tools.sequence as seq


def test_sequence():
    print(seq.Sequence())


if __name__ == "__main__":
    if len(sys.argv) > 1:
        test_sequence()
