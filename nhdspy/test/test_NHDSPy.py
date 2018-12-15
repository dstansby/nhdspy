import nhdspy


def test_smoke():
    electrons = nhdspy.Species(-1, 1 / 1836, 1, 0, 1, 1)
    protons = nhdspy.Species(1, 1, 1, 0, 1, 1)
    input = nhdspy.InputParams([electrons, protons],
                               0.009 - 0.0001 * 1j,
                               0.001, 1e-4, 0.1, 1)
    assert input.total_charge == 0
    assert input.total_current == 0

    output = nhdspy.run(input)
