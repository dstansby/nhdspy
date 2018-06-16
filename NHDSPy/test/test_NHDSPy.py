import NHDSPy


def test_smoke():
    electrons = NHDSPy.Species(-1, 1 / 1836, 1, 0, 1, 1)
    protons = NHDSPy.Species(1, 1, 1, 0, 1, 1)
    input = NHDSPy.InputParams([electrons, protons], 0.001, 1e-4)
    output = NHDSPy.run(input)
