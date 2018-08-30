import nhdspy
import matplotlib.pyplot as plt

electrons = nhdspy.Species(-1, 1 / 1836, 1, 0, 1, 1)
protons = nhdspy.Species(1, 1, 1, 0, 1, 1)
input = nhdspy.InputParams([electrons, protons], 0.009 - 0.0001 * 1j, 0.001, 1e-4)
output = nhdspy.run(input)

fig, ax = plt.subplots()
ax.plot(output.kz, output.omega_real)
ax.plot(output.kz, output.omega_imag)
plt.show()
