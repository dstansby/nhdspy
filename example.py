import nhdspy
import matplotlib.pyplot as plt

electrons = nhdspy.Species(-1, 1 / 1836, 1, 0, 1, 1)
protons = nhdspy.Species(1, 1, 1, 0, 1, 1)

va = 1e-4
theta_kB = 0.001
omega_0 = 0.009 - 0.0001 * 1j
kzmin = 0.01
kzmax = 0.2
input = nhdspy.InputParams([electrons, protons], omega_0, theta_kB,
                           va, kzmin, kzmax)
output = nhdspy.run(input)

fig, ax = plt.subplots()
ax.plot(output.kz, output.omega_real, label='Real part')
ax.plot(output.kz, output.omega_imag, label='Imaginary part')
ax.set_xlabel('$kd_{p}$')
ax.set_ylabel(r'$\omega / \Omega_{p}$')
ax.legend()

plt.show()
