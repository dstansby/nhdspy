import NHDSPy
import matplotlib.pyplot as plt

input = NHDSPy.InputParams([], 0, 1e-4)
print(NHDSPy.format_input_file(input))
output = NHDSPy.run(input)

fig, ax = plt.subplots()
ax.plot(output.kz, output.x_real)
ax.plot(output.kz, output.x_imag)
plt.show()
