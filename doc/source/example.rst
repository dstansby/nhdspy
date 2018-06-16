Examples
========

.. code-block:: python

  import NHDSPy
  import matplotlib.pyplot as plt

  electrons = NHDSPy.Species(-1, 1 / 1836, 1, 0, 1, 1)
  protons = NHDSPy.Species(1, 1, 1, 0, 1, 1)
  propagation_angle = 0.001
  input = NHDSPy.InputParams([electrons, protons], propagation_angle, 1e-4)
  print(NHDSPy.format_input_file(input))
  output = NHDSPy.run(input)

  fig, ax = plt.subplots()
  ax.plot(output.kz, output.omega_real)
  ax.plot(output.kz, output.omega_imag)
  plt.show()
