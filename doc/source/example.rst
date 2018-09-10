Examples
========

.. code-block:: python

  import nhdspy
  import matplotlib.pyplot as plt

  electrons = nhdspy.Species(-1, 1 / 1836, 1, 0, 1, 1)
  protons = nhdspy.Species(1, 1, 1, 0, 1, 1)
  propagation_angle = 0.001
  input = nhdspy.InputParams([protons, electrons], propagation_angle, 1e-4)
  print(nhdspy.format_input_file(input))
  output = nhdspy.run(input)

  fig, ax = plt.subplots()
  ax.plot(output.kz, output.omega_real)
  ax.plot(output.kz, output.omega_imag)
  plt.show()
