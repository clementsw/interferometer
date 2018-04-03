This package contains code for manipulating universal multiport interferometers,
and in particular implements the decomposition algorithm in:
Clements, William R., et al. "Optimal design for universal multiport interferometers." Optica 3.12 (2016): 1460-1465

It introduces an Interferometer class containing a list of beam splitters and phase shifters,
and contains algorithms for decomposing any unitary matrix into an interferometer. This code
also contains a method for viewing an interferometer

Installation ---------

For a local installation, use

$ pip install .

Example 1 ------

import interferometer as itf

I = itf.Interferometer()
I.add_BS([1,2,1,4])
I.add_BS([2,3,3,2])
I.add_BS([3,4,8,9])
I.add_BS([1,2,2,2])
I.add_phase([3,2])
I.draw_interferometer()

Example 2 ------

U = itf.random_unitary(5)
I = itf.square_decomposition(U)
I.draw_interferometer()