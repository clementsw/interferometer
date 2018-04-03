# interferometer

Implements the algorithms in:

Reck, Michael, et al. "Experimental realization of any discrete unitary operator." Physical review letters 73.1 (1994): 58.

Clements, William R., et al. "Optimal design for universal multiport interferometers." Optica 3.12 (2016): 1460-1465.

Install:
```bash
$pip install interferometer
```
Some examples:
```import interferometer as itf

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
```