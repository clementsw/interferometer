# Interferometer package

The ability to implement arbitrary interference between any number of optical inputs is useful for both classical and quantum optics. This package provides implementations of the following two papers:

- Reck, Michael, et al. "Experimental realization of any discrete unitary operator." Physical review letters 73.1 (1994): 58.
- Clements, William R., et al. "Optimal design for universal multiport interferometers." Optica 3.12 (2016): 1460-1465.

These papers show how any interfometer can be implemented by triangular or square meshes of 2x2 beam splitters, which are experimentally simple to fabricate.

This package also includes the following features:
- a function to generate Haar-random unitary matrices of any size
- a method for drawing interferometers using matplotlib
- a method for building up interferometers one beam splitter at a time

### Installation

This package requires python 3, and apart from the tests and the demo notebook only numpy and matplotlib are required. You can install this package with:

```
pip install interferometer
```

### Examples

A demonstration notebook can be found in the notebook folder. Alternatively, you can try the following code snippet that produces a 5x5 random unitary matrix U describing an interferometer, and determines how this interferometer can be implemented using a square mesh of beam splitters.

```
import interferometer as itf

U = itf.random_unitary(5)
I = itf.square_decomposition(U)
I.draw()
```

### Testing

You can run the tests from the root directory with
```
pip install pytest
python -m pytest
```

### Contributing

If you want to add new features to this package, feel free to suggest pull requests!


