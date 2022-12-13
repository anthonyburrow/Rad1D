# Rad1D

A simple 1D, plane-parallel radiative transfer code that implements accelerated
lambda iteration (ALI).

This project was made in a collaboration including Manuel Barrientos, Anthony
Burrow, Adam Moss, and Sarah Stangl.

For more details on this project, see `.doc/FinalProjectReport.pdf`.

## Requirements

This may be compiled and run on either Windows or Linux platforms with an
available C++ compiler. To use Rad1D as a Python-wrapped module, the following
are required:

* Python 3.6+
* `pybind11` package

## Installing the code

Rad1D may be installed with `pip` (setuptools) by:

```
git clone https://github.com/anthonyburrow/Rad1D.git
pip install ./Rad1D/
```

## Running with Python

The `RadModel` class is exposed to Python and instantiated in the following
way:

```python
from Rad1D import RadModel

params = {
    # Parameter key-values here
}

model = RadModel(params)
```

See "./examples/example1.py" for how `RadModel` may be used. Remember to set
all absolute paths in the global parameters with the correct paths for the user
machine.

## Parameters

Below are the currently available global parameters, and their default values:
```
params = {
    'data_dir'       : '../data',   # Absolute path to line list
    'wave_start'     : 4000.,       # Starting wavelength
    'wave_end'       : 7000.,       # Ending wavelength
    'cont_res'       : 0.05,        # Points per angstrom resolved for continuum
    'line_res'       : 3.0,         # Points per angstrom resolved for lines in line list
    'tau_min'        : 1e-8,        # Minimum tau of atmosphere for continuum
    'tau_max'        : 1e6,         # Maximum tau of atmosphere for continuum
    'eps'            : 1e-4,        # Thermalization factor
    'T_eff'          : 6000.,       # Characteristic temperature (K) of atmosphere
    'n_zones'        : 256,         # Number of tau points
    'max_iter'       : 500,         # Maximum number of lambda iterations allowed
    'accelerated'    : True,        # Use accelerated lambda iteration
    'Ng_accelerated' : True,        # Implement Ng acceleration
    'eps_converge'   : 1e-8,        # Factor to determine J is converged
    'n_quad'         : 32,          # Order of Gaussian quadrature integration
    'verbose'        : True,        # Display stdout output
}
```

# Running the Tests

One may run the suite of tests that generates plots in the documentation with
```
pytest -vs ./Rad1D/tests
```
(Point to the tests directory accordingly.)
