# Rad1D

A simple 1D radiative transfer code.

## Requirements

Currently this may be compiled and run on either Windows or Linux platforms
with a. available C++ compiler. To use Rad1D as a Python-wrapped module, the
following are needed:

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

See "./tests/example.py" for a list of input parameters and how `RadModel` may
be used.
