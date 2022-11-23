from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
import glob



__version__ = '0.0.1'


# src_files = (
#     'blackbody',
#     'gaussianQuadrature',
#     'initialize',
#     'io',
#     'lambda',
#     'lambdaIterate',
#     'lineProfiles',
#     'NuModel',
#     'RadModel',
#     'util',
# )

src_files = glob.glob('src/*.cpp')
src_files.extend(glob.glob('tests/test_cpp/*.cpp'))

ext_modules = [
    Pybind11Extension('Rad1D',
        src_files,
        define_macros = [('VERSION_INFO', __version__)],
    ),
]

setup(
    name='Rad1D',
    version=__version__,
    author='Anthony Burrow',
    author_email='',
    url='https://github.com/anthonyburrow/Rad1D',
    description='A simple radiative transfer code.',
    long_description='',
    ext_modules=ext_modules,
    extras_require={'test': 'pytest'},
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
    python_requires='>=3.6',
)
