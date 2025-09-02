# Available at setup time due to pyproject.toml
from glob import glob
from pathlib import Path
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages



__version__ = "25.9.5"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

# SRC_DIR = Path("src").resolve()

def print_sources(sources):
    print(sources)
    return sources

ext_modules = [
    Pybind11Extension("_Sailfish",
        sources = sorted(print_sources(glob("src/*.cpp"))),
        cxx_std = "17",
        extra_objects=[str(x) for x in Path(".").resolve().glob("libs/*") if x.is_file()],
        # extra_compile_args=["-g"],
        # Example: passing in the version to the compiled code
        define_macros = [('VERSION_INFO', __version__)],
        ),
]


setup(
    name="msasim",
    version=__version__,
    author="Elya Wygoda",
    author_email="elya.wygoda@gmail.com",
    url="https://github.com/elyawy/Sailfish-backend",
    description="A fast MSA simulator",
    # long_description="Sailfish is a performant multiple sequence alignment simulator, written in C++, allowing fast generation of large simualted datasets.",
    long_description=open("README.md", 'r').read(),
    long_description_content_type='text/markdown',
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.6",
    packages=find_packages(include=['msasim','tests']) + ['_Sailfish'],
    package_data={
        "_Sailfish":["py.typed","__init__.pyi"],
        },
)
