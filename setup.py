from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# codes to build myfpzip.so 
ext_modules=[
    Extension("myfpzip",
              ["myfpzip.pyx"],
              libraries=["fpzip"],
              include_dirs=["fpzip-custom-version/inc"],
              library_dirs=["fpzip-custom-version/lib"])
               
]

setup(
  name="call c fpzip",
  cmdclass={"build_ext": build_ext},
  ext_modules = ext_modules
)
