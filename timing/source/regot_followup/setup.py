from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension,build_ext
from glob import glob
import os
setup(name="regot_followup",version="0.1",ext_modules=[Pybind11Extension("regot_followup",["adapter.cpp"]+sorted(glob("src/*.cpp")),include_dirs=["src",os.environ["EIGEN3_INCLUDE_DIR"]],extra_compile_args=["-O3","-mavx2","-fopenmp","-DNDEBUG"],extra_link_args=["-fopenmp"])],cmdclass={"build_ext":build_ext})
