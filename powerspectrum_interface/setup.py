from setuptools import setup
from distutils.command.install import install as _install
from distutils.command.clean import clean as _clean
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop

import subprocess

import os
import shutil

LIBS = ["libhmx_wrapper.so",]
HERE = os.path.abspath(os.path.dirname(__file__))


def compile_library(env):
    subprocess.check_call(["make", "debug"], env=env, cwd="fortran_wrapper")
    shutil.copy(os.path.join(HERE, "fortran_wrapper/lib/libhmx_wrapper.so"),
                os.path.join(HERE, "pyhmx/libhmx_wrapper.so"))


def clean_library(env={}):
    subprocess.check_call(["make", "clean"], env=env, cwd="fortran_wrapper")
    try:
        os.remove(os.path.join(HERE, "pyhmx/libhmx_wrapper.so"))
    except FileNotFoundError:
        pass


class build(_build_py):
    def run(self):
        env = os.environ
        compile_library(env)
        super().run()


class develop(_develop):
    def run(self):
        env = os.environ
        compile_library(env)
        super().run()


class install(_install):
    def __init__(self, dist):
        super().__init__(dist)
        self.build_args = {}
        if self.record is None:
            self.record = "install-record.txt"

    def run(self):
        super().run()


class clean(_clean):
    def run(self):
        clean_library()
        super().run()


setup(name              = "pyhmx",
      description       = "Python interface for the HMx library",
      author            = "Tilman Troester",
      author_email      = "tilmantroester@gmail.com",
      packages          = ["pyhmx"],
      package_data      = {"" : LIBS,},
      install_requires  = ['numpy'],
      cmdclass          = {"install"   : install,
                           "develop"   : develop,
                           "build_py"  : build,
                           "clean"     : clean},
)
