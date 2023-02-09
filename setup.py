import sys
from setuptools import setup, find_packages

import src, runscripts, PyRATP

if sys.version_info < (3, 0):
    print('ERROR: requires at least Python 3.0 to run.')
    sys.exit(1)

setup(
    name="LightVegeManager",
    version="0.0.0",
    description="Manage light with several plant models",
    url="https://github.com/mwoussen/lightvegemanager/",
    packages=find_packages(include=["PyRATP","src","runscript"]),
    package_data = {'' : ['*.pyd', '*.so'],}
)
