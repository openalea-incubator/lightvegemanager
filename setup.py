import sys
from setuptools import setup, find_packages

if sys.version_info < (3, 0):
    print('ERROR: requires at least Python 3.0 to run.')
    sys.exit(1)

setup(
    name="LightVegeManager",
    version="0.0.0",
    description="Works only with cn-wheat",
    url="https://github.com/mwoussen/lightvegemanager",
    packages=find_packages(include=["src","runscript"]),
)
