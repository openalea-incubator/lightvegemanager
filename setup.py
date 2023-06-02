import sys
from setuptools import setup, find_packages

if sys.version_info < (3, 0):
    print('ERROR: requires at least Python 3.0 to run.')
    sys.exit(1)

setup(
    name="LightVegeManager",
    version="1.0.0",
    description="Light management tool for plant modelling with tools from OpenAlea platform",
    url="https://github.com/mwoussen/lightvegemanager",

    packages=find_packages('lightvegemanager'),
    package_dir={'' : 'src'}
)
