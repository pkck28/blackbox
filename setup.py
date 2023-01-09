from setuptools import setup, find_packages

setup(
    name='blackbox',
    version='0.0.6',
    packages=find_packages(include=['blackbox*']),
    install_requires=["numpy>=1.21.6", "scipy>=1.7.3", "mpi4py=3.1.3", "pygeo>=1.12.2", "prefoil>=2.0.0", "pyDOE2>=1.3.0"],
)
