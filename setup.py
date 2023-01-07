from setuptools import setup, find_packages

setup(
    name='blackbox',
    version='0.0.6',
    packages=find_packages(include=['blackbox*']),
    install_requires=["numpy", "scipy"],
)
