from setuptools import setup, find_packages

setup(
   name="binomial_cis",
   version='0.0.8',
   author="Joe Vincent",
   description="Confidence intervals for binomial distributions.",
   packages=find_packages(),
   install_requires=['numpy', 'numba', 'scipy', 'matplotlib']
)