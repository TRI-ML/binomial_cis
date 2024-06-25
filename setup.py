from setuptools import setup, find_packages

setup(
   name="binomial_cis",
   version='0.0.10',
   author="Joe Vincent",
   description="Confidence intervals for binomial distributions.",
   packages=find_packages(),
   install_requires=['numpy', 'numba', 'scipy', 'matplotlib'],
   long_description=open('README.md').read(),
   long_description_content_type='text/markdown'
)
