from setuptools import setup

setup(name='gmg',
      version='1.0',
      description='An open source geophysical modelling GUI',
      url='https://github.com/btozer/gmg.git',
      author='Brook Tozer',
      author_email='brookt@earth.ox.ac.uk',
      license='BSD-3-Clause',
      packages=['gmg'],
      install_requires=[
      'scipy',
      'matplotlib',
      'numpy',
      'fatiando',
      'obspy',
      'wxpython'
      ],
      include_package_data=True,
      zip_safe=False)
