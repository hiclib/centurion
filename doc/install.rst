================================================================================
Installation
================================================================================

This package uses distutils, which is the default way of installing
python modules.

The dependencies are:

- python (>= 2.6)
- setuptools
- numpy (>= 1.3)
- scipy (>= 0.7)

All of these dependencies can be installed at once using `Anaconda
<http://docs.continuum.io/anaconda/install.html>`_.

To install in your home directory, use::

    python setup.py install --user

or using pip (if pip is installed)::

    pip install --user centurion

To install for all users on Unix/Linux::

    python setup.py build
    sudo python setup.py install

or using pip::

  pip install centurion

This will install a python package ``centurion``.
