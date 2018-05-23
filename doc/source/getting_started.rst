.. _getting-started:

================================
Getting started with PhyreEngine
================================

.. _installing-via-conda:

Installing via :command:`conda`
===============================

The recommended way to use PhyreEngine is via `Conda <https://conda.io/>`_, a
package manager used for managing software *environments*. A conda environment
is just a directory containing the binaries, libraries, and data files
installed by various software *packages*, along with a small amount of metadata
used internally by conda. The big advantage of conda is that it can be used to
install packages in any language, so it can also be use to manage the various
tools used by PhyreEngine pipelines.

For details about using conda, see the `Conda user guide
<https://conda.io/docs/user-guide/index.html>`_. Very briefly, conda can be
installed and set up with the following commands:

.. code-block:: bash

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -p $CONDA -b
    conda create -n phyreengine

This will install conda to the directory ``$CONDA``, then create an empty
environment called ``phyreengine``. This is where we will install PhyreEngine
and any dependencies.

Conda retrieves packages from online *channels*---most package managers would
call these *repositories*. We want to use the following channels:

.. TODO: Edit sbg-conda URLs when we actually make a release.

`sbg-conda <https://github.com/StefansM/>`_
    Software packaged by (not necessarily developed by!) Imperial College
    London's `Structural Bioinformatics Group <www.sbg.bio.ic.ac.uk>`_.
    This includes PhyreEngine, and the immediate dependencies of PhyreEngine.

`defaults <https://anaconda.org/>`_
    The default package repository developed by Anaconda. This includes
    system tools such as compilers or core utilities. Most systems will
    already include the tools we install from here, but using the Anaconda
    packages allows us to target single versions of these tools.

`bioconda <https://bioconda.github.io/>`_
    A channel containing common bioinformatics software.

`conda-forge <https://conda-forge.org/>`_
    A community-run channel containing many different applications.

Configure these channels for the ``phyreengine`` environment by placing the
following in ``$CONDA/envs/phyreengine/.condarc``:

.. code-block:: yaml

    channels:
    - file:///bmm/phyreengine/software/conda
    - defaults
    - bioconda
    - conda-forge

Finally, activate the ``phryeengine`` environment and install PhyreEngine.
"Activating" a conda environment alters your ``$PATH`` and various other
environment variables, so that the packages installed in that environment
are available as though they were installed normally to your system.

.. code-block:: bash

    source activate phyreengine
    conda install phyreengine

To deactivate the environment, you can run ``source deactivate``.

.. _developing-phyreengine:

Developing PhyreEngine
----------------------

If you plan on developing PhyreEngine, you will probably want to install it as
a *development package*. A development package functions like a normal conda
package, but the source files are kept in their original location, and any
changes to the source files affect the installed package.

First, install ``conda-build`` to the root environment. Then, reactivate
``phyreengine``, move to the PhyreEngine source directory and install the
package:

.. code-block:: bash

    source deactivate
    conda install conda-build
    source activate phyreengine
    cd $PHYRE_ENGINE_SRC
    conda develop .

Installing via :command:`pip` (advanced)
========================================

It is also possible to install PhyreEngine using the :command:`pip` package
manager. PhyreEngine is not available on `PyPI <https://pypi.python.org/pypi>`_
because we strongly recommend using conda, but :command:`pip` may be used from
the source directory of PhyreEngine:

.. code-block:: bash

   pip install .

If you are developing PhyreEngine, you probably want to install it as an
editable package. This will allow you to edit PhyreEngine while python sees it
as a regular installed package.

.. code-block:: bash

    pip install --editable .

You probably won't be able to run :command:`pip install` using the python
version distributed with your operating unless you are an administrator. Even
then, it isn't a great idea to interfere with the system's python installation
because it can cause conflicts with your system's package manager when you come
to update python. In general, the system python installation should only
contain python packages installed using your package manager.

We can work around this problem either by using a *virtualenv*
(virtual environment) or by using pyenv.

Virtualenv
----------

A virtualenv isolates a set of python packages into a "virtual environment,"
and provides tools to use that environment as though it was the default python
installation. It uses your system python as the python interpreter, but allows
you to keep separate installations of different packages for different
projects. Virtualenv is probably already installed on your system or available
in your package manager.

To create and load a new virtual environment for PhyreEngine, run the following:

.. code-block:: bash

    virtualenv phyre_engine
    source phyre_engine/bin/activate

This will install a new virtualenv into the ``phyre_engine`` directory and then
load that virtualenv into your shell. You may then install PhyreEngine using
the commands :ref:`above <getting-started>`. Note that you must the load the
virtualenv using the above ``source`` command each time you use PhyreEngine.


Pyenv
-----

`Pyenv <https://github.com/pyenv/pyenv>`_ is a tool for easily installing and
switching between different python versions. Because pyenv allows for the
installation of entire python interpreters, you are no longer tied to the
version of python distributed with your system. The disadvantage compared to
virtualenv is that you will need a C compiler. Pyenv is the recommended way of
using PhyreEngine.

Detailed instructions are available in the `pyenv README
<https://github.com/pyenv/pyenv/blob/master/README.md>`_. For the impatient,
the following commands will download pyenv, install a new python version and
set that python version as the default:

.. warning::

    These commands will make changes to your ``~/.bash_profile``, an important
    system file. If you are unsure about what each command does, read the pyenv
    README.

.. code-block:: bash

    git clone https://github.com/pyenv/pyenv.git ~/.pyenv
    echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile
    echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile
    echo 'eval "$(pyenv init -)"' >> ~/.bash_profile
    exec $SHELL
    pyenv install 3.6.1
    pyenv global 3.6.1


Running PhyreEngine
===================

Once PhyreEngine is installed the command :command:`phyre_engine` should be
available in your shell. You can check the installation was successful and get
a brief listing of the available command-line options with the following
command:

.. code-block:: bash

    phyre_engine -h

.. code-block:: none

    usage: phyre_engine [-h] [-v] [-e] [-s START] pipeline

    positional arguments:
      pipeline              YAML file describing the pipeline

    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         set verbosity level [default: 1]
      -e, --example         Dump a sample pipeline and exit.
      -s START, --start START
                            Add a value to the initial pipeline state.

The important parameters here are ``pipeline`` and ``--start``. The ``pipeline``
parameter is mandatory, and must point to a `YAML <http://www.yaml.org/>`_ file
describing a pipeline. The ``--start`` option lets you pass data directly to
the pipeline.
