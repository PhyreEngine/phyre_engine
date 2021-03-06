.. _development:

======================
Developing PhyreEngine
======================

Contributions to PhyreEngine are very welcome. This document gives a
description of the code style conventions in use, a recommended workflow for
development, and some information about writing and generating documentation.

Code style
==========

Broadly speaking, PhyreEngine follows the `PEP8
<https://www.python.org/dev/peps/pep-0008/>`_ Python style guides. There is no
need to slavishly stick to the style guide---it is a *guide*, after all---but
try to keep code style consistent. The hard and fast rules are:

1. The only acceptable indentation is four spaces.

2. Modules and functions should be lowercase, with
   ``words_separated_by_underscores``.  Classes should be ``CamelCase``.

Ignoring the other guidelines in PEP8 will cause annoyance, but violating the
rules above will cause confusion and make everyone's life harder.

Development
===========

It is recommended that you use :ref:`conda <installing-via-conda>` to manage
the development of PhyreEngine.  Install the dependencies with ``conda install
phryreengine --only-deps``, and then develop PhyreEngine by running ``conda
develop .`` from within your git clone.

.. TODO: Discuss git hooks, pylint, sphinx

Documentation
=============

PhyreEngine is documented using `Sphinx <http://www.sphinx-doc.org/>`_.  For
our purposes, Sphinx can be thought of as a tool for converting documents
written in `reStructuredText <http://docutils.sourceforge.net/rst.html>` (RST)
into pretty HTML pages, including documents embedded in the docstrings of
PhyreEngine's Python code. Sphinx can do more than that, of course: it can
write several formats, parse Python code, and has a general-purpose extension
framework built in.

Along with long-form documentation such as this document, we have API
documentation generated from docstrings in the project's code. This is grabbed
automatically from the code by Sphinx when it sees an ``.. automodule``, or
``.. autoclass`` (etc) directive, but we still need to tell Sphinx about which
files to document. This is done using the :command:`sphinx-apidoc` command,
which generates RST files with the appropriate ``.. auto*`` directives.

The :file:`Makefile` included in the :file:`doc` subdirectory will take care of
running :command:`sphinx-apidoc` for you. To generate HTML documentation in the
:file:`doc/build/html` directory, simply run :command:`make html` from the
:file:`doc` directory. The :command:`make clean` command may be used to remove
the autogenerated API documentation from the :file:`doc/source` directory.
