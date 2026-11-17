Our Docs
========

Sphinx builds this documentation and we integrate C++ documentation built by Doxygen using the
Breathe extension.

You can build the documentation locally by building the Doxygen .xml files
followed by the Sphinx .html files, e.g.

.. code-block:: console

  $ cd ./docs && mkdir -p build/doxygen
  $ doxygen doxygen/doxygen.dox && uv run make html

which you can then view in your preferred browser e.g.

.. code-block:: console

  $ open build/html/index.html


More Resources
---------------

Sphinx: https://www.sphinx-doc.org/en/master/

Breathe: https://breathe.readthedocs.io/en/latest/

Doxygen: https://www.doxygen.nl


Need More Documentation?
---------------------------
If you would like to contribute to the documentation, pl, have questions or would like clarification or
more detail on a particular topic, please :ref:`contact us! <contact>`.
