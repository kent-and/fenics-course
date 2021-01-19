"""
Converts a FEniCS Python script into an Jupyter notebook file

Usage:

    python ipython_converter.py input.py output.ipynb

"""

import sys
import nbformat

header = """%matplotlib inline
from fenics import *
parameters["plotting_backend"] = 'matplotlib'"""

main = open(sys.argv[1], "r").readlines()
has_pylab = False

# Do some necessary replacments
for i in range(len(main)):
    main[i] = main[i].replace(", interactive=True", "")
    main[i] = main[i].replace("from fenics import *", "")
    main[i] = main[i].replace("interactive()", "")
    main[i] = main[i].replace(", rescale=True", "")
    if "plot(" in main[i]:
        if not has_pylab:
            header += "\nimport pylab"
            has_pylab = True
        leading_spaces = len(main[i]) - len(main[i].lstrip())
        main.insert(i+1, " "*leading_spaces + "pylab.show()\n")


sources = [header, main]
nb = nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(source=src) for src in sources])

nbformat.write(nb, sys.argv[2])
