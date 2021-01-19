__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2012-01-17"
__copyright__ = "Copyright (C) 2012 Anders Logg"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from dolfin.mesh.netgen import *

geometry = LEGO(10, 2, 3)
mesh = Mesh(geometry)

file = File("lego_beam.xml")
file << mesh
