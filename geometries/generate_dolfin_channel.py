__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2012-01-17"
__copyright__ = "Copyright (C) 2012 Anders Logg"
__license__  = "GNU LGPL version 3 or any later version"

from dolfin import *
from dolfin.mesh.netgen import *

mesh = Mesh("dolfin-1.xml")

file = File("dolfin_channel.xml")

x = mesh.coordinates()
for i in range(len(x)):
    xx = 1.0 - x[i][0]
    if abs(xx - 0.0) < DOLFIN_EPS:
        xx = 0.0
    if abs(xx - 1.0) < DOLFIN_EPS:
        xx = 1.0
    x[i][0] = xx

file << mesh
