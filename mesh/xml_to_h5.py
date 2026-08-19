from fenics import *
import sys

fname = sys.argv[1]

print(fname)

mesh = Mesh(fname + ".xml")                                                  
subdomains = MeshFunction("size_t", mesh, fname + "_physical_region.xml")    
boundaries = MeshFunction("size_t", mesh, fname + "_facet_region.xml")
hdf = HDF5File(mesh.mpi_comm(), fname + ".h5", "w")
hdf.write(mesh, "/mesh")
hdf.write(subdomains, "/subdomains")
hdf.write(boundaries, "/boundaries")
