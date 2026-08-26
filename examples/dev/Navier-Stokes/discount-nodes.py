import gmsh

def check(filename):
    gmsh.initialize()
    gmsh.open(filename)

    nodeTags, _, _ = gmsh.model.mesh.getNodes()

    n = len(nodeTags)

    print(filename)
    print("Number of nodes :", n)
    print("Minimum node tag:", min(nodeTags))
    print("Maximum node tag:", max(nodeTags))
    print("Contiguous      :", min(nodeTags) == 1 and max(nodeTags) == n)
    print()

    gmsh.finalize()


check("pipe-cont.geo")
check("pipe-discont.geo")