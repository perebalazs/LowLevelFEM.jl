"""
    line_mesh(; lx=1.0,
                n=10,
                dx=lx/n,
                order=1)

Create a structured 1D mesh of a line domain using Gmsh.

The domain is defined on the interval `[0, lx]`. A transfinite structured
mesh is generated along the line.

# Keyword Arguments
- `lx`: Length of the domain.
- `n`: Default subdivision parameter (used to define `dx`).
- `dx`: Target element size along the line.
- `order`: Polynomial order of the finite elements.

# Physical Groups
The following physical groups are created:
- 0D boundary points: `"left"`, `"right"`
- 1D domain: `"body"`

# Notes
- A transfinite structured mesh is enforced.
- The number of nodes is determined from `ceil(lx/dx) + 1`.
- The mesh is generated but not written to file.
- Initialize Gmsh automatically if needed, clear existing models with
  `gmsh.clear()`, and set Gmsh verbosity to zero. Existing geometry, meshes
  and views are discarded; reconstruct descriptors after remeshing.

Returns `nothing`.

Note: These helper functions are primarily intended for rapid prototyping,
testing and educational examples in LowLevelFEM workflows.
"""
function line_mesh(; lx=1.0, n=10, dx=lx/n, order=1)
    if gmsh.isInitialized() == 0
        gmsh.initialize()
    end
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.clear()
    gmsh.model.add("line_mesh")
    #gmsh.option.setNumber("General.Terminal", 1)

    #gmsh.model.add("line")

    # --- geometry ---
    p1 = gmsh.model.occ.addPoint(0.0, 0.0, 0.0)
    p2 = gmsh.model.occ.addPoint(lx, 0.0, 0.0)

    l = gmsh.model.occ.addLine(p1, p2)

    gmsh.model.occ.synchronize()

    # --- transfinite (structured mesh) ---
    gmsh.model.mesh.setTransfiniteCurve(l, ceil(lx / dx) + 1)

    # --- physical groups ---
    gmsh.model.addPhysicalGroup(0, [p1], 1)  # left
    gmsh.model.setPhysicalName(0, 1, "left")

    gmsh.model.addPhysicalGroup(0, [p2], 2)  # right
    gmsh.model.setPhysicalName(0, 2, "right")

    gmsh.model.addPhysicalGroup(1, [l], 3)   # domain
    gmsh.model.setPhysicalName(1, 3, "body")

    # --- mesh ---
    gmsh.option.setNumber("Mesh.ElementOrder", order)

    gmsh.model.mesh.generate(1)

    return nothing
end

"""
    structured_rect_mesh(; x0=0.0, y0=0.0,
                           lx=1.0, ly=1.0,
                           n=10,
                           dx=lx/n, dy=ly/n,
                           order=1, nx=nothing, ny=nothing)

Create a structured quadrilateral mesh of a rectangular 2D domain using Gmsh.

The domain is defined by its lower-left corner `(x0, y0)` and its dimensions
`lx × ly`. A transfinite structured mesh is generated and recombined into
quadrilateral elements.

# Keyword Arguments
- `x0, y0`: Coordinates of the lower-left corner.
- `lx, ly`: Domain dimensions.
- `n`: Default subdivision parameter (used to define `dx`, `dy`).
- `dx, dy`: Target element sizes in x and y directions.
- `order`: Polynomial order of the finite elements.
- `nx, ny`: Optional positive integer element counts in the x/y directions.
  Each supplied count overrides the size-based subdivision in that direction.
  With neither supplied, the historical common size `min(dx,dy)` is retained.

For a slender rectangle, `nx=25, ny=3` gives 25 longitudinal elements and
3 transverse layers without also refining the long direction. By contrast,
`n=3` retains the common-size rule and can produce many more elements along
an elongated domain.

# Physical Groups
The following physical groups are created:
- 1D boundaries: `"bottom"`, `"right"`, `"top"`, `"left"`
- 2D domain: `"body"`
- 0D corner points: `"leftbottom"`, `"rightbottom"`, `"righttop"`, `"lefttop"`

# Notes
- A transfinite structured mesh is enforced.
- Elements are recombined into quadrilaterals.
- The mesh is generated but not written to file.
- Initialize Gmsh automatically if needed, clear existing models with
  `gmsh.clear()`, and set Gmsh verbosity to zero. Existing geometry, meshes
  and views are discarded; reconstruct descriptors after remeshing.

Returns `nothing`.

Note: These helper functions are primarily intended for rapid prototyping,
testing and educational examples in LowLevelFEM workflows.
# Optional remote points

`leftRemote=false` and `rightRemote=false` optionally create extra point
physical groups `"leftremote"` and `"rightremote"` at the end-section centers.
These provide reference nodes for MPC construction; enabling them alone
creates no rigid coupling or boundary condition.
"""
function structured_rect_mesh(; x0=0.0, y0=0.0, lx=1.0, ly=1.0, n=10, dx=lx / n, dy=ly / n, order=1, nx=nothing, ny=nothing, leftRemote::Bool=false, rightRemote::Bool=false)
    for (name, count) in ((:nx, nx), (:ny, ny))
        count === nothing || (count isa Integer && !(count isa Bool) && count > 0) ||
            throw(ArgumentError("structured_rect_mesh: $name must be a positive integer or nothing."))
    end
    if gmsh.isInitialized() == 0
        gmsh.initialize()
    end
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.clear()
    gmsh.model.add("structured_rect")
    # --------------------------------------------------
    # Geometry
    # --------------------------------------------------
    # Box: origin (0,0,0), size l x l x l
    box = gmsh.model.occ.addRectangle(x0, y0, 0.0, lx, ly)
    if leftRemote
        lr = gmsh.model.occ.addPoint(x0, y0 + ly / 2, 0)
    end
    if rightRemote
        rr = gmsh.model.occ.addPoint(x0 + lx, y0 + ly / 2, 0)
    end

    gmsh.model.occ.synchronize()
    #return

    # --------------------------------------------------
    # Physical groups
    # --------------------------------------------------
    gmsh.model.addPhysicalGroup(1, [1], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [4], -1, "left")

    gmsh.model.addPhysicalGroup(2, [1], -1, "body")

    gmsh.model.addPhysicalGroup(0, [1], -1, "leftbottom")
    gmsh.model.addPhysicalGroup(0, [2], -1, "rightbottom")
    gmsh.model.addPhysicalGroup(0, [3], -1, "righttop")
    gmsh.model.addPhysicalGroup(0, [4], -1, "lefttop")

    if leftRemote
        gmsh.model.addPhysicalGroup(0, [lr], -1, "leftremote")
    end
    if rightRemote
        gmsh.model.addPhysicalGroup(0, [rr], -1, "rightremote")
    end

    # --------------------------------------------------
    # Mesh settings (structured hex mesh)
    # --------------------------------------------------
    #gmsh.option.setNumber("Mesh.Algorithm3D", 1)      # Delaunay for recombination
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    #gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)

    #gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)

    d = min(dx, dy)
    count_x = nx === nothing ? ceil(lx / d) : nx
    count_y = ny === nothing ? ceil(ly / d) : ny

    gmsh.model.mesh.setTransfiniteCurve(1, count_x + 1)
    gmsh.model.mesh.setTransfiniteCurve(3, count_x + 1)
    gmsh.model.mesh.setTransfiniteCurve(2, count_y + 1)
    gmsh.model.mesh.setTransfiniteCurve(4, count_y + 1)

    gmsh.model.mesh.setTransfiniteSurface(1)

    gmsh.option.setNumber("Mesh.ElementOrder", order)

    # Characteristic length control
    #lc = min(dx, dy)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)

    # --------------------------------------------------
    # Generate mesh
    # --------------------------------------------------
    gmsh.model.mesh.generate(2)

    T = [
        1.0, 0.0, 0.0, lx,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    ]

    gmsh.model.mesh.setPeriodic(
        1,
        [2],
        [4],
        T
    )

    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    #gmsh.write("cube.msh")

    #gmsh.finalize()

    return nothing
end

"""
    structured_box_mesh(; x0=0.0, y0=0.0, z0=0.0,
                          lx=1.0, ly=1.0, lz=1.0,
                          n=10,
                          dx=lx/n, dy=ly/n, dz=lz/n,
                          order=1)

Create a structured hexahedral mesh of a rectangular 3D domain using Gmsh.

The domain is defined by its lower-front-left corner `(x0, y0, z0)` and its
dimensions `lx × ly × lz`. A fully transfinite structured mesh is generated
and recombined into hexahedral elements.

# Keyword Arguments
- `x0, y0, z0`: Coordinates of the origin corner.
- `lx, ly, lz`: Domain dimensions.
- `n`: Default subdivision parameter.
- `dx, dy, dz`: Target element sizes in each direction.
- `order`: Polynomial order of the finite elements.

# Physical Groups
The following physical groups are created:
- 2D boundary faces: `"left"`, `"right"`, `"front"`, `"rear"`,
  `"bottom"`, `"top"`
- 3D domain: `"body"`

# Notes
- A transfinite structured mesh is enforced in all directions.
- Elements are recombined into hexahedra.
- The mesh is generated but not written to file.
- Initialize Gmsh automatically if needed, clear existing models with
  `gmsh.clear()`, and set Gmsh verbosity to zero. Existing geometry, meshes
  and views are discarded; reconstruct descriptors after remeshing.

Returns `nothing`.

Note: These helper functions are primarily intended for rapid prototyping,
testing and educational examples in LowLevelFEM workflows.
# Optional remote points

`leftRemote=false` and `rightRemote=false` optionally create extra point
physical groups `"leftremote"` and `"rightremote"` at the end-section centers.
These provide reference nodes for MPC construction; enabling them alone
creates no rigid coupling or boundary condition.
"""
function structured_box_mesh(; x0=0.0, y0=0.0, z0=0.0, lx=1.0, ly=1.0, lz=1.0, n=10, dx=lx / n, dy=ly / n, dz=lz / n, order=1, leftRemote::Bool=false, rightRemote::Bool=false)
    if gmsh.isInitialized() == 0
        gmsh.initialize()
    end
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.clear()
    gmsh.model.add("structured_box")

    # --------------------------------------------------
    # Geometry
    # --------------------------------------------------
    # Box: origin (0,0,0), size l x l x l
    box = gmsh.model.occ.addBox(x0, y0, z0, lx, ly, lz)
    
    if leftRemote
        lr = gmsh.model.occ.addPoint(x0, y0 + ly / 2, z0 + lz / 2)
    end
    if rightRemote
        rr = gmsh.model.occ.addPoint(x0 + lx, y0 + ly / 2, z0 + lz / 2)
    end

    gmsh.model.occ.synchronize()
    #return

    # --------------------------------------------------
    # Physical groups
    # --------------------------------------------------
    gmsh.model.addPhysicalGroup(2, [1], -1, "left")
    gmsh.model.addPhysicalGroup(2, [2], -1, "right")
    gmsh.model.addPhysicalGroup(2, [6], -1, "front")
    gmsh.model.addPhysicalGroup(2, [5], -1, "rear")
    gmsh.model.addPhysicalGroup(2, [3], -1, "bottom")
    gmsh.model.addPhysicalGroup(2, [4], -1, "top")

    gmsh.model.addPhysicalGroup(3, [1], -1, "body")

    if leftRemote
        gmsh.model.addPhysicalGroup(0, [lr], -1, "leftremote")
    end
    if rightRemote
        gmsh.model.addPhysicalGroup(0, [rr], -1, "rightremote")
    end

    # --------------------------------------------------
    # Mesh settings (structured hex mesh)
    # --------------------------------------------------
    #gmsh.option.setNumber("Mesh.Algorithm3D", 1)      # Delaunay for recombination
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    #gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)

    #gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)

    d = min(dx, dy, dz)

    gmsh.model.mesh.setTransfiniteCurve(9, ceil(lx / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(10, ceil(lx / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(11, ceil(lx / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(12, ceil(lx / d) + 1)

    gmsh.model.mesh.setTransfiniteCurve(2, ceil(ly / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(4, ceil(ly / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(6, ceil(ly / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(8, ceil(ly / d) + 1)

    gmsh.model.mesh.setTransfiniteCurve(1, ceil(lz / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(3, ceil(lz / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(5, ceil(lz / d) + 1)
    gmsh.model.mesh.setTransfiniteCurve(7, ceil(lz / d) + 1)

    gmsh.model.mesh.setTransfiniteSurface(1)
    gmsh.model.mesh.setTransfiniteSurface(2)
    gmsh.model.mesh.setTransfiniteSurface(3)
    gmsh.model.mesh.setTransfiniteSurface(4)
    gmsh.model.mesh.setTransfiniteSurface(5)
    gmsh.model.mesh.setTransfiniteSurface(6)

    gmsh.model.mesh.setTransfiniteVolume(1)

    #gmsh.model.mesh.setTransfiniteAutomatic()

    gmsh.option.setNumber("Mesh.ElementOrder", order)

    # Characteristic length control
    #lc = min(dx, dy, dz)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    #gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)

    # --------------------------------------------------
    # Generate mesh
    # --------------------------------------------------
    gmsh.model.mesh.generate(3)
    
    T = [
        1.0, 0.0, 0.0, lx,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    ]

    gmsh.model.mesh.setPeriodic(
        2,
        [2],
        [1],
        T
    )

    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    #gmsh.write("cube.msh")

    #gmsh.finalize()

    return nothing
end

"""
    openGeometry(name::String; verbose::Int=0)

Open a geometry or mesh file in the current Gmsh session.

If Gmsh is not yet initialized, the function calls `gmsh.initialize()`.
If the file `name` does not exist and has extension `.geo`, a new geometry
file is created and initialized with a minimal header. Otherwise an error
is thrown. `verbose` sets the verbosity level of Gmsh.

Before opening the file the current Gmsh model is cleared with `gmsh.clear()`,
so any previously loaded geometry or mesh is removed.

# Arguments
- `name::String`: Path to a geometry (`.geo`) or mesh file (`.msh`, `.step`,
  `.brep`, etc.) that should be opened in Gmsh.

# Behavior
- Existing files are opened directly.
- Missing `.geo` files are automatically created.
- Missing files with other extensions result in an error.

# Notes
This function only loads the geometry into Gmsh. The model name used later
by `Problem` is obtained from `gmsh.model.getCurrent()`.

# Examples
```julia
openGeometry("beam.geo")     # create if missing and open
openGeometry("mesh.msh")     # open existing mesh
```
"""
function openGeometry(name::String; verbose::Int=0)

    if gmsh.isInitialized() == 0
        gmsh.initialize()
    end
    
    gmsh.option.setNumber("General.Verbosity", verbose)

    if !isfile(name)
        ext = lowercase(splitext(name)[2])

        if ext == ".geo"
            open(name, "w") do io
                println(io, "//+ LowLevelFEM")
            end
        else
            error("File $name does not exist")
        end
    end

    gmsh.clear()
    gmsh.open(name)

end

"""
    functionName(arg1, arg2; kw1=..., kw2=...)

Short one-sentence description of what the function does.

# Arguments
- `arg1::Type`: description
- `arg2::Type`: description

# Keyword arguments
- `kw1`: description
- `kw2`: description

# Returns
Description of the returned object.

# Notes
Optional implementation notes or behavior details.

# Examples
```julia
example code
```
"""
