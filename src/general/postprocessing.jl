"""
    showDoFResults(q::Union{ScalarField,VectorField,TensorField}, comp::Symbol; 
                   name=comp, visible=false, factor=0)

Loads nodal results into a View in Gmsh. `q` is the field to show, `comp` is
the component of the field (:vector, :uvec, :ux, :uy, :uz, :vvec, :vx, :vy, :vz,
:qvec, :qx, :qy, :qz, :T, :p, :qn, :s, :sx, :sy, :sz, :sxy, :syx, :syz,
:szy, :szx, :sxz, :e, :ex, :ey, :ez, :exy, :eyx, :eyz, :ezy, :ezx, :exz, :x, :y, :z, :xy, :yz, :zx, :yx, :zy, :xz, :seqv, :scalar, :tensor),
`name` is a title to display and `visible` is a Boolean to toggle the initial visibility in Gmsh on or off.
If `q` has more columns, then a sequence of results will be shown (e.g., as an animation).
`factor` multiplies the DoF result to increase for better visibility.
This function returns the tag of the View.

Returns: `tag`

Types:
- `q`: ScalarField, VectorField or TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showDoFResults(q::Union{ScalarField,VectorField,TensorField}, comp::Symbol; name=comp, visible=false, ff = 0, factor=0)
    problem = q.model
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    t = q.t
    #display(q)
    if q.a == [;;]
        q = elementsToNodes(q)
    end
    if q.a == [;;]
        error("showDoFResults: No data")
    end
    edim = problem.dim
    if problem.type == :Truss
        edim = 1
    end
    pdim = problem.pdim
    pdim = div(size(q.a,1), problem.non)
    #display("showDoFResults: dim=$pdim, type=$(q.type).")
    nodeTags = []
    u = []
    ##############################################################################
    #if problem.type == :Reynolds || problem.type == :NavierStokes
    #    phName = problem.material.phName
    #    tag = getTagForPhysicalName(phName)
    #    nT, coords = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
    #    append!(nodeTags, nT)
    ##############################################################################
    #else #########################################################################
    for ipg in 1:length(problem.material)
        phName = problem.material[ipg].phName
        tag = getTagForPhysicalName(phName)
        edim = getDimForPhysicalName(phName)
        nT, coords = gmsh.model.mesh.getNodesForPhysicalGroup(edim, tag)
        append!(nodeTags, nT)
    end
    if problem.geometry.nameVolume ≠ "" && problem.geometry.nameGap ≠ ""
        phName = problem.geometry.nameVolume
        tag = getTagForPhysicalName(phName)
        edim = getDimForPhysicalName(phName)
        nT, coords = gmsh.model.mesh.getNodesForPhysicalGroup(edim, tag)
        append!(nodeTags, nT)
    end
    #end #########################################################################
    
    #nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim, -1, true)
    non = length(nodeTags)
    uvec = gmsh.view.add(name)
    if size(q.a, 2) != length(t)
        error("showDoFResults: number of time steps missmatch ($(size(q.a,2)) <==> $(length(t))).")
    end
    if comp == :s || comp == :e || q isa TensorField
        pdim = 9
    end
    for j in 1:length(t)
        k = 1im
        if comp == :uvec || comp == :vvec || comp == :qvec || comp == :vector
            nc = 3
            u = zeros(3 * non)
            for i in 1:length(nodeTags)
                u[3i-2] = q.a[pdim*nodeTags[i]-(pdim-1), j]
                u[3i-1] = pdim > 1 ? q.a[pdim*nodeTags[i]-(pdim-2), j] : 0
                u[3i-0] = pdim == 3 ? q.a[pdim*nodeTags[i]-(pdim-3), j] : 0
            end
        elseif comp == :s || comp == :e || comp == :tensor
            nc = 9
            u = zeros(9 * non)
            for i in 1:length(nodeTags)
                u[9i-8] = q.a[pdim*nodeTags[i]-(pdim-1), j]
                u[9i-7] = q.a[pdim*nodeTags[i]-(pdim-2), j]
                u[9i-6] = q.a[pdim*nodeTags[i]-(pdim-3), j]
                u[9i-5] = q.a[pdim*nodeTags[i]-(pdim-4), j]
                u[9i-4] = q.a[pdim*nodeTags[i]-(pdim-5), j]
                u[9i-3] = q.a[pdim*nodeTags[i]-(pdim-6), j]
                u[9i-2] = q.a[pdim*nodeTags[i]-(pdim-7), j]
                u[9i-1] = q.a[pdim*nodeTags[i]-(pdim-8), j]
                u[9i-0] = q.a[pdim*nodeTags[i]-(pdim-9), j]
            end
        else
            nc = 1
            if comp == :ux || comp == :vx || comp == :p || comp == :T || comp == :qx  || comp == :qn || comp == :sx || comp == :ex || comp == :scalar || (comp == :x && q isa VectorField) || ((comp == :x || comp == :xx) && q isa TensorField)
                k = 1
            elseif comp == :uy || comp == :vy || comp == :qy || comp == :syx || comp == :eyx || (comp == :y && q isa VectorField) || (comp == :yx && q isa TensorField)
                k = 2
            elseif comp == :uz || comp == :vz || comp == :qz || comp == :szx || comp == :ezx || (comp == :z && q isa VectorField) || (comp == :zx && q isa TensorField)
                k = 3
            elseif comp == :sxy || comp == :exy || (comp == :xy && q isa TensorField)
                k = 4
            elseif comp == :sy || comp == :ey || ((comp == :y || comp == :yy) && q isa TensorField)
                k = 5
            elseif comp == :szy || comp == :ezy || (comp == :zy && q isa TensorField)
                k = 6
            elseif comp == :sxz || comp == :exz || (comp == :xz && q isa TensorField)
                k = 7
            elseif comp == :syz || comp == :eyz || (comp == :yz && q isa TensorField)
                k = 8
            elseif comp == :sz || comp == :ez || ((comp == :z || comp == :zz) && q isa TensorField)
                k = 9
            elseif comp == :seqv
                k = 10
            else
                error("ShowDisplacementResults: component is $comp ????")
            end
            u = zeros(non)
            for i in 1:length(nodeTags)
                null = k <= 9 ? 0 : √(0.5 * ((q.a[pdim*nodeTags[i]-(pdim-1), j]-q[pdim*nodeTags[i]-(pdim-5), j])^2+(q[pdim*nodeTags[i]-(pdim-5), j]-q[pdim*nodeTags[i]-(pdim-9), j])^2+(q[pdim*nodeTags[i]-(pdim-9), j]-q[pdim*nodeTags[i]-(pdim-1), j])^2 + 6*(q[pdim*nodeTags[i]-(pdim-2), j]^2+q[pdim*nodeTags[i]-(pdim-3), j]^2+q[pdim*nodeTags[i]-(pdim-6), j]^2)))
                u[i] = k > pdim ? null : q.a[pdim*nodeTags[i]-(pdim-k), j]
            end
        end
        gmsh.view.addHomogeneousModelData(uvec, j-1, problem.name, "NodeData", nodeTags, u, t[j], nc)
    end
    
    gmsh.view.option.setNumber(uvec, "GlyphLocation", 2)
    #gmsh.view.option.setNumber(uvec, "DisplacementFactor", 0)
    eigenmode = q.kind in (:eigen_real, :eigen_imag)
    scale = factor != 0 ? factor : (eigenmode ? 1.0 : 0.0)
    gmsh.view.option.setNumber(uvec, "DisplacementFactor", scale)
    if ff == 1 || ff == 2
        gmsh.view.option.setNumber(uvec, "AdaptVisualizationGrid", 1)
    else
        gmsh.view.option.setNumber(uvec, "AdaptVisualizationGrid", 0)
    end
    gmsh.view.option.setNumber(uvec, "TargetError", -1e-4)
    gmsh.view.option.setNumber(uvec, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(uvec, "Visible", 0)
    elseif visible == true
        gmsh.view.option.setNumber(uvec, "Visible", 1)
    end
    gmsh.view.option.setNumber(uvec, "ShowTime", 0)
    if q.kind == :eigen_real
    gmsh.view.option.setNumber(uvec, "ShowTime", 6)

    elseif q.kind == :eigen_imag
        gmsh.view.option.setNumber(uvec, "ShowTime", 7)

    elseif ff == 1 || ff == 2
        gmsh.view.option.setNumber(uvec, "ShowTime", 6)

    elseif ff == 0 && length(t) > 1
        gmsh.view.option.setNumber(uvec, "ShowTime", 1)

    else
        gmsh.view.option.setNumber(uvec, "ShowTime", 0)
    end
    if q isa VectorField &&
        (q.kind in (:eigen_real, :eigen_imag) || ff in (1, 2))
        gmsh.view.option.setNumber(uvec, "VectorType", 5)
    end
    #=if ff == 0 && length(t) > 1
        gmsh.view.option.setNumber(uvec, "ShowTime", 1)
    elseif ff == 1 || ff == 2 || q.kind == :eigen_real
        gmsh.view.option.setNumber(uvec, "ShowTime", 6)
        gmsh.view.option.setNumber(uvec, "VectorType", 5)
    end=#
    if factor != 0
        gmsh.view.option.setNumber(uvec, "VectorType", 5)
        gmsh.view.option.setNumber(uvec, "DisplacementFactor", factor)
    end
    return uvec
end

function showDoFResults(q::Union{ScalarField,VectorField,TensorField}; name=q.type, visible=false, ff = 0, factor=0)
    if q isa ScalarField
        showDoFResults(q, :scalar, name=name, visible=visible, ff=ff, factor=factor)
    elseif q isa VectorField
        showDoFResults(q, :vector, name=name, visible=visible, ff=ff, factor=factor)
    elseif q isa TensorField
        showDoFResults(q, :tensor, name=name, visible=visible, ff=ff, factor=factor)
    else
        error("showDoFResults: argument must be a ScalarField, VectorField or TensorField.")
    end
end

"""
    showModalResults(Φ::Eigen; name="", visible=false)

Loads modal results into a View in Gmsh. `Φ` is an `Eigen` struct. `name` is a title to display and
`visible` is a Boolean to toggle the initial visibility in Gmsh on or off. Click on ▷| to change
the results. This function returns the tag of the View.

Returns: `tag`

Types:
- `Φ`: Eigen
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showModalResults(Φ::Eigen; name=:modal, visible=false, ff=1)
    return showDoFResults(VectorField([], Φ.ϕ, Φ.f, [], length(Φ.f), :v3D, Φ.model), :uvec, name=name, visible=visible, ff=ff)
end

"""
    showBucklingResults(Φ::Eigen; name="", visible=false)

Loads buckling results into a View in Gmsh. `Φ` is an `Eigen` struct. `name` is a title to display and
`visible` is a Boolean to toggle the initial visibility in Gmsh on or off. Click on ▷| to change
the results. This function returns the tag of the View.

Returns: `tag`

Types:
- `Φ`: Eigen
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showBucklingResults(Φ::Eigen; name="buckling", visible=false, ff=2)
    return showDoFResults(VectorField([], Φ.ϕ, Φ.f, [], length(Φ.f), :v3D, Φ.model), :uvec, name=name, visible=visible, ff=ff)
end

"""
    showStrainResults(E::TensorField, comp::Symbol; 
                      name=comp, visible=false, smooth=true)

Loads strain results into a View in Gmsh. `E` is a strain field to show, `comp` is
the component of the field (:e, :ex, :ey, :ez, :exy, :eyz, :ezx, :x, :y, :z, :xy, :yz, :zx, :yx, :zy, :xz),
`name` is a title to display, `visible` is a Boolean to toggle the initial visibility in Gmsh on or
off and `smooth` is a Boolean to toggle
smoothing the stress field on or off. If `E` contains more than one time steps, then a 
sequence of results will be shown (e.g., as an animation). This function returns
the tag of the View.

Returns: `tag`

Types:
- `E`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
# Display configuration

The keyword `factor` is accepted by this view constructor; the current call
signature is:

    showStrainResults(E0, comp; name=comp, visible=false, smooth=true)

View options affect visualization only, not the stored finite-element field.
"""
function showStrainResults(E0, comp; name=comp, visible=false, smooth=true)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    E = nodesToElements(E0)
    problem = E.model
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    t = E.t
    if E.nsteps != length(t)
        error("showStrainResults: number of time steps missmatch ($(E.nsteps) <==> $(length(t))).")
    end
    EE = gmsh.view.add(name)
    if E.A == []
        error("showStrainResults: No data")
    end
    ε = E.A
    numElem = E.numElem
    for jj in 1:length(t)
        
        k = 1im
        if comp == :e || comp == :tensor
            εcomp = [ε[i][:,jj] for i in 1:length(E.numElem)]
            nc = 9
        else
            nc = 1
            if comp == :ex || ((comp == :x || comp == :xx) && E isa TensorField)
                k = 8
            elseif comp == :ey || ((comp == :y || comp == :yy) && E isa TensorField)
                k = 4
            elseif comp == :ez || ((comp == :z || comp == :zz) && E isa TensorField)
                k = 0
            elseif comp == :exy || comp == :xy
                k = 5
            elseif comp == :eyz || comp == :yz
                k = 1
            elseif comp == :ezx || comp == :zx
                k = 6
            elseif comp == :eyx || comp == :yx
                k = 7
            elseif comp == :ezy || comp == :zy
                k = 3
            elseif comp == :exz || comp == :xz
                k = 2
            else
                error("ShowStrainResults: component is $comp ????")
            end
            εcomp = []
            sizehint!(εcomp, length(numElem))
            for i in 1:length(E.numElem)
                ex = zeros(div(size(ε[i], 1), 9))
                for j in 1:(div(size(ε[i], 1), 9))
                    ex[j] = ε[i][9j-k, jj]
                end
                push!(εcomp, ex)
            end
        end
        gmsh.view.addModelData(EE, jj-1, problem.name, "ElementNodeData", numElem, εcomp, t[jj], nc)
    end
    
    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end
    
    gmsh.view.option.setNumber(EE, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(EE, "TargetError", -1e-4)
    gmsh.view.option.setNumber(EE, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(EE, "Visible", 0)
    end
    gmsh.view.option.setNumber(EE, "ShowTime", 0)
    if length(t) > 1
        gmsh.view.option.setNumber(EE, "ShowTime", 1)
    end
    #display("$comp..ok")
    return EE
end

"""
    showElementResults(F::Union{ScalarField,VectorField,TensorField}, comp=:Symbol; 
                       name=comp, visible=false, smooth=false, factor=0)

    showElementResults(F, comp; name=..., visible=..., smooth=...)

Same as `ShowStressResults` or `showStrainResults`, depending on the type of `F` data field.

Return: `tag`

Types:
- `F`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
# Display configuration

The keyword `ff` is accepted by this view constructor; the current call
signature is:

    showElementResults(F::Union{ScalarField,VectorField,TensorField}, comp; name=comp, visible=false, smooth=false, factor=0)

View options affect visualization only, not the stored finite-element field.
# Display configuration

The keyword `factor` is accepted by this view constructor; the current call
signature is:

    showStressResults(q::TensorField; name="StressField", visible=false, smooth=false, ff = 0, factor=0)

View options affect visualization only, not the stored finite-element field.
"""
function showElementResults(F::Union{ScalarField,VectorField,TensorField}, comp; name=comp, visible=false, smooth=false, factor=0)
    if F.type == :e || F.type == :tensor && isElementwise(F)
        return showStrainResults(F, comp, name=name, visible=visible, smooth=smooth)
    elseif F.type == :s && isElementwise(F)
        return showStressResults(F, comp, name=name, visible=visible, smooth=smooth)
    elseif F isa VectorField && isElementwise(F)
        return showHeatFluxResults(F, comp, name=name, visible=visible, smooth=smooth, factor=factor)
    elseif F isa ScalarField && isElementwise(F)
        return showScalarResults(F, name=name, visible=visible, smooth=smooth, factor=factor)
    elseif isNodal(F)
        return showDoFResults(F, comp, name=name, visible=visible, factor=factor)
    else
        error("showElementResults: type is '$(F.type)'")
    end
end

function showElementResults(q::Union{ScalarField,VectorField,TensorField}; name="__field__", visible=false, smooth=false, ff = 0, factor=0)
    if q isa ScalarField
        showElementResults(q, :scalar, name=name == "__field__" ? :scalar : name, visible=visible, smooth=smooth)
    elseif q isa VectorField
        showElementResults(q, :vector, name=name == "__field__" ? :vector : name, visible=visible, smooth=smooth, factor=factor)
    elseif q isa TensorField
        showElementResults(q, :tensor, name=name == "__field__" ? :tensor : name, visible=visible, smooth=smooth)
    else
        error("showElementResults: argument must be a ScalarField, VectorField or TensorField.")
    end
end

function showStressResults(q::TensorField; name="StressField", visible=false, smooth=false, ff = 0, factor=0)
    if q isa TensorField
        showStressResults(q, :s, name=name, visible=visible, smooth=smooth)
    else
        error("showStressResults: argument must be a TensorField.")
    end
end

function showStrainResults(q::TensorField; name="StrainField", visible=false, smooth=false, ff = 0, factor=0)
    if q isa TensorField
        showStrainResults(q, :e, name=name, visible=visible, smooth=smooth)
    else
        error("showStrainResults: argument must be a TensorField.")
    end
end

function showHeatFluxResults(q::VectorField; name="VectorField", visible=false, smooth=false, ff = 0, factor=0)
    if q isa VectorField
        showHeatFluxResults(q, :vector, name=name, visible=visible, smooth=smooth)
    else
        error("showHeatFluxResults: argument must be a VectorField.")
    end
end

"""
    showStressResults(S0::TensorField, comp::Symbol; 
                      name=comp, visible=false, smooth=false)

Loads stress results into a View in Gmsh. `S` is a stress field to show, `comp` is
the component of the field (:s, :sx, :sy, :sz, :sxy, :syz, :szx, :seqv, :x, :y, :z, :xy, :yz, :zx, :yx, :zy, :xz),
`name` is a title to display, `visible` is a Boolean to toggle the initial visibility in Gmsh on or
off, and `smooth` is a Boolean to toggle smoothing the stress field on or off. If `S` contains more
than one time step, a sequence of results will be shown (e.g., as an animation). This function returns
the tag of the View.

Returns: `tag`

Types:
- `S`: TensorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showStressResults(S0::TensorField, comp; name=comp, visible=false, smooth=false)
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    S = nodesToElements(S0)
    problem = S.model
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    t = S.t
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    if S.nsteps != length(t)
        error("showStressResults: number of time steps missmatch ($(S.nsteps) <==> $(length(t))).")
    end
    SS = gmsh.view.add(name)
    if S.A == []
        error("showStressResults: No data")
    end
    σ = S.A
    numElem = S.numElem
    for jj in 1:length(t)
        
        k = 1im
        if comp == :s || comp == :tensor
            σcomp = [σ[i][:,jj] for i in 1:length(S.numElem)]
            nc = 9
        elseif comp == :seqv
            nc = 1
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                seqv = zeros(div(size(σ[i], 1), 9))
                for j in 1:(div(size(σ[i], 1), 9))
                    sx = σ[i][9j-8, jj]
                    sy = σ[i][9j-4, jj]
                    sz = σ[i][9j-0, jj]
                    sxy = σ[i][9j-7, jj]
                    syz = σ[i][9j-3, jj]
                    szx = σ[i][9j-6, jj]
                    seqv[j] = √(((sx-sy)^2+(sy-sz)^2+(sz-sx)^2)/2. + 3*(sxy^2+syz^2+szx^2))
                end
                push!(σcomp, seqv)
            end
        else
            nc = 1
            if comp == :sx || comp == :x
                k = 8
            elseif comp == :sy || comp == :y
                k = 4
            elseif comp == :sz || comp == :z
                k = 0
            elseif comp == :sxy || comp == :syx || comp == :xy || comp == :yx
                k = 7                                                       
            elseif comp == :syz || comp == :szy || comp == :yz || comp == :zy
                k = 3                                                       
            elseif comp == :szx || comp == :sxz || comp == :zx || comp == :xz
                k = 6
            else
                error("ShowStressResults: component is $comp ????")
            end
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                sx = zeros(div(size(σ[i], 1), 9))
                for j in 1:(div(size(σ[i], 1), 9))
                    sx[j] = σ[i][9j-k, jj]
                end
                push!(σcomp, sx)
            end
        end
        gmsh.view.addModelData(SS, jj-1, problem.name, "ElementNodeData", numElem, σcomp, t[jj], nc)
    end
    
    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end
    
    gmsh.view.option.setNumber(SS, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(SS, "TargetError", -1e-4)
    gmsh.view.option.setNumber(SS, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(SS, "Visible", 0)
    end
    gmsh.view.option.setNumber(SS, "ShowTime", 0)
    if length(t) > 1
        gmsh.view.option.setNumber(SS, "ShowTime", 1)
    end
    #display("$comp..ok")
    return SS
end

"""
    showHeatFluxResults(S0::VectorField, comp::Symbol; 
                        name=comp, visible=false, smooth=true, factor=0)

Loads heat flux results into a View in Gmsh. `Q` is a heat flux field to show, `comp` is
the component of the field (:qvec, :qx, :qy, :qz, :q, :x, :y, :z),
`name` is a title to display, `visible` is a Boolean to toggle the initial visibility in Gmsh on or
off, and `smooth` is a Boolean to toggle smoothing on or off. If `Q` contains more than one time step,
a sequence of results will be shown (e.g., as an animation). This function returns the tag of the View.

Returns: `tag`

Types:
- `S`: VectorField
- `comp`: Symbol
- `name`: String
- `visible`: Boolean
- `smooth`: Boolean
- `tag`: Integer
"""
function showHeatFluxResults(S0::VectorField, comp; name=comp, visible=false, smooth=true, factor=0)
    S = nodesToElements(S0)
    problem = S.model
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    if S.type == :v3D
        dim = 3
    elseif S.type == :v2D
        dim = 2
    else
        error("showHeatFluxResults: type is $(S.type)")
    end
    t = S.t
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    if S.nsteps != length(t)
        error("showStressResults: number of time steps missmatch ($(S.nsteps) <==> $(length(t))).")
    end
    SS = gmsh.view.add(name)
    if S.A == []
        error("showHeatFluxResults: No data")
    end
    σ = S.A
    numElem = S.numElem
    for jj in 1:length(t)
        
        k = 1im
        if comp == :qvec || comp == :vector
            σcomp = []
            for i in eachindex(S.numElem)
                es = length(σ[i][:,jj]) ÷ dim
                aa = zeros(3es)
                aa[1:3:3es] .= σ[i][1:dim:dim*es, jj]
                aa[2:3:3es] .= σ[i][2:dim:dim*es, jj]
                aa[3:3:3es] .= dim == 3 ? σ[i][3:dim:dim*es, jj] : zeros(es,1)
                push!(σcomp, aa)
            end
            #σcomp = [σ[i][:,jj] for i in 1:length(S.numElem)]
            nc = 3
        elseif comp == :q
            nc = 1
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                seqv = zeros(div(size(σ[i], 1), dim))
                for j in 1:(div(size(σ[i], 1), dim))
                    sx = σ[i][dim*j-(dim-1), jj]
                    sy = σ[i][dim*j-(dim-2), jj]
                    sz = dim == 3 ? σ[i][dim*j-(dim-3), jj] : 0
                    seqv[j] = √(sx^2+sy^2+sz^2)
                end
                push!(σcomp, seqv)
            end
        else
            nc = 1
            if comp == :qx || comp == :x
                k = 1
            elseif comp == :qy || comp == :y
                k = 2
            elseif comp == :qz || comp == :z && dim == 3
                k = 3
            else
                error("ShowHeatFluxResults: component is $comp ????")
            end
            σcomp = []
            sizehint!(σcomp, length(numElem))
            for i in 1:length(S.numElem)
                ss = zeros(div(size(σ[i], 1), dim))
                for j in 1:(div(size(σ[i], 1), dim))
                    ss[j] = σ[i][dim*j-(dim-k), jj]
                end
                push!(σcomp, ss)
            end
        end
        gmsh.view.addModelData(SS, jj-1, problem.name, "ElementNodeData", numElem, σcomp, t[jj], nc)
    end
    
    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end
    
    gmsh.view.option.setNumber(SS, "GlyphLocation", 2)
    gmsh.view.option.setNumber(SS, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(SS, "TargetError", -1e-4)
    gmsh.view.option.setNumber(SS, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(SS, "Visible", 0)
    end
    gmsh.view.option.setNumber(SS, "ShowTime", 0)
    if length(t) > 1
        gmsh.view.option.setNumber(SS, "ShowTime", 1)
    end
    gmsh.view.option.setNumber(SS, "DisplacementFactor", factor)
    #display("$comp..ok")
    return SS
end

function showScalarResults(S0; name="ScalarField", visible=false, smooth=false, factor=0)
    S = nodesToElements(S0)
    problem = S.model
    #gmsh.fltk.openTreeItem("0Modules/Post-processing")
    gmsh.model.setCurrent(problem.name)
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)
    dim = problem.dim
    dim = 1
    t = S.t
    #elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, -1)
    #elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(elemTypes[1])
    if S.nsteps != length(t)
        error("showStressResults: number of time steps missmatch ($(S.nsteps) <==> $(length(t))).")
    end
    SS = gmsh.view.add(name)
    if S.A == []
        error("showHeatFluxResults: No data")
    end
    σ = S.A
    numElem = S.numElem
    for jj in 1:length(t)
        
        nc = 1
        σcomp = []
        sizehint!(σcomp, length(numElem))
        for i in 1:length(S.numElem)
            #=
            sc = zeros(div(size(σ[i], 1), dim))
            for j in 1:(div(size(σ[i], 1), dim))
                sc[j] = σ[i][j, jj]
            end
            push!(σcomp, sc)
            =#
            push!(σcomp, σ[i][:, jj])
        end
        gmsh.view.addModelData(SS, jj-1, problem.name, "ElementNodeData", numElem, σcomp, t[jj], nc)
    end
    
    if smooth == true
        gmsh.plugin.setNumber("Smooth", "View", -1)
        gmsh.plugin.run("Smooth")
    end
    
    gmsh.view.option.setNumber(SS, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(SS, "TargetError", -1e-4)
    gmsh.view.option.setNumber(SS, "MaxRecursionLevel", 1)
    if visible == false
        gmsh.view.option.setNumber(SS, "Visible", 0)
    end
    gmsh.view.option.setNumber(SS, "ShowTime", 0)
    if length(t) > 1
        gmsh.view.option.setNumber(SS, "ShowTime", 1)
    end
    gmsh.view.option.setNumber(SS, "DisplacementFactor", factor)
    #display("$comp..ok")
    return SS
end

"""
    plotOnPath(pathName, field; points=100, step=nothing, plot=false,
               name="field [\$field] on \$pathName", visible=false,
               offsetX=0, offsetY=0, offsetZ=0)

Sample a Gmsh model-data view along the curves of a physical group.
`step` is a zero-based view step; `nothing` selects all steps. All selected
steps are retrieved in one probe per sampling point. Vector values are shown
as magnitudes and tensor values as von Mises equivalents.
Return the new view tag, or `(tag, (x, y))` when `plot=true`; `y` contains the
first selected step. Samples are uniform in curve parameter, as in the
original implementation.
"""
function plotOnPath(pathName, field; points=100, step=nothing, plot=false,
    name="field [$field] on $pathName", visible=false,
    offsetX=0, offsetY=0, offsetZ=0)
    points isa Integer && points >= 2 ||
        throw(ArgumentError("plotOnPath: points must be an integer >= 2."))
    dimTags = gmsh.model.getEntitiesForPhysicalName(pathName)
    isempty(dimTags) && error("plotOnPath: physical group '$pathName' was not found.")
    all(dt -> dt[1] == 1, dimTags) ||
        error("plotOnPath: physical group '$pathName' must contain only curves.")
    nsteps = Int(gmsh.view.option.getNumber(field, "NbTimeStep"))
    nsteps > 0 || error("plotOnPath: the source view has no time steps.")
    if step !== nothing
        step isa Integer && step >= 0 ||
            throw(ArgumentError("plotOnPath: step must be a nonnegative integer or nothing."))
        if step >= nsteps
            @warn "plotOnPath: step exceeds the available range; using the last step." step nsteps
        end
    end
    probe_step = step === nothing ? -1 : min(step, nsteps - 1)
    count = step === nothing ? nsteps : 1
    _, _, _, _, nc = gmsh.view.getHomogeneousModelData(field, 0)
    nc in (1, 3, 9) || error("plotOnPath: unsupported component count $nc.")
    nsamples = points * length(dimTags)
    stride = 3 + count
    data = Vector{Float64}(undef, nsamples * stride)
    x = plot ? Vector{Float64}(undef, nsamples) : Float64[]
    y = plot ? Vector{Float64}(undef, nsamples) : Float64[]
    origin = zeros(3)
    previous = zeros(3)
    sample = 0
    for (_, path) in dimTags
        bounds = gmsh.model.getParametrizationBounds(1, path)
        lo, hi = bounds[1][1], bounds[2][1]
        for i in 1:points
            xyz = gmsh.model.getValue(1, path, [lo + (i - 1) * (hi - lo) / (points - 1)])
            sample += 1
            sample == 1 && (origin .= xyz)
            base = (sample - 1) * stride
            for k in 1:3
                data[base + k] = xyz[k] - origin[k]
            end
            values, _ = gmsh.view.probe(field, xyz[1] + offsetX,
                xyz[2] + offsetY, xyz[3] + offsetZ, probe_step, -1, false, -1)
            length(values) == nc * count ||
                error("plotOnPath: incomplete probe result at sample $sample.")
            for j in 1:count
                vbase = (j - 1) * nc
                value = if nc == 1
                    values[vbase + 1]
                elseif nc == 3
                    sqrt(sum(abs2, @view values[vbase+1:vbase+3]))
                else
                    a,b,c = values[vbase+1], values[vbase+5], values[vbase+9]
                    sqrt(0.5*((a-b)^2+(b-c)^2+(c-a)^2+
                        6*(values[vbase+2]^2+values[vbase+3]^2+values[vbase+6]^2)))
                end
                data[base + 3 + j] = value
            end
            if plot
                x[sample] = sample == 1 ? 0.0 : x[sample-1] + norm(xyz - previous)
                y[sample] = data[base + 4]
                previous .= xyz
            end
        end
    end
    tag = gmsh.view.add(name)
    gmsh.view.addListData(tag, "SP", nsamples, data)
    gmsh.view.option.setNumber(tag, "Type", 2)
    gmsh.view.option.setNumber(tag, "Axes", 1)
    gmsh.view.option.setNumber(tag, "AdaptVisualizationGrid", 0)
    gmsh.view.option.setNumber(tag, "Visible", visible ? 1 : 0)
    return plot ? (tag, (x, y)) : tag
end

"""
    showOnSurface(field::Number, phName::String; 
                  grad=false, component=:x, offsetX=0, offsetY=0, offsetZ=0, 
                  name=phName, visible=false)

Shows the values of a scalar field on a surface with physical name `phName`.
`field` is the tag of a View in Gmsh. The values of the field are calculated at the
intersection with the surface. `grad` is a Boolean to toggle the gradient of the field on or off.
`component` is the component of the gradient of `field` (:x, :y, :z) to be shown. `offsetX`, `offsetY`, `offsetZ`
are the offsets in the x, y, and z directions where the values are sampled. `name` is a title to display,
and `visible` is a Boolean to toggle the initial visibility in Gmsh on or off.

Returns: `tag`

Types:
- `field`: Integer
- `phName`: String
- `grad`: Boolean
- `component`: Symbol
- `offsetX`: Float64
- `offsetY`: Float64
- `offsetZ`: Float64
- `name`: String
- `visible`: Boolean
- `tag`: Integer
"""
function showOnSurface(field::Number, phName::String; grad=false, component=:x, offsetX=0, offsetY=0, offsetZ=0, name=phName, visible=false)
    SS = gmsh.view.add(name)
    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)
    nodeTags, ncoord, parametricCoord = gmsh.model.mesh.getNodes(-1, -1, true, false)
    ncoord2 = similar(ncoord)
    ncoord2[nodeTags*3 .- 2] .= ncoord[1:3:length(ncoord)]
    ncoord2[nodeTags*3 .- 1] .= ncoord[2:3:length(ncoord)]
    ncoord2[nodeTags*3 .- 0] .= ncoord[3:3:length(ncoord)]
    if component == :x
        comp0 = 1
    elseif component == :y
        comp0 = 2
    elseif component == :z
        comp0 = 3
    else
        error("resultant: invalid component '$component'")
    end
    ret2 = []
    ret3 = []
    ret4 = []
    el2 = 0
    el3 = 0
    el4 = 0
    x = [0.0, 0.0, 0.0]
    for idm in 1:length(dimTags)
        dimTag = dimTags[idm]
        edim = dimTag[1]
        etag = dimTag[2]
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(edim, etag)
        for i in 1:length(elemTypes)
            et = elemTypes[i]
            elementName, dim, order, numNodes::Int64, localNodeCoord, numPrimaryNodes = gmsh.model.mesh.getElementProperties(et)
            nn = zeros(numPrimaryNodes * 4)
            for j in 1:length(elemTags[i])
                elem = elemTags[i][j]
                for l in 1:numPrimaryNodes
                    for k in 1:3
                        x[k] = ncoord2[elemNodeTags[i][j*numNodes-(numNodes-l)]*3-(3-k)]
                        nn[k*numPrimaryNodes-numPrimaryNodes + l] = x[k]
                    end
                    k = 4
                    f, d = gmsh.view.probe(field, x[1]+offsetX, x[2]+offsetY, x[3]+offsetZ, -1, -1, grad, -1)
                    nn[4*numPrimaryNodes-numPrimaryNodes + l] = f[comp0]
                end
                if numPrimaryNodes == 3
                    append!(ret3,nn)
                elseif numPrimaryNodes == 4
                    append!(ret4,nn)
                elseif numPrimaryNodes == 2
                    append!(ret2,nn)
                end
            end
            if numPrimaryNodes == 3
                el3 += length(elemTags[i])
            elseif numPrimaryNodes == 4
                el4 += length(elemTags[i])
            elseif numPrimaryNodes == 2
                el2 += length(elemTags[i])
            end
        end
    end
    if el3 > 0
        gmsh.view.addListData(SS, "ST", el3, ret3)
    end
    if el4 > 0
        gmsh.view.addListData(SS, "SQ", el4, ret4)
    end
    if el2 > 0
        gmsh.view.addListData(SS, "SL", el2, ret2)
    end
    if visible == false
        gmsh.view.option.setNumber(SS, "Visible", 0)
    end
    return SS
end

function showOnSurface(field::ScalarField, phName::String; name=phName)
    s = elementsToNodes(field)
    s = nodesToElements(s, onPhysicalGroup=phName)
    view = showElementResults(s, name=name)
    return view
end

"""
    plotOnBeam(
        phName::String,
        s::ScalarField;
        step=nothing,
        plot=false,
        name="field on \$phName",
        visible=false,
        reverse=false
    )

Create a diagram of a scalar field along a one-dimensional finite element path.

The physical group `phName` must define a single connected, non-branching,
open chain of 1D finite elements. The element order is reconstructed from the
mesh topology, independently of the order in which Gmsh returns the elements.

The horizontal coordinate of the graph is the cumulative length of the
discretized beam path.

Both nodal and elementwise `ScalarField` representations are supported.
For elementwise fields, common element-end nodes are repeated in the graph,
so discontinuities between adjacent elements are preserved.

If the field contains multiple time steps and `step=nothing`, all time steps
are added to the Gmsh view and can be animated. If `step` is specified, it is
interpreted as a zero-based time-step index, consistently with `plotOnPath`.

With `reverse=true`, traversal starts from the opposite end of the beam path.

# Arguments

- `phName`: name of a 1D Gmsh physical group defining the beam path.
- `s`: scalar field to plot.

# Keyword arguments

- `step=nothing`: plot all time steps, or select one zero-based time-step index.
- `plot=false`: if `true`, also return `(x, y)` data for the selected or first step.
- `name="field on \$phName"`: name of the generated Gmsh view.
- `visible=false`: initial visibility of the generated Gmsh view.
- `reverse=false`: reverse the traversal direction.

# Returns

If `plot=false`:

    beamView

If `plot=true`:

    beamView, (x, y)

where `x` is the cumulative beam length and `y` contains the values of the
selected time step, or the first time step if all steps are shown.

# Notes

The physical group must contain exactly one open path. Branching, disconnected
paths and closed loops are rejected because they do not define a unique beam
coordinate.
"""
function plotOnBeam(
    phName::String,
    ss::ScalarField;
    step=nothing,
    plot=false,
    name="field on $phName",
    visible=false,
    reverse=false
    )

    problem = ss.model
    gmsh.model.setCurrent(problem.name)

    s = elementsToElements(ss, onPhysicalGroup=phName)

    # ------------------------------------------------------------------
    # 1) Get the 1D physical group
    # ------------------------------------------------------------------

    dimTags = gmsh.model.getEntitiesForPhysicalName(phName)

    isempty(dimTags) &&
        error("plotOnBeam: physical group \"$phName\" does not exist.")

    for (dim, _) in dimTags
        dim == 1 ||
            error(
                "plotOnBeam: physical group \"$phName\" must contain " *
                "only one-dimensional entities."
            )
    end

    # ------------------------------------------------------------------
    # 2) Collect all finite elements of the physical group
    # ------------------------------------------------------------------

    elementTags  = Int[]
    elementTypes = Int[]
    elementNodes = Vector{Vector{Int}}()

    endNode1 = Int[]
    endNode2 = Int[]

    seenElements = Set{Int}()

    for (_, entityTag) in dimTags

        elemTypes, elemTags, elemNodeTags =
            gmsh.model.mesh.getElements(1, entityTag)

        for itype in eachindex(elemTypes)

            et = Int(elemTypes[itype])

            _, edim, _, numNodes, _, numPrimaryNodes =
                gmsh.model.mesh.getElementProperties(et)

            edim == 1 ||
                error(
                    "plotOnBeam: element type $et is not one-dimensional."
                )

            numPrimaryNodes >= 2 ||
                error(
                    "plotOnBeam: element type $et has fewer than two " *
                    "primary nodes."
                )

            tags = elemTags[itype]
            conn = elemNodeTags[itype]

            length(conn) == length(tags) * numNodes ||
                error(
                    "plotOnBeam: inconsistent Gmsh connectivity for " *
                    "element type $et."
                )

            for ie in eachindex(tags)

                elem = Int(tags[ie])

                # Avoid accidental duplication if entities overlap.
                elem in seenElements && continue
                push!(seenElements, elem)

                first = (ie - 1) * numNodes + 1

                nodes = Int.(
                    conn[first:first + numNodes - 1]
                )

                push!(elementTags, elem)
                push!(elementTypes, et)
                push!(elementNodes, nodes)

                # For Gmsh Lagrange line elements the first two primary
                # nodes are the two element endpoints.
                push!(endNode1, nodes[1])
                push!(endNode2, nodes[2])
            end
        end
    end

    nelem = length(elementTags)

    nelem > 0 ||
        error(
            "plotOnBeam: physical group \"$phName\" contains no 1D mesh elements."
        )

    # ------------------------------------------------------------------
    # 3) Build endpoint topology
    # ------------------------------------------------------------------

    adjacency = Dict{Int,Vector{Int}}()

    for ie in 1:nelem

        n1 = endNode1[ie]
        n2 = endNode2[ie]

        push!(get!(adjacency, n1, Int[]), ie)
        push!(get!(adjacency, n2, Int[]), ie)
    end

    # A beam path cannot contain branching nodes.
    for (node, elems) in adjacency
        if length(elems) > 2
            error(
                "plotOnBeam: physical group \"$phName\" contains a " *
                "branching node (node $node is connected to " *
                "$(length(elems)) beam elements)."
            )
        end
    end

    endNodes = sort([
        node
        for (node, elems) in adjacency
        if length(elems) == 1
    ])

    if length(endNodes) == 0
        error(
            "plotOnBeam: physical group \"$phName\" forms a closed path. " *
            "An open beam path with two endpoints is required."
        )
    elseif length(endNodes) != 2
        error(
            "plotOnBeam: physical group \"$phName\" does not form one " *
            "connected open beam path. Found $(length(endNodes)) endpoints."
        )
    end

    startNode =
        reverse ? endNodes[2] : endNodes[1]

    targetEnd =
        reverse ? endNodes[1] : endNodes[2]

    # ------------------------------------------------------------------
    # 4) Traverse the beam graph
    # ------------------------------------------------------------------

    visited = falses(nelem)

    orderedElements = Int[]
    orderedStartNode = Int[]

    currentNode = startNode

    while true

        candidates = Int[]

        for ie in adjacency[currentNode]
            !visited[ie] && push!(candidates, ie)
        end

        isempty(candidates) && break

        length(candidates) == 1 ||
            error(
                "plotOnBeam: ambiguous traversal at node $currentNode."
            )

        ie = candidates[1]

        push!(orderedElements, ie)
        push!(orderedStartNode, currentNode)

        visited[ie] = true

        if endNode1[ie] == currentNode
            currentNode = endNode2[ie]
        elseif endNode2[ie] == currentNode
            currentNode = endNode1[ie]
        else
            error("plotOnBeam: internal topology error.")
        end
    end

    all(visited) ||
        error(
            "plotOnBeam: physical group \"$phName\" contains disconnected " *
            "beam segments."
        )

    currentNode == targetEnd ||
        error(
            "plotOnBeam: beam traversal did not terminate at the expected endpoint."
        )

    # ------------------------------------------------------------------
    # 5) Select time steps
    # ------------------------------------------------------------------

    s.nsteps >= 1 ||
        error("plotOnBeam: ScalarField contains no time steps.")

    stepRange = Int[]

    if step === nothing

        stepRange = collect(1:s.nsteps)

    else

        step isa Integer ||
            error("plotOnBeam: step must be an integer or nothing.")

        step >= 0 ||
            error("plotOnBeam: step must be non-negative.")

        if step >= s.nsteps

            @warn(
                "plotOnBeam: step is greater than max. number of steps " *
                "(max. number is chosen) $step <==> $(s.nsteps - 1)!"
            )

            stepRange = [s.nsteps]

        else

            # User-facing step numbering follows Gmsh: 0,1,2,...
            stepRange = [Int(step) + 1]
        end
    end

    nplotsteps = length(stepRange)

    # ------------------------------------------------------------------
    # 6) Global node coordinates
    # ------------------------------------------------------------------

    nodeTags, coordinates, _ =
        gmsh.model.mesh.getNodes()

    coord = Dict{Int,NTuple{3,Float64}}()
    sizehint!(coord, length(nodeTags))

    for (i, tag0) in enumerate(nodeTags)

        tag = Int(tag0)

        coord[tag] = (
            Float64(coordinates[3i - 2]),
            Float64(coordinates[3i - 1]),
            Float64(coordinates[3i])
        )
    end

    # ------------------------------------------------------------------
    # 7) Local node ordering for each Gmsh line element type
    # ------------------------------------------------------------------

    localPermutation =
        Dict{Int,Vector{Int}}()

    function get_local_permutation(et::Int, numNodes::Int)

        return get!(localPermutation, et) do

            _, edim, _, nnode, localNodeCoord, _ =
                gmsh.model.mesh.getElementProperties(et)

            edim == 1 ||
                error("plotOnBeam: expected a 1D element.")

            nnode == numNodes ||
                error(
                    "plotOnBeam: inconsistent number of nodes for element type $et."
                )

            stride, rem =
                divrem(length(localNodeCoord), numNodes)

            rem == 0 && stride >= 1 ||
                error(
                    "plotOnBeam: cannot interpret local coordinates " *
                    "of element type $et."
                )

            ξ = Vector{Float64}(undef, numNodes)

            for a in 1:numNodes
                ξ[a] =
                    localNodeCoord[(a - 1) * stride + 1]
            end

            sortperm(ξ)
        end
    end

    # ------------------------------------------------------------------
    # 8) Field lookup
    # ------------------------------------------------------------------

    elementFieldIndex = Dict{Int,Int}()

    if isElementwise(s)

        sizehint!(elementFieldIndex, length(s.numElem))

        for (i, elem) in enumerate(s.numElem)
            elementFieldIndex[Int(elem)] = i
        end

    elseif !isNodal(s)

        error(
            "plotOnBeam: ScalarField is neither nodal nor elementwise."
        )
    end

    # ------------------------------------------------------------------
    # 9) Build ordered graph points
    # ------------------------------------------------------------------

    pointCoordinates =
        NTuple{3,Float64}[]

    pointValues =
        Vector{Vector{Float64}}()

    for (ipath, ie) in enumerate(orderedElements)

        elem = elementTags[ie]
        et = elementTypes[ie]
        nodes = elementNodes[ie]

        numNodes = length(nodes)

        perm =
            copy(get_local_permutation(et, numNodes))

        firstNode = orderedStartNode[ipath]

        # Orient the local element node sequence consistently with
        # the global beam traversal.
        if nodes[perm[1]] == firstNode

            # already correctly oriented

        elseif nodes[perm[end]] == firstNode

            reverse!(perm)

        else

            error(
                "plotOnBeam: could not orient element $elem along the beam path."
            )
        end

        if isElementwise(s)

            fieldIndex =
                get(elementFieldIndex, elem, 0)

            fieldIndex != 0 ||
                error(
                    "plotOnBeam: ScalarField contains no data for element $elem."
                )

            Ae = s.A[fieldIndex]

            # Usually an elementwise ScalarField contains one value at
            # each local element node. A single-row field is also
            # accepted as a constant value over the element.
            constantElement = size(Ae, 1) == 1

            if !constantElement && size(Ae, 1) != numNodes
                error(
                    "plotOnBeam: field data for element $elem has " *
                    "$(size(Ae,1)) rows, but the element has $numNodes nodes."
                )
            end

            size(Ae, 2) >= maximum(stepRange) ||
                error(
                    "plotOnBeam: field data for element $elem contains " *
                    "too few time steps."
                )

            # Keep both endpoints of every element. At a common node
            # this intentionally creates two points at the same beam
            # coordinate, preserving jumps in elementwise quantities.
            for localIndex in perm

                node = nodes[localIndex]

                xyz =
                    get(coord, node, nothing)

                xyz === nothing &&
                    error(
                        "plotOnBeam: coordinates of node $node were not found."
                    )

                row =
                    constantElement ? 1 : localIndex

                values =
                    Float64[
                        Ae[row, it]
                        for it in stepRange
                    ]

                push!(pointCoordinates, xyz)
                push!(pointValues, values)
            end

        else

            # Nodal quantities are continuous in the FE representation.
            # Skip the first node of every element after the first one
            # to avoid duplicating common endpoints.
            firstLocal =
                ipath == 1 ? 1 : 2

            for q in firstLocal:length(perm)

                localIndex = perm[q]
                node = nodes[localIndex]

                node <= size(s.a, 1) ||
                    error(
                        "plotOnBeam: nodal field does not contain node $node."
                    )

                xyz =
                    get(coord, node, nothing)

                xyz === nothing &&
                    error(
                        "plotOnBeam: coordinates of node $node were not found."
                    )

                values =
                    Float64[
                        s.a[node, it]
                        for it in stepRange
                    ]

                push!(pointCoordinates, xyz)
                push!(pointValues, values)
            end
        end
    end

    isempty(pointCoordinates) &&
        error("plotOnBeam: no plot points were generated.")

    # ------------------------------------------------------------------
    # 10) Convert to Gmsh list-data representation
    # ------------------------------------------------------------------

    p0 = pointCoordinates[1]

    CoordValue = Float64[]

    recordLength = 3 + nplotsteps

    sizehint!(
        CoordValue,
        recordLength * length(pointCoordinates)
    )

    for i in eachindex(pointCoordinates)

        p = pointCoordinates[i]

        push!(CoordValue, p[1] - p0[1])
        push!(CoordValue, p[2] - p0[2])
        push!(CoordValue, p[3] - p0[3])

        append!(CoordValue, pointValues[i])
    end

    # ------------------------------------------------------------------
    # 11) Create Gmsh graph view
    # ------------------------------------------------------------------

    beamView = gmsh.view.add(name)

    gmsh.view.addListData(
        beamView,
        "SP",
        length(pointCoordinates),
        CoordValue
    )

    gmsh.view.option.setNumber(
        beamView,
        "Type",
        2
    )

    gmsh.view.option.setNumber(
        beamView,
        "Axes",
        1
    )

    gmsh.view.option.setNumber(
        beamView,
        "AdaptVisualizationGrid",
        0
    )

    if !visible
        gmsh.view.option.setNumber(
            beamView,
            "Visible",
            0
        )
    end

    # ------------------------------------------------------------------
    # 12) Optional x-y data
    # ------------------------------------------------------------------

    if plot

        np = length(pointCoordinates)

        x = zeros(Float64, np)
        y = zeros(Float64, np)

        y[1] = pointValues[1][1]

        for i in 2:np

            p0i = pointCoordinates[i - 1]
            p1i = pointCoordinates[i]

            dx = p1i[1] - p0i[1]
            dy = p1i[2] - p0i[2]
            dz = p1i[3] - p0i[3]

            x[i] =
                x[i - 1] +
                sqrt(dx^2 + dy^2 + dz^2)

            y[i] = pointValues[i][1]
        end

        return beamView, (x, y)

    else

        return beamView
    end
end

