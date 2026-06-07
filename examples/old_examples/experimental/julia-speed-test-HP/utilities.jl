import Gmsh: gmsh

#==========================================================================================================================
                                                I. STRUCTURES
===========================================================================================================================#

# A structure containing the material properties:
mutable struct Material
    name::String  # name of the material
    E::Float64    # Young's modulus
    ν::Float64    # Poisson's ratio
    G::Float64    # shear modulus
    ρ::Float64    # mass density
    k::Float64    # heat conductivity
    c::Float64    # specific heat
    α::Float64    # heat expansion coefficient

    # constructing default material (Steel):
    function Material(name = "mat1"; E=2.0e5, ν=0.3, ρ=7.85e-9, k=45, c=4.2e8, α=1.2e-5)
        G = E/(2*(1+ν))
        return new(name, E, ν, G, ρ, k, c, α)
    end
end

# A structure containig Solid3D section properties:
mutable struct Solid3D 
    scope::Any         # scope of the section
    mat::Material      # material of the section
    dim::Int64         # spatial dimension of the section
    pdim::Int64        # physical dimension of the section
    D::Matrix{Float64} # material matrix
    intid::Int64       # integration rule identifier

    # constructing default Solid section:
    function Solid3D(mat=Material(); scope=:v, dim=3, pdim=3, D=[], intid=1)
        if isempty(D)
            E = mat.E
            ν = mat.ν
            D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0 0 0;
                                            ν 1-ν ν 0 0 0;
                                            ν ν 1-ν 0 0 0;
                                            0 0 0 (1-2ν)/2 0 0;
                                            0 0 0 0 (1-2ν)/2 0;
                                            0 0 0 0 0 (1-2ν)/2]  
        end
        return new(scope, mat, dim, pdim, D, intid)
    end
end

# A structure containig Solid2D section properties:
mutable struct Solid2D
    scope::Any         # scope of the section
    type::Symbol       # type of the 2D Solid
    mat::Material      # material of the section
    dim::Int64         # spatial dimension of the section
    pdim::Int64        # physical dimension of the section
    width::Float64     # width of the section (in case of Plane Stress)
    D::Matrix{Float64} # material matrix
    intid::Int64       # integration rule identifier

    # constructing default 2D Solid section:
    function Solid2D(type=:PlaneStress, mat=Material(); scope=:a, dim=2, pdim=2, width=1, D=[], intid=1)
        if isempty(D)
            E = mat.E
            ν = mat.ν
            if type == :PlaneStrain  # 2D planestarin element
                D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν 0;
                                                ν 1-ν 0;
                                                0 0 (1-2ν)/2]   

            elseif type == :PlaneStress # 2D planestress element
                D = E / (1 - ν^2) * [1 ν 0;
                                    ν 1 0;
                                    0 0 (1-ν)/2]   

            elseif type == :AxiSymmetricX || type == :AxiSymmetricY  # axisymmetric element    
                D = E / ((1 + ν) * (1 - 2ν)) * [1-ν ν ν 0;
                                                ν 1-ν ν 0;
                                                ν ν 1-ν 0;
                                                0 0 0 (1-2ν)/2]   
            end  
        end          
        return new(scope, type, mat, dim, pdim, width, D, intid)
    end
end

# A structure containig Shell section properties:
mutable struct Shell
    scope::Any         # scope of the section
    mat::Material      # material of the section
    dim::Int64         # spatial dimension of the section
    pdim::Int64        # physical dimension of the section
    width::Float64     # width of the section
    ks::Float64        # shear deformation correction
    kt::Float64        # torsional defomration factor
    D::Matrix{Float64} # material matrix
    intid::Int64       # integration rule identifier
    intm::Int64        # integration rule identifier for mebrane part
    ints::Int64        # integration rule identifier for shear part
    intt::Int64        # integration rule identifier for torsional part

    # constructing default Shell section:
    function Shell(mat=Material();scope=:a, dim=3, pdim=6, width=1, ks=5/6, kt=0.006, D=[], intid=1, intm=1, ints=1, intt=1)
        if isempty(D)
            E = mat.E
            ν = mat.ν
            D = (E/(1-ν^2))*[1 ν 0       0          0;
                             ν 1 0       0          0;
                             0 0 (1-ν)/2 0          0;
                             0 0 0       ks*(1-ν)/2 0;
                             0 0 0       0          ks*(1-ν)/2]
        end         

        return new(scope, mat, dim, pdim, width, ks, kt, D, intid, intm, ints, intt)
    end
end

# A structure containig PointMass section properties:
mutable struct PointMass
    scope::Any           # scope of the section
    dim::Int64           # spatial dimension of the section
    pdim::Int64          # physical dimension of the section
    m::Float64           # mass
    mx::Float64          # mass in x direction 
    my::Float64          # mass in y direction 
    mz::Float64          # mass in z direction 
    Jx::Float64          # moment of inertia about axis x 
    Jy::Float64          # moment of inertia about axis y 
    Jz::Float64          # moment of inertia about axis z 
    intid::Int64         # integration rule identifier

    # constructing default PintMass section:
    function PointMass(scope="pilot1", dim=3, pdim=6; m=0, mx=0, my=0, mz=0, Jx=0, Jy=0, Jz=0, intid=1)
        if m != 0
            mx = my = mz = m
        end
        return new(scope, dim, pdim, m, mx, my, mz, Jx, Jy, Jz, intid)
    end
end

# A structure containig element shape, type, default GaussPoint properties
struct Gauss
    shape::Symbol
    type::Symbol
    GaussString::String
end

# A structure containing element properties:
struct Element
    name::String                 # Name of the element
    dim::Int32                   # Dimension of the element (0,1,2 or 3)
    order::Int32                 # Order of the element
    non::Int32                   # Number of nodes
    lnc::Vector{Float64}         # Local (ξηζ) nodal coordinates
    nopn::Int32                  # Number of primary nodes
    tags::Vector{Int32}          # Element tag vector
    nn::Matrix{Int32}            # Connectivity matrix
    GP::Vector{Matrix{Float64}}  # Gauss Point coordinates in ξηζ csys
    GW::Vector{Vector{Float64}}  # Gauss Weigths
    HG::Vector{Matrix{Float64}}  # Shape function values at Gauss points
    HξG::Vector{Matrix{Float64}} # Shape function derivatives with respect to ξ
    HηG::Vector{Matrix{Float64}} # Shape function derivatives with respect to η
    HζG::Vector{Matrix{Float64}} # Shape function derivatives with respect to ζ
    HN::Matrix{Float64}          # Shape function values at Nodes
    HξN::Matrix{Float64}         # Shape function derivatives with respect to ξ
    HηN::Matrix{Float64}         # Shape function derivatives with respect to η
    HζN::Matrix{Float64}         # Shape function derivatives with respect to ζ
    CP::Matrix{Float64}          # center point xyz coordinates        
end

# A structure containing node properties:
struct Node
    xyz::Matrix           # Coordinate matrix of the nodes
    name::String          # name of the physical group of the nodes
    nTags::Vector{Int64}  # Vector of the node tags
    phTag::Int64          # tag of the Physical group containing the nodes
    function Node(xyz,name="nodeset1";nTags=[0],phTag=0)
        n = size(xyz,1)  # number of defined nodes
        if size(xyz,2) == 2
            xyz = [xyz zeros(n,1)]
        end
        nTags isa Array ? nTags = nTags : nTags = [nTags]
        return new(xyz,name,nTags,phTag)
    end
end

# A structure containing selection properties:
struct Selection
    name::String                      # Name of the selection
    dim::Int32                        # Dimension of the selection
    id::Int32                         # Id of the selection
    phnTags::Vector{Int32}            # Node tags of the selection
    elems::Dict{Int16, Vector{Int32}} # Element types and tags of the selection
end

# A structure containing selector properties:
struct Selector
    object::Symbol
    seldef::Symbol
    selpar::Any
    presel::Any
    id::Int32
    name::String
    tol::Float64
end

# A structure describing a finite element mesh:
mutable struct Mesh
    xyz::Matrix{Float64}             # Coordinate matrix
    elem::Dict{Int16, Element}       # Elemet property dictionary
    set::Dict{Tuple{Int32, Int32}, Selection}  # sets
    findset::Dict{Any, Tuple{Int32, Int32}}    # set dictionary
    nTags::Vector{Int32}             # Node tags
    non::Int64                       # number of nodes
    name::String

    function Mesh(args...;approx="Lagrange", intrule=:default, renumber=:none, selready=true, addnodes=:none)
        xyz = Matrix{Float64}(undef, 0, 0)  # initialization of coordinate matrix
        elem = Dict{Int16, Element}()       # initialization of elemet property dictionary
        set = Dict{Tuple{Int32, Int32}, Selection}() # dictionary of sets
        findset = Dict{Any, Tuple{Int32, Int32}}()   # set dictionary

        nargs = length(args)  # number of arguments
        if endswith(args[1], ".geo") || args[1] == "selgeo" || args[1] == "Gmsh" # Gmsh mesh
            
            # Opening Gmsh file ================================================

            gmsh.initialize()  # initialization of gmsh
            gmsh.option.setNumber("General.Verbosity", 0)  # turn off verbosity

            # name = gmsh.model.getCurrent()

            if endswith(args[1], ".geo")
                gmsh.open(args[1])
            elseif endswith(args[1], "selgeo")
                selected_file = open_dialog("Pick a *.geo file", GtkNullContainer(), ("*.geo",), select_multiple=false)
                if selected_file == "" || selected_file === nothing
                    return  
                end
                gmsh.open(selected_file)
            end

            name = gmsh.model.getCurrent()  # name of the gmsh file

            if renumber != :none
                # Reducing the bandwith by node renumbering:
                oldTags, newTags = gmsh.model.mesh.computeRenumbering(renumber)
                gmsh.model.mesh.renumberNodes(oldTags, newTags)
            end

            # adding extra nodes if requested:
            if addnodes != :none
                addnodes isa Array && !isempty(addnodes) ? addnodes = addnodes : addnodes = [addnodes]
                if !isempty(addnodes)  # if additional nodes are provided

                    for (nodeKey, nodeData) in pairs(addnodes)  # for each nodeset
                        xyz    = nodeData.xyz     # xyz coordinates of the current node set
                        phName   = nodeData.name  # physical group name of the current node set
                        nTags = nodeData.nTags    # node tags of the current nodes
                        phTag  = nodeData.phTag   # physical group tag of the current node set

                        n     = size(xyz,1)       # number of nodes in the current node set

                        if nTags[1] == 0          # no node tags were provided
                            mnt = get_max_node_tag()      # get the maximum existing node tag
                            nTags = collect(mnt+1:mnt+n)  # creating the node tags vector
                        end

                        if phTag == 0             # no physical group tag was provided
                            mpht = get_max_physical_group_tag()  # get maximum existing physical group tag
                            phTag = mpht+1        # physical group tag of the current nodeset (max tag + 1)
                        end

                        ment = get_max_entity_tag(0)     # largest point entity tag
                        enTags = collect(ment+1:ment+n)   # vector of entity tags

                        met = get_max_element_tag()  # largest element tag

                        for i = 1:n
                            enTag = enTags[i]  # entity tag of the ith point
                            nTag  = [nTags[i]]   # node tag of ith node
                            xyzi = xyz[i,:]    # coordinates of the ith node
                            gmsh.model.addDiscreteEntity(0, enTag) # add a discrete 0D entity (which will hold the extra node)
                            gmsh.model.mesh.addNodes(0, enTag, nTag, xyzi)  # assign node to the ith point entity
                            gmsh.model.mesh.addElements(0, enTag, [15], [[met+i]], [[nTag[1]]])
                        end

                        gmsh.model.addPhysicalGroup(0, enTags, phTag)  # assign physical group to the current node set
                        gmsh.model.setPhysicalName(0, phTag, phName)   # assign name to the physical group
                    end
                end
            end

            # Coordinate Matrix ================================================

            nTags, nCoords, _ = gmsh.model.mesh.getNodes()
            xyz = reshape(nCoords, 3, length(nTags))' # Coordinate matrix
            xyz[nTags, :] .= round.(xyz, digits = 8)

            # Element Dictionary ===============================================

            # Obtaining elemet Types, Element Tags and Node tags:
            eTypes, eTags, enTags = gmsh.model.mesh.getElements()
            elem = Dict()  # initialization of Element Dictionary

            for i in eachindex(eTypes)     # for all element types
                eTypei = eTypes[i]         # element Type number of the ith element type
                eTagsi = Int64.(eTags[i])  # element Tags of the ith element type

                # Obtaining element properties for the ith element type:
                namei, dimi, ordi, noni, lnci, pni = gmsh.model.mesh.getElementProperties(eTypei)

                nodeCoord = zeros(noni * 3)  # local nodal coordinate vector
                for k in 1:dimi, j = 1:noni
                    nodeCoord[k+(j-1)*3] = lnci[k+(j-1)*dimi]
                end

                # connectivity matrix of the ith element type
                nni = Int32.(reshape(Int64.(enTags[i]), Int64(noni), length(eTagsi)))'

                # Gauss Points and Weights:
                shape = getIntprop(eTypei).shape
                type  = getIntprop(eTypei).type

                if intrule isa NamedTuple && haskey(intrule, type)
                    rules = to_vector(intrule[type])
                    nir = length(rules)   # number of integration rules
                    GaussStringVec = Vector{String}(undef, nir)
                    for j = 1:nir
                        rule = rules[j]
                        if rule isa Integer
                            shape = Symbol(match(r"^\D+", String(type)).match)
                            GaussStringVec[j] = getIntprop(shape,rule)
                        elseif rule == :default
                            GaussStringVec[j] = getIntprop(eTypei).GaussString
                        elseif rule isa Symbol
                            GaussStringVec[j] = getIntprop(rule)
                        end
                    end
                elseif intrule isa NamedTuple && haskey(intrule, shape)
                    rules = to_vector(intrule[shape])
                    nir = length(rules)   # number of integration rules
                    GaussStringVec = Vector{String}(undef, nir)
                    for j = 1:nir
                        rule = rules[j]
                        if rule isa Integer
                            GaussStringVec[j] = getIntprop(shape,rule)
                        elseif rule == :default
                            GaussStringVec[j] = getIntprop(eTypei).GaussString
                        elseif rule isa Symbol
                            GaussStringVec[j] = getIntprop(rule)
                        end
                    end
                else                    
                    nir = 1   # number of integration rules
                    GaussStringVec = [getIntprop(eTypei).GaussString] # string to add to GMSH
                end

                GP  = [Array{Float64}(undef, 0, 0) for _ in 1:nir]
                GW  = [Vector{Float64}(undef, 0) for _ in 1:nir]
                HG  = [Array{Float64}(undef, 0, 0) for _ in 1:nir]
                HξG = [Array{Float64}(undef, 0, 0) for _ in 1:nir]
                HηG = [Array{Float64}(undef, 0, 0) for _ in 1:nir]
                HζG = [Array{Float64}(undef, 0, 0) for _ in 1:nir]

                for (j, GaussString) in enumerate(GaussStringVec)
                    GPj, GWj = gmsh.model.mesh.getIntegrationPoints(eTypei,GaussString)
                    nGP = Int64.(length(GWj))  # number of Gauss Points
                    GP[j] = reshape(GPj,3,nGP)'  # Gauss Point matrix
                    GW[j] = GWj
                    # Shape functions at Gauss points:
                    _, basisf, _ = gmsh.model.mesh.getBasisFunctions(eTypei, GPj, approx)
                    HG[j] = reshape(basisf,Int64(noni),nGP)
                    # Shape function derivatives at Gauss points:
                    _, gradbasisf, _ = gmsh.model.mesh.getBasisFunctions(eTypei, GPj, "Grad$approx")
                    HξG[j] = reshape(gradbasisf[1:3:3*noni*nGP-2],Int64(noni),nGP)
                    HηG[j] = reshape(gradbasisf[2:3:3*noni*nGP-1],Int64(noni),nGP)
                    HζG[j] = reshape(gradbasisf[3:3:3*noni*nGP  ],Int64(noni),nGP)   
                end                 

                # Shape functions at nodes:
                _, basisf, _ = gmsh.model.mesh.getBasisFunctions(eTypei, nodeCoord, approx)
                HN = reshape(basisf,Int64(noni),Int64(noni))

                # Shape function derivatives at nodes:
                _, gradbasisf, _ = gmsh.model.mesh.getBasisFunctions(eTypei, nodeCoord, "Grad$approx")
                HξN = reshape(gradbasisf[1:3:3*noni*noni-2],Int64(noni),Int64(noni))
                HηN = reshape(gradbasisf[2:3:3*noni*noni-1],Int64(noni),Int64(noni))
                HζN = reshape(gradbasisf[3:3:3*noni*noni  ],Int64(noni),Int64(noni))

                if selready # if selection capability (center point coordinates) is required:
                    # center point in ξηζ coordinate system:
                    GPC, _ = gmsh.model.mesh.getIntegrationPoints(eTypei,"Gauss0")
                    # center point in xyz coordinate system:
                    _, _, coords = gmsh.model.mesh.getJacobians(eTypei,GPC)
                    CP = reshape(coords, 3, length(eTagsi))' # Coordinate matrix
                    CP .= round.(CP, digits = 8)
                else
                    CP = Matrix{Float64}(undef, 0, 0)
                end

                # Creating Element structure for the ith element:
                elem[eTypei] = Element(namei,dimi,ordi,noni,lnci,pni,eTagsi,nni,GP,GW,HG,HξG,HηG,HζG,HN,HξN,HηN,HζN,CP)
            end

            # Physical group Dictionary ========================================
            
            findset = Dict{Any, Tuple{Int32, Int32}}()   # initialization of set finder Dictionary
            set = Dict{Tuple{Int32, Int32}, Selection}() # initialization of set Dictionary
            
            # creating all sets:
            dimPhtags = gmsh.model.getPhysicalGroups()  # obtainig all Physical groups
            for i in eachindex(dimPhtags) # for each Physaical group
                dim = dimPhtags[i][1]     # dimension of the ith Physical group
                phTag = dimPhtags[i][2]   # Physicsal tag of the ith Physical group
                phName  = gmsh.model.getPhysicalName(dim,phTag)  # name of the ith Physical group
                phnTags,_ = gmsh.model.mesh.getNodesForPhysicalGroup(dim,phTag)  # node Tags of the ith Physical group
                enTags = gmsh.model.getEntitiesForPhysicalGroup(dim,phTag)  # entity tags for the ith Physical group
                etags = Dict{Int16, Vector{Int32}}()  # initialization of element tags Dictionary
                for j in eachindex(enTags) # for each entity belonging tho the ith Physical group
                    eTypes, eTags, _ = gmsh.model.mesh.getElements(dim,enTags[j]) # element types and tags of the jth entity
                    eTags = [Int64.(v) for v in eTags]  # converting element tags to Int format
                    for k in eachindex(eTags)           # for each element type
                        append!(get!(etags, eTypes[k], Int[]), eTags[k]) # add etags to kth etype
                    end
                end
                dim = Int32(dim)                        # convert dim to Int32
                phnTags = Int32.(phnTags)               # convert phnTags to Int32
                findset[(dim,phTag)] = (dim,phTag)      # creating set finder Dictionary
                findset[phTag] = (dim,phTag)            # creating set finder Dictionary
                findset[phName] = (dim,phTag)           # creating set finder Dictionary
                set[(dim,phTag)] = Selection(phName,dim,phTag,phnTags,etags) # adding the ith set to set Dictionary              
            end
            
            non = size(xyz,1) # number of nodes
            # gmsh.finalize()
        end
        return new(xyz, elem, set, findset, nTags, non, name)
    end
end

# A structure containing displacement definitions:
mutable struct Disp
    scope::Vector{Any}      # vector of scoped geometries
    predef::Vector{Symbol}  # vector of predefined displacements
    dir::Vector{Any}        # vecotr of direction definition
    ux::Vector{Any}         # vecotr of displacement in x direction
    uy::Vector{Any}         # vecotr of displacement in y direction
    uz::Vector{Any}         # vecotr of displacement in z direction
    rx::Vector{Any}         # vecotr of rotation about x axis
    ry::Vector{Any}         # vecotr of rotation about y axis  
    rz::Vector{Any}         # vecotr of rotation about z axis
    td::Vector{Any}         # vecotr of time dependency

    # constructing default displacement:
    function Disp(scope, predef=:none; dir=:gsys, ux=:free, uy=:free, uz=:free, rx=:free, ry=:free, rz=:free, td=:none)
        
        multisel = any(x -> x isa Vector, (predef, dir, ux, uy, uz, rx, ry, rz, td))

        scope = scopeProcessor(scope,multisel)
 
        n = length(scope)
        dir = to_vector2(dir,n)
        predef = to_vector2(predef,n)
        ux = to_vector2(ux,n)
        uy = to_vector2(uy,n)
        uz = to_vector2(uz,n)
        rx = to_vector2(rx,n)
        ry = to_vector2(ry,n)
        rz = to_vector2(rz,n)
        td = to_vector2(td,n)

        return new(scope,predef,dir,ux,uy,uz,rx,ry,rz,td)

    end
end

# A structure containing loading definitions:
mutable struct Load
    scope::Vector{Any}      # vector of scoped geometries
    dir::Vector{Symbol}     # vecotor of direction definition
    type::Vector{Symbol}    # type of the loads (e.g. distributed/concentrated)
    fx::Vector{Any}         # vecotor of loads in x direction
    fy::Vector{Any}         # vecotor of loads in y direction
    fz::Vector{Any}         # vecotor of loads in z direction
    mx::Vector{Any}         # moments about x axis
    my::Vector{Any}         # moments about y axis
    mz::Vector{Any}         # moments about z axis
    td::Vector{Any}         # vecotor of time dependencies
    width::Vector{Float64}  # vector of width of the loaded surfaces

    # constructing default loading:
    function Load(scope; dir=:gsys, type=:default, fx=0, fy=0, fz=0, mx=0, my=0, mz=0, td=:none, width=1)

        multisel = any(x -> x isa Vector, (dir, type, fx, fy, fz, mx, my, mz, td, width))

        scope = scopeProcessor(scope,multisel)

        n = length(scope)
        dir = to_vector2(dir,n)
        type = to_vector2(type,n)
        fx = to_vector2(fx,n)
        fy = to_vector2(fy,n)
        fz = to_vector2(fz,n)
        mx = to_vector2(mx,n)
        my = to_vector2(my,n)
        mz = to_vector2(mz,n)
        td = to_vector2(td,n)
        width = to_vector2(width,n)

        return new(scope,dir,type,fx,fy,fz,mx,my,mz,td,width)

    end
end

# A structure containing elastic support definitions:
mutable struct Esup
    scope::Vector{Any}      # vector of scoped geometries
    dir::Vector{Symbol}     # vecotor of direction definition
    k::Vector{Any}          # stiffness of the elastic support
    c::Vector{Any}          # damping of the elastic support
    ux::Vector{Any}         # vecotor of loads in x direction
    uy::Vector{Any}         # vecotor of loads in y direction
    uz::Vector{Any}         # vecotor of loads in z direction
    rx::Vector{Any}         # moments about x axis
    ry::Vector{Any}         # moments about y axis
    rz::Vector{Any}         # moments about z axis
    td::Vector{Any}         # vecotor of time dependencies
    width::Vector{Float64}  # vector of width of the supported surfaces   

    function Esup(scope; dir=:gsys, ux=:free, uy=:free, uz=:free, rx=:free, ry=:free, rz=:free, td=:none, width=1)

        multisel = any(x -> x isa Vector, (dir, ux, uy, uz, rx, ry, rz, td, width))

        scope = scopeProcessor(scope,multisel)

        n = length(scope)
        dir = to_vector2(dir,n)
        ux = to_vector2(ux,n)
        uy = to_vector2(uy,n)
        uz = to_vector2(uz,n)
        rx = to_vector2(rx,n)
        ry = to_vector2(ry,n)
        rz = to_vector2(rz,n)
        td = to_vector2(td,n)   
        width = to_vector2(width,n)  
        
        return new(scope,dir,ux,uy,uz,rx,ry,rz,td,width)
    end  

end

# A structure containing initial displacement definitions:
mutable struct IniDisp
    scope::Vector{Any}      # vector of scoped geometries
    dir::Vector{Any}        # vector of direction definition
    ux::Vector{Any}         # vector of displacement in x direction
    uy::Vector{Any}         # vector of displacement in y direction
    uz::Vector{Any}         # vector of displacement in z direction
    rx::Vector{Any}         # vector of rotation about x axis
    ry::Vector{Any}         # vector of rotation about y axis  
    rz::Vector{Any}         # vector of rotation about z axis

    # constructing default initial isplacement:
    function IniDisp(scope; dir=:gsys, ux=:free, uy=:free, uz=:free, rx=:free, ry=:free, rz=:free)

        multisel = any(x -> x isa Vector, (dir, ux, uy, uz, rx, ry, rz))

        scope = scopeProcessor(scope,multisel)

        n = length(scope)
        dir = to_vector2(dir,n)
        ux = to_vector2(ux,n)
        uy = to_vector2(uy,n)
        uz = to_vector2(uz,n)
        rx = to_vector2(rx,n)
        ry = to_vector2(ry,n)
        rz = to_vector2(rz,n)      

        return new(scope,dir,ux,uy,uz,rx,ry,rz)
    end
end

# A structure containing initial velocity definitions:
mutable struct IniVelo
    scope::Vector{Any}      # vector of scoped geometries
    dir::Vector{Any}        # vector of direction definition
    vx::Vector{Any}         # vector of velocity in x direction
    vy::Vector{Any}         # vector of velocity in y direction
    vz::Vector{Any}         # vector of velocity in z direction
    wx::Vector{Any}         # vector of angular velocity about x axis
    wy::Vector{Any}         # vector of angular velocity y axis  
    wz::Vector{Any}         # vector of angular velocity z axis

    # constructing default initial isplacement:
    function IniVelo(scope; dir=:gsys, ux=:free, uy=:free, uz=:free, rx=:free, ry=:free, rz=:free)

        multisel = any(x -> x isa Vector, (dir, vx, vy, vz, wx, wy, wz))

        scope = scopeProcessor(scope,multisel)

        n = length(scope)
        dir = to_vector2(dir,n)
        vx = to_vector2(vx,n)
        vy = to_vector2(vy,n)
        vz = to_vector2(vz,n)
        wx = to_vector2(wx,n)
        wy = to_vector2(wy,n)
        wz = to_vector2(wz,n)      

        return new(scope,dir,vx,vy,vz,wrx,wy,wz)
    end
end

# A structure containing multi point constraint definitions:
mutable struct MPC
    scope1::Any      # vector of 1st scoped geometry
    scope2::Any      # vector of 2nd scoped geometry

    pairdef::Any     # vector of pairing definitions
    contype::Any     # vector of connection types

    show::Any        # vector of mpc visibility

    # constructing default mpc struct:
    function MPC(scope1,scope2;pairdef=:default,contype=:rigid,show=true)

         scope1 = scopeProcessor(scope1,true)
         scope2 = scopeProcessor(scope2,true)

      #  scope1 = to_vector(scope1)
      #  scope2 = to_vector(scope2)

        n1 = length(scope1)
        n2 = length(scope2)

        n1 == 1 ? n = n2 : n = n1

        pairdef = to_vector2(pairdef,n)
        contype = to_vector2(contype,n)
        show = to_vector2(show,n)
        
        return new(scope1, scope2, pairdef, contype, show)
    end
end      

# A structure containing spring definitions:
mutable struct Spring
    scope1::Any      # vector of 1st scoped geometry
    scope2::Any      # vector of 2nd scoped geometry

    k::Any           # vector of vector of stiffnesses
    c::Any           # vector of vector of dampings

    pairdef::Any     # vector of pairing definitions
    orientation::Any # vector of orientations

    # constructing default spring struct:
    function Spring(scope1=:none,scope2=:none;k,c=[],pairdef=:default,orientation=:default)

        scope1 = scopeProcessor(scope1,true)
        scope2 = scopeProcessor(scope2,true)

      #  scope1 = to_vector(scope1)
      #  scope2 = to_vector(scope2)

        n1 = length(scope1)
        n2 = length(scope2)

        n1 == 1 ? n = n2 : n = n1

        scope1 = to_vector2(scope1,n)
        scope2 = to_vector2(scope2,n)
        k = to_vectorvector(k,n)
        c = to_vectorvector(c,n)
        pairdef = to_vector2(pairdef,n)
        orientation = to_vector3(orientation,n)
        
        return new(scope1, scope2, k, c, pairdef, orientation)
    end
end  

# A structure containing structural damping definitions:
mutable struct Damping
    type::Symbol
    alpha::Float64
    beta::Any
    xi::Any
    zeta::Any
    maxmode::Any

    function Damping(type=:Raylight,alpha=0.1,beta=0.1,zeta=[],maxmode=0)
        return new(type,alpha,beta,zeta,maxmode)
    end
end

#=
# A structure containing analysis step properties:
mutable struct Step
    name::String          # name of the step
    type::Symbol          # type of the step (e.g. Static/Modal/Dynamic)
    dura::Float64         # duration of the step
    dt::Float64           # timestep of the step
    tsamp::Int32          # sampling rate of the step
    deco::Symbol          # decomposition applied in the step
    solver::Symbol        # solver type for the step (e.g. Newmark/CDM0...)
    disp::Dict{Any, Disp} # displacement definitions active in the step
    load::Dict{Any, Load} # loading definitions active in the step
    field::Vector{Symbol} # requested fields of the step

    # constructing default step (a static step with 1 timestep)
    function Step(type=:Static;name="step1",dura=1,dt=0.00001,tsamp=100,
                  deco=:none,solver=:Newmark,
                  disp = Dict{Any, Disp}(),load = Dict{Any, Load}(),field=[:q,:v,:a])
            return new(name,type,dura,dt,tsamp,deco,solver,disp,load,field)
    end
end

# A structure containing Model properties:
mutable struct Model
    name::String             # name of the model
    sec::Dict{Any, Section}  # sections of the model
    mesh::Mesh               # finite element mesh of the model
    disp::Dict{Any, Disp}    # displacemenet definitions of the model
    load::Dict{Any, Load}    # loading definitions of the model
    step::Dict{Any, Step}    # step definitions of the model

    # constructing default model (Solid):
    function Model()
        name = "Model1"               # default name of the model (Model1)
        sec  = Dict(1 => Section())   # initialization of section dictionary
        mesh = Mesh("")               # an empty mesh
        disp = Dict{Any, Disp}()      # initialization of displacement definitions dictionary
        load = Dict{Any, Load}()      # initialization of loading definitions dictionary
        step = Dict(1 => Step())      # initialization of step dictionary  

        return new(name, sec, mesh, disp, load, step)
    end
end
=#

# A structure containing system matrices:
mutable struct SystemMatrices
    Ke::SparseMatrixCSC{Float64, Int}  # elastic support stiffness matrix
    Ce::SparseMatrixCSC{Float64, Int}  # elastic support damping matrix
    Ks::SparseMatrixCSC{Float64, Int}  # spring stiffness matrix
    Cs::SparseMatrixCSC{Float64, Int}  # spring damping matrix
    C::SparseMatrixCSC{Float64, Int}   # structural damping matrix
    G::SparseMatrixCSC{Float64, Int}   # MPC matrix
    
    # constructing default (empty) system matrices:
    function SystemMatrices()
        Ke     = spzeros(Float64, 0, 0)
        Ce     = spzeros(Float64, 0, 0)
        Ks     = spzeros(Float64, 0, 0)
        Cs     = spzeros(Float64, 0, 0)
        C      = spzeros(Float64, 0, 0)
        G      = spzeros(Float64, 0, 0)
        return new(Ke,Ce,Ks,Cs,C,G)
    end
end

# A structure containing boundary condition definitions:
mutable struct BoundaryConditions
    Xset::Vector{Int32}     # set of fixed dofs
    Kset::Vector{Int32}     # set of kinematically loaded dofs
    Nset::Vector{Int32}     # set of free dofs
    fconst::Vector{Float64} # constant nodal load vector
    fbaseA::Matrix{Float64} # base nodal load vector for function type amplitude curve
    fvarA::Vector{Function} # function type amplitude curves
    fbaseB::Matrix{Float64} # base nodal load vector for numeric type amplitude curve
    fvarB::Matrix{Float64}  # numeric type amplitude curve  
    qconst::Vector{Float64} # constant nodal displacement vector
    qbaseA::Matrix{Float64} # base nodal displacement vector for function type amplitude curve
    qvarA::Vector{Function} # function type amplitude curves
    qbaseB::Matrix{Float64} # base nodal displacement vector for numeric type amplitude curve
    qvarB::Matrix{Float64}  # numeric type amplitude curve  
    
    # constructing default (empty) boundary conditions:
    function BoundaryConditions()
        Xset   = Vector{Int32}(undef, 0)
        Kset   = Vector{Int32}(undef, 0) 
        Nset   = Vector{Int32}(undef, 0)
        fconst = Vector{Float64}(undef, 0)
        qconst = Vector{Float64}(undef, 0)
        fbaseA = Matrix{Float64}(undef, 0, 0)
        qbaseA = Matrix{Float64}(undef, 0, 0)
        fbaseB = Matrix{Float64}(undef, 0, 0)
        qbaseB = Matrix{Float64}(undef, 0, 0)
        fvarB  = Matrix{Float64}(undef, 0, 0)
        qvarB  = Matrix{Float64}(undef, 0, 0)
        fvarA  = Vector{Function}(undef, 0)
        qvarA  = Vector{Function}(undef, 0)

        return new(Xset,Kset,Nset,
                   fconst,fbaseA,fvarA,fbaseB,fvarB,
                   qconst,qbaseA,qvarA,qbaseB,qvarB)
    end
end

# A structure containing solution data:
mutable struct Solution
    t::Vector{Float64}  # time vector
    q::Matrix{Float64}  # displacement matrix
    v::Matrix{Float64}  # velocity matrix
    a::Matrix{Float64}  # acceleration matrix
    ls::Matrix{Float64} # last step data

    # constructing default (empty) solution data:
   # function Solution(;t = Vector{Float64}(undef, 0),
   #                    q = Matrix{Float64}(undef, 0, 0),
   #                    v = Matrix{Float64}(undef, 0, 0),
   #                    a = Matrix{Float64}(undef, 0, 0),
   #                    ls = Matrix{Float64}(undef, 0, 0))
   #     return new(t,q,v,a,ls)
   # end
end

# A structure containing elemental (unaveraged) tensor field data:
mutable struct ElementNodeData
    data::Array{Float64}  # 3D matrix of unaveraged tensor field
    eTags::Vector{Int32}    # vector of element tags
    nonpe::Vector{Int32}    # vector of number of nodes per elements
    name::Symbol

     # constructing default (empty) unaveraged tensor field:
     function ElementNodeData()
        data  = Array{Float64, 3}(undef, 0, 0, 0)
        eTags = Vector{Int32}(undef, 0)
        nonpe = Vector{Int32}(undef, 0)
        name  = :general
        return new(data, eTags, nonpe, name)
     end
end

# A structure containing nodal (averaged) tensor field data:
mutable struct NodeData
    data::Array{Float64}  # 3D matrix of unaveraged tensor field
    name::Symbol

     # constructing default (empty) unaveraged tensor field:
     function NodeData()
        data = Array{Float64, 3}(undef, 0, 0, 0)
        name = :general
        return new(data, name)
     end
end

# A structure containing nodal (averaged) tensor field data:
mutable struct ElementData
    data::Array{Float64}  # 3D matrix of unaveraged tensor field
    eTags::Vector{Int32}    # vector of element tags
    name::Symbol
    
     # constructing default (empty) unaveraged tensor field:
     function ElementData()
        data = Array{Float64, 3}(undef, 0, 0, 0)
        eTags = Vector{Int32}(undef, 0)
        name = :general
        return new(data, eTags, name)
     end
end

# A structure containing results data:
mutable struct Result
    Sn::NodeData
    Se::ElementNodeData
    En::NodeData
    Ee::ElementNodeData
    Un::NodeData
    Ue::ElementNodeData
end


#==========================================================================================================================
                                            II. LOW LEVEL FUNCTIONS
===========================================================================================================================#

# A function for selecting nodes
function nodeSelelection(mesh,seldef=:all,params=[],presel=:none,name="nset1",id=0,tol=0.001)

    dim = 0;
    elems = Dict{Int16, Vector{Int32}}()

    # get coordinates to search from:
    if presel == :none
        xyz = mesh.xyz
        nTags = collect(1:size(xyz,1))
    else
        presel = scope2selection(mesh,presel) # scoped geometry set
        nTags = presel.phnTags
        xyz = mesh.xyz[nTags,:]
    end

    if params isa Number
        params = to_vector(params)
    end

    if seldef == :all
        nset = nTags
    elseif seldef ∈ (:point, :p)
        if params isa Vector
            params = reshape(params, 1, :)
        end
        if size(params,2) == 2
            params = [params zeros(size(params,1),1)]
        end
        nset = nTags[dsearchn(xyz,params)]
    elseif seldef ∈ (:x, :y, :z)
        col = seldef == :x ? 1 : seldef == :y ? 2 : 3
        vals = xyz[:, col]
        if size(params,2) == 1
            nset = []
            for (i, row) in enumerate(eachrow(params))
                append!(nset,nTags[vals .≈ row[1]])
            end
        elseif length(params) == 2
            nset = []
            for (i, row) in enumerate(eachrow(params))
                append!(nset,nTags[(vals .>= row[1]-tol) .& (vals .<= row[2]+tol)])
            end
        end  
    elseif seldef ∈ (:line, :l)
        # not yet implemented              
    elseif seldef ∈ (:box, :b, :boxboundary, :boxb, :bb, :boxinverse, :boxi, :bi)    
        cmp = seldef ∈ (:box, :b)                   ? (<) :
                seldef ∈ (:boxboundary, :boxb, :bb) ? (≈) :
                seldef ∈ (:boxinverse, :boxi, :bi)  ? (>) :
                error("Invalid seldef: $seldef")
        nset = []
        for (i, row) in enumerate(eachrow(params))
            append!(nset,nTags[cmp.(dbox(xyz, row), tol)])
        end
    end

    return Selection(name,dim,id,nset,elems)
end

# A function for selecting elements
function elementSelelection(mesh,obj,seldef=:all,params=[],presel=:none,name="eset1",id=0,tol=0.001)

    elems = Dict{Int16, Vector{Int32}}()  # initialization of empty elems for Selection
    nset = Int32[]                        # initialization of empty nodes for Selection

    dim = obj == :p ? 0 : obj == :l ? 1 : obj ∈ (:s, :a, :f) ? 2 : obj == :v ? 3 : error("Invalid object: $object")

    if params isa Number
        params = to_vector(params)
    end

    for (eType, eData) in mesh.elem  # for all element types of the mesh
        if eData.dim == dim    # if the dimension of the current element type corresponds to the required
            xyz = eData.CP   # center point coordinates of the current element type
            eTags = eData.tags # element tags of the current element type
            nnET  = eData.nn
            if presel != :none
                presel   = scope2selection(mesh,presel) # scoped geometry set
                if haskey(presel.elems,eType)
                    pseTags = spresel.elems[eType]
                    psid = in.(eTags, Ref(pseTags))
                    xyz = xyz[psid,:]
                    eTags = eTags[psid]
                    nnET = nnET[psid,:]
                end
            end
            if seldef == :all
                selid = trues(size(eTags))
            elseif seldef ∈ (:point, :p)
                if params isa Vector
                    params = reshape(params, 1, :)
                end
                if size(params,2) == 2
                    params = [params zeros(size(params,1),1)]
                end
                selid = dsearchn(xyz,params)
            elseif seldef ∈ (:x, :y, :z)
                col = seldef == :x ? 1 : seldef == :y ? 2 : 3
                vals = xyz[:, col]
                n = size(xyz,1);
                selid = falses(n)
                if size(params,2) == 1
                    for (i, row) in enumerate(eachrow(params))
                        selid = selid .| (vals .≈ row[1])
                    end
                elseif size(params,2) == 2
                    for (i, row) in enumerate(eachrow(params))
                        selid = selid .| (vals .>= row[1]-tol) .& (vals .<= row[2]+tol)
                    end
                end  
            elseif seldef ∈ (:line, :l)
                # not yet implemented
            elseif seldef ∈ (:box, :b, :boxboundary, :boxb, :bb, :boxinverse, :boxi, :bi)    
                cmp = seldef ∈ (:box, :b)                 ? (<) :
                      seldef ∈ (:boxboundary, :boxb, :bb) ? (≈) :
                      seldef ∈ (:boxinverse, :boxi, :bi)  ? (>) :
                      error("Invalid seldef: $seldef")
                n = size(xyz,1);
                selid = falses(n)
                for (i, row) in enumerate(eachrow(params))
                    selid = selid .| cmp.(dbox(xyz, params), tol)
                end

            end
            elems[eType] = eTags[selid]
            append!(nset,nnET[selid,:][:])
        end
    end
    nset = unique(nset)

    return Selection(name,dim,id,nset,elems)
end

function sel(obj,seldef=:all,selpar=[];presel=:none,name="set1",id=0,tol=0.001)
    Selector(obj,seldef,selpar,presel,id,name,tol)
end

function nsel(seldef=:all,selpar=[];presel=:none,name="nset1",id=0,tol=0.001)
    Selector(:n,seldef,selpar,presel,id,name,tol)
end

function psel(seldef=:all,selpar=[];presel=:none,name="pset1",id=0,tol=0.001)
    Selector(:p,seldef,selpar,presel,id,name,tol)
end

function lsel(seldef=:all,selpar=[];presel=:none,name="lset1",id=0,tol=0.001)
    Selector(:l,seldef,selpar,presel,id,name,tol)
end

function asel(seldef=:all,selpar=[];presel=:none,name="aset1",id=0,tol=0.001)
    Selector(:a,seldef,selpar,presel,id,name,tol)
end

function vsel(seldef=:all,selpar=[];presel=:none,name="vset1",id=0,tol=0.001)
    Selector(:v,seldef,selpar,presel,id,name,tol)
end

# A function creating degree of freedoam table
function dofMatrix(mesh::Mesh,sec)
    n = mesh.non     # number of nodes
    s = length(sec)  # number of sections

    dimlist = zeros(Int32,s)  # list of dimensions associated with sections
    pdimlist = zeros(Int32,s)  # list of physical dimensions associated with sections
    i = 1

    for (secKey, secData) in pairs(sec)  # for each section
        dim1 = secData.dim  # spatial dimension of the current section
        #pdim1 = secData.pdim  # physical dimension of the current section
        dimlist[i] = dim1   # adding current dimension to dimension list
        #pdimlist[i] = pdim1  # adding current physical dimension to physical dimension list
        i += 1
    end

    # checking dimension consistency and establishing global dim and pdim values
    if all(x -> x == dimlist[1], dimlist)
        dim = dimlist[1]          # spatial dimesnion of the model
        # pdim = maximum(pdimlist)  # physical dimesion of the model
    else
        error("Dimension inconsistency!")
    end

    dim == 2 ? pdim = 3 : pdim = 6

    dofTable01 = falses(pdim,n)  # initialization of boolean dof table

    for (secKey, secData) in pairs(sec)  # for each section
        pdim1 = secData.pdim         # physical dimension of the current section
        scope = secData.scope       # name of the current section
        setData = scope2selection(mesh,scope)  # convert scope to selection structure
      #  setData = set(mesh,secName)  # retrieve current set data from current section name              
        nTags  = setData.phnTags       # node tags of the current set
        dofTable01[1:pdim1,nTags] .= true  # activate corresponding dofs
    end

    true_indices = findall(dofTable01)                # vector of true indices
    dofTable = zeros(Int32, size(dofTable01))         # create a matrix of zeros with same size as dofTable01
    dofTable[true_indices] = 1:length(true_indices)   # fill in the sequential numbers at the true positions
    dofTable = dofTable'

    return Matrix(dofTable)
end

# A function merging system matrices:
function mergeSystemMatrices(sm_list::Union{Vector{SystemMatrices}, Dict{Any, SystemMatrices}})
    sm_values = values(sm_list) isa Vector ? sm_list : collect(values(sm_list))
    
    merged = SystemMatrices()

    if all(isempty, [getfield(bc, :Ke) for bc in sm_values if hasfield(typeof(bc), :Ke)])
        merged.Ke = Matrix{Float64}(undef, 0, 0)
    else
        merged.Ke = sum([bc.Ke for bc in sm_values if hasfield(typeof(bc), :Ke) && !isempty(bc.Ke)])
    end

    if all(isempty, [getfield(bc, :Ce) for bc in sm_values if hasfield(typeof(bc), :Ce)])
        merged.Ce = Matrix{Float64}(undef, 0, 0)
    else
        merged.Ce = sum([bc.Ce for bc in sm_values if hasfield(typeof(bc), :Ce) && !isempty(bc.Ce)])
    end

    if all(isempty, [getfield(bc, :Ks) for bc in sm_values if hasfield(typeof(bc), :Ks)])
        merged.Ks = Matrix{Float64}(undef, 0, 0)
    else
        merged.Ks = sum([bc.Ks for bc in sm_values if hasfield(typeof(bc), :Ks) && !isempty(bc.Ks)])
    end

    if all(isempty, [getfield(bc, :Cs) for bc in sm_values if hasfield(typeof(bc), :Cs)])
        merged.Cs = Matrix{Float64}(undef, 0, 0)
    else
        merged.Cs = sum([bc.Cs for bc in sm_values if hasfield(typeof(bc), :Cs) && !isempty(bc.Cs)])
    end

    merged.G = isempty([bc.G for bc in sm_values if hasfield(typeof(bc), :G) && !isempty(bc.G)]) ? 
            Matrix{Float64}(undef, 0, 0) : 
            vcat([bc.G for bc in sm_values if hasfield(typeof(bc), :G) && !isempty(bc.G)]...)

    return merged
end

# A function merging boundary condition definitions:
function mergeBoundaryConditions(bc_list::Union{Vector{BoundaryConditions}, Dict{Any, BoundaryConditions}})
    bc_values = values(bc_list) isa Vector ? bc_list : collect(values(bc_list))
    
    merged = BoundaryConditions()

    merged.Xset = unique(vcat([bc.Xset for bc in bc_values if !isempty(bc.Xset)]...))
    merged.Kset = unique(vcat([bc.Kset for bc in bc_values if !isempty(bc.Kset)]...))
    merged.Nset = unique(vcat([bc.Nset for bc in bc_values if !isempty(bc.Nset)]...))
    
    merged.fconst = sum([bc.fconst for bc in bc_values if !isempty(bc.fconst)])

   #=qconst1 = isempty([bc.qconst for bc in bc_values if !isempty(bc.qconst)]) ? Matrix{Float64}(undef, 0, 0) : hcat([bc.qconst for bc in bc_values if !isempty(bc.qconst)]...)
    if !isempty(qconst1)
        merged.qconst = sum(qconst1, dims=2)
    else
        merged.qconst = []
    end =#
    
   #= for i in eachindex(bc_values)
        qconst1 = bc_values[i].qconst
        println("qconst1= ",qconst1)
        if !isempty(qconst1)
            addit = true
            break
        else
            addit = false
        end
        println("addit1= ",addit)
    end
    println("addit2= ",addit)
    if addit
        merged.qconst = sum([bc.qconst for bc in bc_values if !isempty(bc.qconst)])
    else
        merged.qconst = [];
    end =#

    merged.qconst = []
    for i in eachindex(bc_values)
        qconst1 = bc_values[i].qconst
        if !isempty(qconst1) && isempty(merged.qconst)
            merged.qconst = qconst1
        elseif !isempty(qconst1) && !isempty(merged.qconst)
            merged.qconst += qconst1
        end
    end

    # merged.qconst = reduce(+, [bc.qconst for bc in bc_values if !isempty(bc.qconst)]; init=Float64[])

    merged.fbaseA = isempty([bc.fbaseA for bc in bc_values if !isempty(bc.fbaseA)]) ? Matrix{Float64}(undef, 0, 0) : hcat([bc.fbaseA for bc in bc_values if !isempty(bc.fbaseA)]...)
    merged.qbaseA = isempty([bc.qbaseA for bc in bc_values if !isempty(bc.qbaseA)]) ? Matrix{Float64}(undef, 0, 0) : hcat([bc.qbaseA for bc in bc_values if !isempty(bc.qbaseA)]...)
    merged.fbaseB = isempty([bc.fbaseB for bc in bc_values if !isempty(bc.fbaseB)]) ? Matrix{Float64}(undef, 0, 0) : hcat([bc.fbaseB for bc in bc_values if !isempty(bc.fbaseB)]...)
    merged.qbaseB = isempty([bc.qbaseB for bc in bc_values if !isempty(bc.qbaseB)]) ? Matrix{Float64}(undef, 0, 0) : hcat([bc.qbaseB for bc in bc_values if !isempty(bc.qbaseB)]...)

    merged.fvarA = isempty([bc.fvarA for bc in bc_values if !isempty(bc.fvarA)]) ? Vector{Function}(undef, 0) : vcat([bc.fvarA for bc in bc_values if !isempty(bc.fvarA)]...)
    merged.qvarA = isempty([bc.qvarA for bc in bc_values if !isempty(bc.qvarA)]) ? Vector{Function}(undef, 0) : vcat([bc.qvarA for bc in bc_values if !isempty(bc.qvarA)]...)
    merged.fvarB = isempty([bc.fvarB for bc in bc_values if !isempty(bc.fvarB)]) ? Matrix{Float64}(undef, 0, 0) : vcat([bc.fvarB for bc in bc_values if !isempty(bc.fvarB)]...)
    merged.qvarB = isempty([bc.qvarB for bc in bc_values if !isempty(bc.qvarB)]) ? Matrix{Float64}(undef, 0, 0) : vcat([bc.qvarB for bc in bc_values if !isempty(bc.qvarB)]...)

    return merged
end

# A function returning multi point constraint matrix:
function mpcMatrix(mesh, mpc, dofTable)

    sys  = SystemMatrices() # initialization of an empty boundary condition structure

    mpc = to_vector(mpc)        # if mpc is not a vector/dictionary, convert to a vector
    xyz = mesh.xyz              # coordinate matrix

    dof  = maximum(dofTable)    # number of degrees of freedom of the model
    pdim = size(dofTable,2)     # physical dimension of the model

    ridx = []  # initialization of G matrix row index vector
    cidx = []  # initialization of G matrix column index vector
    data = []  # initialization of G matrix data vector
    data = convert(Vector{Float64}, data)
    row = 0    # number of rows of G matrix

    for (mpcKey, mpcData) in pairs(mpc)  # for each mpc set
        scope1vec = mpcData.scope1
        scope2vec = mpcData.scope2
        pairdefvec = mpcData.pairdef
        contypevec = mpcData.contype
        showvec = mpcData.show

        n = length(scope1vec)           # number of mpcs in the current mpc set

        for i = 1:n
            scope1 = scope1vec[i]
            scope2 = scope2vec[i]
            pairdef = pairdefvec[i]
            contype = contypevec[i]
            show    = showvec[i]

            if !isa(scope1,Selection)
                scope1 = set(mesh,scope1)  # retrieve set data from scope1
            end
            nTags1 = scope1.phnTags
            nCrds1 = xyz[nTags1,:]
            numNodes1 = length(nTags1)
            
            if !isa(scope2,Selection)
                scope2 = set(mesh,scope2)  # retrieve set data from scope2
            end
            nTags2 = scope2.phnTags
            nCrds2 = xyz[nTags2,:]
            numNodes2 = length(nTags2)

            # Setting default pairing definitions:
            if pairdef == :default
                if numNodes1 == 1
                    pairdef = :masterslave   
                    mCrd = nCrds1
                    sCrd = nCrds2
                    mSet = nTags1
                    sSet = nTags2
                elseif numNodes2 == 1    
                    pairdef = :masterslave   
                    mCrd = nCrds2
                    sCrd = nCrds1
                    mSet = nTags2
                    sSet = nTags1                
                end
            end

            if show
                mpcDepicter(mSet,sSet,mCrd,sCrd)
            end

            if pdim == 3
                mux = dofTable[mSet,1]
                muy = dofTable[mSet,2]
                mrz = dofTable[mSet,3]     
                
                nummpc = length(sSet) 

                for j = 1:nummpc  # for each mpc
                    nTag2 = sSet[j]
                    # Extract coordinates
                    x12 = sCrd[j, 1] - mCrd[1, 1]
                    y12 = sCrd[j, 2] - mCrd[1, 2] 

                    sux = dofTable[nTag2,1]
                    suy = dofTable[nTag2,2] 
                    
                    append!(ridx, row .+ [1;1;1;2;2;2])
                    append!(cidx, [mux;mrz;sux;muy;mrz;suy])
                    append!(data, [-1;y12;1;-1;-x12;1])        
                    
                    row+=2                 
                end

            elseif pdim == 6
                mux = dofTable[mSet,1]
                muy = dofTable[mSet,2]
                muz = dofTable[mSet,3]
                mrx = dofTable[mSet,4]
                mry = dofTable[mSet,5]
                mrz = dofTable[mSet,6]

                nummpc = length(sSet) 

                for j = 1:nummpc  # for each mpc

                    nTag2 = sSet[j]

                    # Extract coordinates
                    x12 = sCrd[j, 1] - mCrd[1, 1]
                    y12 = sCrd[j, 2] - mCrd[1, 2]
                    z12 = sCrd[j, 3] - mCrd[1, 3]

                    sux = dofTable[nTag2,1]
                    suy = dofTable[nTag2,2]
                    suz = dofTable[nTag2,3]

                    append!(ridx, row .+ [1;1;1;1;2;2;2;2;3;3;3;3])
                    append!(cidx, [mux;mry;mrz;sux;muy;mrx;mrz;suy;muz;mrx;mry;suz])
                    append!(data, [-1;-z12;y12;1;-1;z12;-x12;1;-1;-y12;x12;1])

                    row += 3
                end
            end
        end
    end

    sys.G = sparse(ridx, cidx, data, row, dof)
    dropzeros!(sys.G)

    return sys
end

# A function returning spring matrix:
function springMatrix(mesh, spring, dofTable)

    sys  = SystemMatrices() # initialization of an empty boundary condition structure

    spring = to_vector(spring)  # if spring is not a vector/dictionary, convert to a vector
    xyz = mesh.xyz              # coordinate matrix
    
    dof  = maximum(dofTable)    # number of degrees of freedom of the model
    epdim = size(dofTable,2)    # extended physical dimension of the model   

    Ik = []  # initialization of stiffness row index vector
    Jk = []  # initialization of stiffness column index vector
    Vk = []  # initialization of stiffness value vector
    Vk = convert(Vector{Float64}, Vk)
    Ic = []  # initialization of damping row index vector
    Jc = []  # initialization of damping column index vector    
    Vc = []  # initialization of damping value vector
    Vc = convert(Vector{Float64}, Vc)
    
    for (springKey, springData) in pairs(spring)  # for each spring set
        scope1vec      = springData.scope1
        scope2vec      = springData.scope2
        kvec           = springData.k
        cvec           = springData.c
        pairdefvec     = springData.pairdef
        orientationvec = springData.orientation

        n = length(scope1vec)           # number of springs in the current mpc set
            
        for i = 1:n
            scope1 = scope1vec[i]
            scope2 = scope2vec[i]
            k      = kvec[i]
            c      = cvec[i]
            pairdef = pairdefvec[i]
            orientation = orientationvec[i]

            pdimS = max(length(k),length(c)) # physical dimension of spring element

            if epdim == 3 
                if pdimS == 1 || pdimS == 2
                    pdim = 2
                else
                    pdim = 3
                end
            elseif epdim == 6
                if pdimS <= 3
                    pdim = 3
                else
                    pdim = 6
                end
            end

            # Extracting nodes and coordinates from the 1st set:
            if scope1 == :none
                nTags1 = mesh.nTags
                nCrds1 = xyz
            else
                if !isa(scope1,Selection)
                    scope1 = set(mesh,scope1)  # retrieve set data from scope1
                end
                nTags1 = scope1.phnTags
                nCrds1 = xyz[nTags1,:]
            end
            numNodes1 = length(nTags1)   

            # Extracting nodes and coordinates from the 2nd set:

            if scope2 == :fixed   # spring is connected to the ground -> no nodes from here
                fixed = true
                nTags2 = []
                nCrds2 = []
                numNodes2 = []
                Iidx = zeros(Int, pdim, pdim)
                Jidx = zeros(Int, pdim, pdim)
                for i in 1:pdim, j in 1:pdim
                    Iidx[i, j] = j
                    Jidx[i, j] = i
                end
            else
                fixed = false
                if scope2 == :none   # 2nd set is empty -> no nodes from here, calculate from 1st set
                    nTags2 = []
                    nCrds2 = []
                    numNodes2 = []  
                else 
                    if !isa(scope2,Selection)
                        scope2 = set(mesh,scope2)  # retrieve set data from scope1
                    end
                    nTags2 = scope2.phnTags
                    nCrds2 = xyz[nTags2,:]
                    numNodes2 = length(nTags2)    
                end
                Iidx = zeros(Int, 2 * pdim, 2 * pdim)
                Jidx = zeros(Int, 2 * pdim, 2 * pdim)
                for i in 1:2*pdim, j in 1:2*pdim
                    Iidx[i, j] = j
                    Jidx[i, j] = i
                end
            end

            # Setting default pairing definitions:
            if pairdef == :default || pairdef == :masterslave
                if numNodes1 == numNodes2
                    pairdef = :sequential
                elseif numNodes1 == 1
                    pairdef = :masterslave
                    masterCrd = nCrds1
                elseif numNodes2 == 1
                    pairdef = :masterslave
                    masterCrd = nCrds2
                else
                    pairdef = :coincident
                end
            end 
            
            # Setting default orientation definitions:
            if orientation == :default
                if pairdef == :masterslave
                    orientation = masterCrd
                elseif pairdef isa Tuple
                    orientation = :nodalax
                else
                    orientation = :gcsys
                end
            end   
            
            # Base stiffness and damping matrix in gcsys:
            if !isempty(k)
                K0 = Diagonal(k)
                if !fixed
                    K0 = [ K0 -K0;
                          -K0  K0]   
                end
            else
                K0 = []
            end

            if !isempty(c)
                C0 = Diagonal(c)
                if !fixed
                    C0 = [ C0 -C0;
                          -C0  C0]  
                end 
            else
                C0 = []
            end

            # Create node sets based on pairing definition:
            nSet1 = []
            nSet2 = []       
            
            pairDict = Dict("samex" => 1, "samey" => 2, "samez" => 3,
                            "samexy" => [1,2], "samexz" => [1,3], "sameyz" => [2,3],
                            "samexyz" => [1,2,3], "coincident" => [1,2,3]) 
                            
            if scope2 == :none
                if pairdef == :coincident
                    for (i, row1) in enumerate(eachrow(nCrds1[1:end-1,:]))
                        for (j, row2) in enumerate(eachrow(nCrds1[i+1:end,:]))
                            if row1 ≈ row2
                                push!(nSet1, nTags1[i])
                                push!(nSet2, nTags1[j])
                            end
                        end
                    end                     
                elseif pairdef isa Tuple 
                    if pairdef[1] == :lessthan
                        tol = pairdef[2]
                        for (i, row1) in enumerate(eachrow(nCrds1[1:end-1,:]))
                            for (j, row2) in enumerate(eachrow(nCrds1[i+1:end,:]))
                                distance = norm(row1 .- row2)
                                if distance < tol
                                    push!(nSet1, nTags1[i])
                                    push!(nSet2, nTags1[j])
                                end
                            end
                        end
                    elseif pairdef[1] == :equals
                        dist = paird[2]
                        for (i, row1) in enumerate(eachrow(nCrds1[1:end-1,:]))
                            for (j, row2) in enumerate(eachrow(nCrds1[i+1:end,:]))
                                distance = norm(row1 .- row2)
                                if distance ≈ dist
                                    push!(nSet1, nTags1[i])
                                    push!(nSet2, nTags1[j])
                                end
                            end
                        end
                    end
                elseif haskey(pairdef,pairDict)
                    d = pairDict(paird)
                    for (i, row1) in enumerate(eachrow(nCrds1[1:end-1,:]))
                        for (j, row2) in enumerate(eachrow(nCrds1[i+1:end,:]))
                            if row1[d] ≈ row2[d]
                                push!(nSet1, nTags1[i])
                                push!(nSet2, nTags1[j])
                            end
                        end
                    end 
                end                
            elseif !fixed
                if pairdef == :masterslave
                    if numNodes1 == 1
                        nSet1 = fill(nTags1[1],numNodes2)
                        nSet2 = nTags2
                    elseif numNodes2 == 1
                        nSet1 = nTags1
                        nSet2 = fill(nTags2[1],numNodes1)
                    end
                elseif pairdef == :sequential
                    nSet1 = nTags1
                    nSet2 = nTags2
                elseif pairdef == :coincident
                    for (i, row1) in enumerate(eachrow(nCrds1))
                        for (j, row2) in enumerate(eachrow(nCrds2))
                            if row1 ≈ row2
                                push!(nSet1, nTags1[i])
                                push!(nSet2, nTags2[j])
                            end
                        end
                    end                     
                elseif pairdef isa Tuple 
                    if pairdef[1] == :lessthan
                        tol = pairdef[2]
                        for (i, row1) in enumerate(eachrow(nCrds1))
                            for (j, row2) in enumerate(eachrow(nCrds2))
                                distance = norm(row1 .- row2)
                                if distance < tol
                                    push!(nSet1, nTags1[i])
                                    push!(nSet2, nTags2[j])
                                end
                            end
                        end
                    elseif pairdef[1] == :equals
                        dist = pairdef[2]
                        for (i, row1) in enumerate(eachrow(nCrds1))
                            for (j, row2) in enumerate(eachrow(nCrds2))
                                distance = norm(row1 .- row2)
                                if distance ≈ dist
                                    push!(nSet1, nTags1[i])
                                    push!(nSet2, nTags2[j])
                                end
                            end
                        end
                    end
                elseif haskey(pairdef,pairDict)
                    d = pairDict(pairdef)
                    for (i, row1) in enumerate(eachrow(nCrds1))
                        for (j, row2) in enumerate(eachrow(nCrds2))
                            if row1[d] ≈ row2[d]
                                push!(nSet1, nTags1[i])
                                push!(nSet2, nTags2[j])
                            end
                        end
                    end 
                end
            else
                nSet1 = nTags1
            end  
            
            # Calculate coordinates based on orientation:

            if orientation isa Array  # masterpoint case
                mCrd = orientation
                # Crds1 = repeat(mCrd,numNodes1,1)
                if numNodes1 == 1
                    Crds1 = repeat(mCrd,numNodes2,1)
                    Crds2 = nCrds2
                else
                    Crds1 = repeat(mCrd,numNodes1,1)
                    Crds2 = nCrds1
                end
            elseif orientation == :masterCrd
                if numNodes1 == 1
                    Crds1 = repeat(nCrds1,numNodes2,1)
                    Crds2 = nCrds2
                elseif numNodes2 == 1
                    Crds1 = nCrds1
                    Crds2 = repeat(nCrds2,numNodes1,1)       
                end                
            elseif orientation == :nodalax
                Crds1 = nCrds1[nSet1,:]
                Crds2 = nCrds2[nSet2,:]
            else
                Crds1 = nCrds1
                Crds2 = nCrds2
            end 
            
            numSprings = length(nSet1)         # number of Springs

            if orientation isa Array || orientation == :nodalax
                for j = 1:numSprings # for each spring element:
                    # Extract coordinates
                    Crd1 = SVector{3, Float64}(Crds1[j, :])
                    Crd2 = SVector{3, Float64}(Crds2[j, :])

                    # Calculation of the normal unit direction vector:
                    en = normalize(Crd2 - Crd1)

                    if pdimS >= 2
                        # Calculation of the tangential unit direction vector:
                        ey = SVector(0.0, 1.0, 0.0)     # unit vector in y direction
                        if en ≈ ey     # slave node is on y axis
                            et = SVector(-1.0, 0.0, 0.0)
                        elseif en ≈ -ey
                            et = SVector(1.0, 0.0, 0.0)
                        else
                            dp = dot(ey, en)   # dot product
                            if dp < 1e-10 
                                et = SVector(0.0, 1.0, 0.0)
                            else
                                et = normalize(ey - dp * en) # tangential direction
                            end
                        end                                            
                    else
                        et = nothing
                    end

                    if pdimS == 3
                        # Calculation of the binormal unit direction vector:
                        eb = normalize(cross(en, et))
                    else
                        eb = nothing
                    end

                    T = hcat(filter(x -> x !== nothing, [en, et, eb])...)
                    
                    if pdim == 6
                        T = [T zeros(3,3)
                             zeros(3,3) T]
                    else
                        T = T[1:pdim,:]  
                    end

                    if !fixed
                        T = [T zeros(pdim,pdimS);
                            zeros(pdim,pdimS) T]
                    end

                    nTag1 = nSet1[j]

                    if fixed
                        nn2 = dofTable[nTag1,1:pdim]
                    else
                        nTag2 = nSet2[j]
                        nn2 = [dofTable[nTag1,1:pdim]; dofTable[nTag2,1:pdim]]
                    end

                    if !isempty(k)
                        append!(Ik, nn2[Iidx[:]])
                        append!(Jk, nn2[Jidx[:]])
                        K1 = T * K0 * T'
                        append!(Vk, K1[:])
                    end

                    if !isempty(c)
                        append!(Ic, nn2[Iidx[:]])
                        append!(Jc, nn2[Jidx[:]])
                        C1 = T * C0 * T'
                        append!(Vc, C1[:])
                    end
                end
            else
                if orientation == :normalto
                    p1 = Crds1[1, :]
                    p2 = Crds1[2, :]
                    p3 = Crds1[3, :]             
                    v1 = SVector(p2 - p1)
                    v2 = SVector(p3 - p1)
                    cp = cross(v1, v2)
                    if norm(cp) < 1e-10
                        j = 4
                        while norm(cp) <= 1e-10
                            p3 = Crds1[j, :] 
                            v2 = SVector(p3 - p1)
                            cp = cross(v1, v2)
                            j += 1
                        end
                    end
                    en = normalize(cp)

                    if pdimS >= 2
                        # Calculation of the tangential unit direction vector:
                        ey = SVector(0.0, 1.0, 0.0)     # unit vector in y direction
                        if en ≈ ey     # slave node is on y axis
                            et = SVector(-1.0, 0.0, 0.0)
                        elseif en ≈ -ey
                            et = SVector(1.0, 0.0, 0.0)
                        else
                            dp = dot(ey, en)   # dot product
                            if dp < 1e-10 
                                et = SVector(0.0, 1.0, 0.0)
                            else
                                et = normalize(ey - dp * en) # tangential direction
                            end
                        end
                    else
                        et = nothing
                    end

                    if pdimS == 3
                        # Calculation of the binormal unit direction vector:
                        eb = normalize(cross(en, et))
                    else
                        eb = nothing
                    end                

                    T = hcat(filter(x -> x !== nothing, [en, et, eb])...)[1:pdim,:]

                    if !fixed
                        T = [T zeros(pdim,pdimS);
                            zeros(pdim,pdimS) T]
                    end

                    if !isempty(k)
                        K1 = T * K0 * T'
                    end
                    if !isempty(c)
                        C1 = T * C0 * T'
                    end
                elseif orientation == :gcsys
                    K1 = K0
                    C1 = C0
                end

                for j = 1:numSprings
                    nTag1 = nSet1[j]

                    if fixed
                        nn2 = dofTable[nTag1,1:pdim]
                    else
                        nTag2 = nSet2[j]
                        nn2 = [dofTable[nTag1,1:pdim]; dofTable[nTag2,1:pdim]]
                    end

                    if !isempty(k)
                        append!(Ik, nn2[Iidx[:]])
                        append!(Jk, nn2[Jidx[:]])
                        append!(Vk, K1[:])
                    end
                    if !isempty(c)
                        append!(Ic, nn2[Iidx[:]])
                        append!(Jc, nn2[Jidx[:]])                   
                        append!(Vc, C1[:])
                    end
                end 
            end
        end     
    end

    if !isempty(Ik)
        sys.Ks = sparse(Ik, Jk, Vk, dof, dof)
        dropzeros!(sys.Ks)
    end

    if !isempty(Ic)   
        sys.Cs = sparse(Ic, Jc, Vc, dof, dof)
        dropzeros!(sys.Cs)
    end

    return sys
end

# A function returning kinematic boundary conditions ---------------------------------------------------------------------
function constraint(mesh::Mesh,K::SparseMatrixCSC{Float64,Int},disp,dofTable::Matrix{Int32})

    bcs  = BoundaryConditions() # initialization of an empty boundary condition structure
    dof  = maximum(dofTable)    # number of degrees of freedom of the model
  #  bcs.qconst = zeros(dof)     # initialization of constant displacement vector
    pdim = size(dofTable,2)

    disp = to_vector(disp)      # if displacement is not a vector/dictionary, convert to a vector

    for (dispKey, dispData) in pairs(disp)  # for each displacement set
        scopevec  = dispData.scope   # vector of scoped geometry identifiers in the current displacement set
        predefvec = dispData.predef  # vector of predefined displacement definitions in the current displacement set
        dirvec    = dispData.dir     # vector of direction definitions in the current displacement set
        uxvec     = dispData.ux      # vector of displacements in x direction in the current displacement set
        uyvec     = dispData.uy      # vector of displacements in y direction in the current displacement set
        uzvec     = dispData.uz      # vector of displacements in z direction in the current displacement set
        rxvec     = dispData.rx      # vector of rotations about x axis in the current displacement set
        ryvec     = dispData.ry      # vector of rotations about y axis in the current displacement set
        rzvec     = dispData.rz      # vector of rotations about z axis in the current displacement set
        tdvec     = dispData.td      # vector of time dependency function in the current displacement set

        n = length(scopevec)         # number of defined displacements in the current displacement set

        for i = 1:n                  # for each displacement in the current displacement set
            scope   = scopevec[i]    # scoped geometry of the ith displacement
            predef  = predefvec[i]   # predefined displacement definition of the ith displacement
            dir     = dirvec[i]      # direction definition of the ith displacement
            ux      = uxvec[i]       # displacement in x direction
            uy      = uyvec[i]       # displacement in y direction
            uz      = uzvec[i]       # displacement in z direction
            rx      = rxvec[i]       # rotation about x axis
            ry      = ryvec[i]       # rotation about y axis
            rz      = rzvec[i]       # rotation about z axis
            td      = tdvec[i]       # time dependency function of the ith displacement

            scope = scope2selection(mesh,scope)  # convert scope to selection structure

           #= if scope isa Selector
                if scope.obj == :n
                    scope = nodeSelelection(mesh,scope.seldef,scope.pardef,scope.presel,scope.name,scope.id,scope.tol)
                else
                    scope =  elementSelelection(mesh,scope.obj,scope.seldef,scope.pardef,scope.presel,scope.name,scope.id,scope.tol)
                end
            elseif scope isa Tuple && scope[1] isa Symbol
                scope = nsel(mesh,scope[1],scope[2])
            elseif !isa(scope,Selection)
                scope = set(mesh,scope)       # (dim, phTag) pair of the entity on which displacement is applied (scoped geometry)       
            end =#
            nset = scope.phnTags          # node ids belonging to the scoped geometry
            non  = length(nset)           # number of nodes belonging to the scoped geometry

            if predef == :fixed              # if predefined displacement is fixed
                append!(bcs.Xset,filter(!=(0), vec(dofTable[nset, :])))  # adding corresponding dofs to Xset matrix
            else                             # if no predefined displacenet is defined
                if pdim == 3                 # if dimension is 2
                    U = [ux,uy,rz]           # general displacement vector of a node
                elseif pdim == 6             # if dimension is 3
                    U = [ux,uy,uz,rx,ry,rz]  # general displacement vector of a node
                end
               # Kset = []                 # initialization of Kset dof vector (dofs on which nonzero displacement is defined)
               # Kload = []                # initialization of displacement vector belonging to Kset
                Kset = Vector{Int32}(undef, 0) 
                Kload = Vector{Float64}(undef, 0) 
                for j in eachindex(U)     # for each possible displacement coordinate
                    if U[j] == 0                    # if the jth displacement coordinate is zero
                        append!(bcs.Xset, dofTable[nset,j])   # add the corresponding dof number to Xset dof vector
                    elseif U[j] isa Number          # if the jth displacement coordinate is a number (but it is not zero)
                        append!(Kset, dofTable[nset,j])       # add the corresponding dof number to Kset dof vector    
                        append!(Kload,fill(U[j],non))                 # add the given displacement to Kload vevtor
                    elseif U[j] isa Function        # if the jth displacement coordinate is a function
                        append!(Kset, dofTable[nset,j])       # add the corresponding dof number to Kset dof vector                     
                        xyz = mesh.xyz              # coordinate matrix
                        x   = xyz[nset,1]           # x coordinates
                        y   = xyz[nset,2]           # y coordinates
                        z   = xyz[nset,3]           # z coordinates
                        append!(Kload,U[j](x,y,z))  # add the numeric values (with substituted x,y,z coordinates) to Kload vector
                    end
                end
                append!(bcs.Kset,Kset)  # add the values of the jth Kset vector to the global Kset vector

                if td == :none && !isempty(Kset)  # if amplitude variation is not defined and Kset is not empty
                    if isempty(bcs.fconst)
                        bcs.fconst = -K[:,Kset]*Kload
                        bcs.qconst = zeros(dof)
                        bcs.qconst[Kset] = Kload
                    else
                        bcs.fconst -= K[:,Kset]*Kload  # add kinematic load to constant nodal load vector
                        bcs.qconst[Kset] = Kload
                    end
                 #   bcs.qconst[Kset] = Kload       # add displacement to constant nodal displacement vector
                elseif td isa Function && !isempty(Kset)  # if amplitude variation is a function and Kset is not empty
                    qbaseA1 = zeros(dof)           # initialization of base displacement vector for function for current disp def.  
                    qbaseA1[Kset] = Kload          # inserting displacements to base displacement vector for functions
                    if isempty(bcs.fbaseA)         # if global base load vector for functions is empty
                        bcs.fbaseA = -K[:,Kset]*Kload  # global base load matrix for functions = currnet base load vector for functions
                        bcs.qbaseA = qbaseA1           # global base displacement matrix for functions = current base displacement vector for functions    
                        bcs.qvarA  = td           # global displacement variation matrix for functions = current amplitude function        
                    else                           # if global base load vector for functions is not empty
                        bcs.fbaseA = [bcs.fbaseA -K[:,Kset]*Kload] # place current base load vector for functions next to current global base load matrix for functions
                        bcs.qbaseA = [bcs.qbaseA qbaseA1] # place current base displacement vector for functions next to current global base displacement matrix for functions
                        bcs.qvarA  = [bcs.qvarA;td]  # place current displacement variation function under current global displacement variation function matrix      
                    end
                elseif td isa Matrix && !isempty(Kset)  # if amplitude variation is a matrix and Kset is not empty
                    qbaseB1 = zeros(dof)     # initialization of base displacement vector for numeric variation for current disp def. 
                    qbaseB1[Kset] = Kload    # inserting displacements to base displacement vector for numeric variation
                    if isempty(bcs.fbaseB)   # if global base load vector for numeric variation is empty
                        bcs.fbaseB = -K[:,Kset]*Kload # global base load matrix for numeric variation = currnet base load vector for numeric variation
                        bcs.qbaseB = qbaseB1 # global base displacement matrix for numeric variation = current base displacement vector for numeric variation   
                        bcs.qvarB  = td     # global displacement variation matrix for numeric variation = current amplitude numeric variation          
                    else                     # if global base load vector for numeric variation is not empty
                        bcs.fbaseB = [bcs.fbaseB -K[:,Kset]*Kload] # place current base load vector for numeric variation next to current global base load matrix for functions
                        bcs.qbaseB = [bcs.qbaseB qbaseB1] # place current base displacement vector for numeric variation next to current global base displacement matrix for functions
                        bcs.qvarB  = [bcs.qvarB;td] # place current numeric displacement variation under current global numeric displacement variation matrix              
                    end
                end
            end
        end
    end
    return bcs   
end

# A function returning dynamic boundary conditions ---------------------------------------------------------------------
function loadVector(mesh,load,dofTable,type=:none)

    bcs = BoundaryConditions()  # initialization of an empty boundary condition structure

    xyz = mesh.xyz              # coordinate matrix
    xcrd   = xyz[:,1]           # x coordinates
    ycrd   = xyz[:,2]           # y coordinates
    zcrd   = xyz[:,3]           # z coordinates

    dof  = maximum(dofTable)    # number of degrees of freedom of the model
    bcs.fconst = zeros(dof)     # initialization of constant load vector

    load = to_vector(load)      # if load is not a vector/dictionary, convert to a vector

    for (loadKey, loadData) in pairs(load)  # for each load set
        scopevec = loadData.scope  # vector of scoped geometry identifiers in the current load set
        dirvec   = loadData.dir    # vector of direction definitions in the current load set
        typevec  = loadData.type   # vector of the types of the current load set
        fxvec    = loadData.fx     # vector of the loads in x direction in the current load set
        fyvec    = loadData.fy     # vector of the loads in y direction in the current load set
        fzvec    = loadData.fz     # vector of the loads in z direction in the current load set
        mxvec    = loadData.mx     # vector of the moments about x axis in the current load set
        myvec    = loadData.my     # vector of the moments about y axis in the current load set
        mzvec    = loadData.mz     # vector of the moments about z axis in the current load set
        tdvec    = loadData.td     # vector of time dependency function in the current load set
        widthvec = loadData.width  # vector of the widths in the current load set

        for (i,scope) in pairs(scopevec) # for each load in the current load set
            dir     = dirvec[i]      # direction definition of the ith load
            Ltype   = typevec[i]     # load type of the ith load
            fx      = fxvec[i]       # load in x direction
            fy      = fyvec[i]       # load in y direction
            fz      = fzvec[i]       # load in z direction
            mx      = mxvec[i]       # moment about x axis
            my      = myvec[i]       # moment about y axis
            mz      = mzvec[i]       # moment about z axis
            td      = tdvec[i]       # time dependency of the ith load
            width   = widthvec[i]    # width of the surface for the ith load

            scope = scope2selection(mesh,scope)  # convert scope to selection structure

           # if !isa(scope,Selection)
           #     scope   = set(mesh,scope) # scoped geometry set
           # end
            nset    = scope.phnTags     # node tags in the set
            elems   = scope.elems       # elemet type => tags Dictionary
            setDim  = scope.dim         # dimension of the scoped geometry

            non = length(nset)          # number of nodes of the ith load
            f = zeros(dof)              # initialization of load vector for the ith load

            # defining load type for default type:
            if Ltype == :default
                if non==1 || mx!=0 || my!=0 || mz!=0  # if the load is applied on 1 node or moment is defined
                    Ltype = :concentrated             # then the load type is defined as concentrated
                else
                    Ltype = :distributed              # in other cases the load type is defined as distributed
                end
            end

            if Ltype == :concentrated  # if the load type is concentrated
                pdim = size(dofTable,2)      # number of physical dofs

                if pdim == 3                 # if physical dimension is 3
                    FM = [fx,fy,mz]          # concnetrated loads for 2D model
                elseif pdim == 6             # if physical dimension is 3
                    FM = [fx,fy,fz,mx,my,mz] # concnetrated loads for 3D model
                end    
                
                for j in eachindex(FM)       # for each possible load component
                    if FM[j] != 0            # if the jth load component is nonzero
                        f[dofTable[nset,j]] .= FM[j]  # put it into the global nodal load vector
                    end
                end

            elseif Ltype == :distributed    # if the load type is distributed
                pdime = size(dofTable,2)    # number of extended physical dofs
                if pdime == 3
                    dim = 2
                    pdim = 2
                elseif pdime == 6
                    dim = 3
                    pdim = 3
                end

                # initialization of a load vector on one node
                if pdim == 3           # if physical dimension is 3 
                    f0 = [.0, .0, .0]
                elseif pdim == 2       # if physical dimension is 2
                    f0 = [.0, .0]
                elseif pdim == 1       # if physical dimension is 1
                    f0 = [.0]
                else
                    error("loadVector: dimension of the problem is $(problem.dim).")
                end
                
                for (eType, eTags) in elems # for each element type
                    nnET    = mesh.elem[eType].nn    # coonectivity matrix of current element type
                    eTagsET = mesh.elem[eType].tags  # element Tags of current element type
                    nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix of selected elements
                    nonpe   = mesh.elem[eType].non        # number of nodes per elements
                    noe     = length(eTags);              # number of selected elements

                    nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix
                    for i = 1:pdim
                        nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                    end
            
                    # properties of current element type:
                    # rows: shape functions, colums: Gauss points
                    HG  = mesh.elem[eType].HG[1]  # shape functions at Gauss points    
                    HξG = mesh.elem[eType].HξG[1] # ∂h/∂ξ derivatives at Gauss points
                    HηG = mesh.elem[eType].HηG[1] # ∂h/∂η derivatives at Gauss points
                    HζG = mesh.elem[eType].HζG[1] # ∂h/∂ζ derivatives at Gauss points

                    W   = mesh.elem[eType].GW[1]  # Gauss weigths      
                    gpe = length(W)            # number of Gauss points per elements

                    # elemental coordinate matrices:
                    # rows -> elements, colums -> element nodes
                    X = xcrd[nn]      # elemental x coordinate matrix
                    Y = ycrd[nn]      # elemental y coordinate matrix
                    Z = zcrd[nn]      # elemental z coordinate matrix

                    # Jacobi matrix elements arranged in matrices:
                    # rows -> elements, colums -> Gauss points
                    Xξ = X * HξG   # ∂x/∂ξ matrix
                    Yξ = Y * HξG   # ∂y/∂ξ matrix
                    Zξ = Z * HξG   # ∂z/∂ξ matrix
                    Xη = X * HηG   # ∂x/∂η matrix
                    Yη = Y * HηG   # ∂y/∂η matrix
                    Zη = Z * HηG   # ∂z/∂η matrix
                    Xζ = X * HζG   # ∂x/∂ζ matrix
                    Yζ = Y * HζG   # ∂y/∂ζ matrix
                    Zζ = Z * HζG   # ∂z/∂ζ matrix

                    # Coordinates at Gauss points arranged in matrices:
                    # rows -> elements, columns -> Gauss points
                    XG = X * HG   # x coordinates at Gauss ponts of the elements
                    YG = Y * HG   # y coordinates at Gauss ponts of the elements
                    ZG = Z * HG   # z coordinates at Gauss ponts of the elements

                    # initialization of transpose of approximation matrix:
                    AT = zeros(Float64, nonpe * pdim, pdim)
            
                    for e in 1:noe # for each selected element from current element type
                        f1 = zeros(Float64, pdim*nonpe) # elemental nodal load vector ini
                        for g in 1:gpe    # for each Gauss Point

                            # extracting Jacobi matrix elements:
                            xξ = Xξ[e,g]  # ∂x/∂ξ
                            yξ = Yξ[e,g]  # ∂y/∂ξ
                            zξ = Zξ[e,g]  # ∂z/∂ξ
                            xη = Xη[e,g]  # ∂x/∂η
                            yη = Yη[e,g]  # ∂y/∂η
                            zη = Zη[e,g]  # ∂z/∂η
                            xζ = Xζ[e,g]  # ∂x/∂ζ
                            yζ = Yζ[e,g]  # ∂y/∂ζ
                            zζ = Zζ[e,g]  # ∂z/∂ζ
                                        
                            # filling up aprroximation matrix:
                            for k ∈ 1:nonpe, l ∈ 1:pdim
                                AT[pdim*k-(pdim-l),l] = HG[k,g]
                            end   
                    
                            # extracting coordinates & filling up nodal load components
                            x = XG[e,g] # x coordinate at Gauss point g of element e
                            f0[1] = isa(fx, Function) ? fx(x, y, z) : fx
                            if pdim > 1
                                y = YG[e,g] # y coordinate at Gauss point g of element e
                                f0[2] = isa(fy, Function) ? fy(x, y, z) : fy
                                if pdim > 2
                                    z = ZG[e,g] # z coordinate at Gauss point g of element e
                                    f0[3] = isa(fz, Function) ? fz(x, y, z) : fz
                                end
                            end

                            # construction of Ja:
                            # volume load on 3D solid or shell:
                            if dim == 3 && (setDim == 3 || (setDim == 2 && type == :shell))
                                Ja = xξ * (yη * zζ - yζ * zη) -
                                     yξ * (xη * zζ - xζ * zη) +
                                     zξ * (xη * yζ - xζ * yη)
                            elseif dim == 3 && setDim == 2  # surface load on 3D solid or shell
                                xy = xξ*yη - yξ*xη
                                yz = yξ*zη - zξ*yη
                                zx = zξ*xη - xξ*zη                        
                                Ja = √(xy^2 + yz^2 + zx^2)      
                            elseif dim == 3 && setDim == 1  # line load on 3D solid or shell
                                Ja = √(xξ^2+yξ^2+zξ^2)
                            elseif dim == 3 && setDim == 0  # point load on 3D solid or shell
                                Ja = 1
                            elseif setDim == 2 && type == :axisymX  # volume load on axisymX                
                                Ja = √(xξ * yη - xη * yξ) * 2π * y
                            elseif setDim == 2 && type == :axisymY  # volume load on axisymY       
                                Ja = √(xξ * yη - xη * yξ) * 2π * x   
                            elseif dim == 2 && setDim == 2          # volume load on 2D solid
                                Ja = √(xξ * yη - xη * yξ) * width    
                            elseif setDim == 1 && type == :axisymX  # surface load on axisymX    
                                Ja = √(xξ^2+yξ^2) * 2π * y  
                            elseif setDim == 1 && type == :axisymY  # surface load on axisymY
                                Ja = √(xξ^2+yξ^2) * 2π * x
                            elseif dim == 2 && setDim == 1          # surface load on 2D solid
                                Ja = √(xξ^2+yξ^2) * width
                            elseif setDim == 0 && type == :axisymX  # line load on axisymX 
                                Ja = 2π * y
                            elseif setDim == 0 && type == :axisymX  # line load on axisymY
                                Ja = 2π * x
                            elseif dim == 2 && setDim == 0          # line load on 2D solid   
                                Ja = width 
                            end

                            f1 += AT * f0 * Ja * W[g] # elemental nodal load vector
                        end
                        f[nn2[e,:]] += f1 # inserting f1 to global nodal load vector
                    end
                end
            end

            # constructing bcs output:

            if td == :none             # if there is no amplitude variation
                bcs.fconst += f        # adding global nodal load vector to global const. load vector
            elseif td isa Function     # if amplitude variation is a function
                if isempty(bcs.fbaseA)  # if global base load vector for functions is empty
                    bcs.fbaseA = reshape(f, :, 1) # global base load matrix for functions = currnet base load vector for functions
                    bcs.fvarA = [td]   # global load variation matrix for functions = current time dependency function
                else                    # if  global base load vector for functions is not empty
                    bcs.fbaseA = [bcs.fbaseA reshape(f, :, 1)] # place current base load vector for functions next to current global base load matrix for functions
                    bcs.fvarA = append!(bcs.fvarA,[td]) # place current load variation function under current global load variation function matrix      
                end
            elseif td isa Matrix       # if amplitude variation is row matrix
                if isempty(bcs.fbaseB)  # if global base load matrix for numeric variation is empty
                    bcs.fbaseB = f      # global base load matrix for numeric variation = currnet base load vector for numeric variation
                    bcs.fvarB = td     # global load variation matrix for numeric variation = current numeric amplitude variation
                else                    # if  global base load vector for numeric variation is not empty
                    bcs.fbaseB = [bcs.fbaseB f] # place current base load vector for numeric variation next to current global base load matrix for numeric variation
                    bcs.fvarB = append!(bcs.fvarB,td) # place current numeric load variation under current global numeric load variation matrix      
                end      
            end
        end
    end
    return bcs
end           

function elasticSupport(mesh,esup,DoFs)

    sys  = SystemMatrices() # initialization of an empty system matrix structure
    bcs  = BoundaryConditions() # initialization of an empty boundary condition structure

    esup = to_vector(esup)  # if esup is not a vector/dictionary, convert to a vector
    xyz  = mesh.xyz          # coordinate matrix
    
    dof  = maximum(DoFs)    # number of degrees of freedom of the model
    dim  = size(DoFs,2)/3+1  # dimension of the model

    ridx  = Vector{Int64}()   # initializing row index vector
    cidx  = Vector{Int64}()   # initializing col index vector
    datak = Vector{Float64}() # initializing data vector
    datac = Vector{Float64}() # initializing data vector
    
    for (esupKey, esupData) in pairs(esup)  # for each spring set
        scopevec  = esupData.scope
        dirvec    = esupData.dir     # vector of direction definitions in the current displacement set
        kxvec     = esupData.kx
        kyvec     = esupData.ky
        kzvec     = esupData.kz
        cxvec     = esupData.cx
        czvec     = esupData.cy
        cxvec     = esupData.cz
        tdvec     = esupData.td      # vector of time dependency function in the current displacement set
        widthvec  = esupData.width

        n = length(scopevec)         # number of defined esup in the current set
        ii = 1    
        for i = 1:n
            scope = scopevec[i]
            dir   = dirvec[i]      # direction definition of the ith displacement
            kx    = kxvec[i]
            ky    = kyvec[i]
            kz    = kzvec[i]
            cx    = cxvec[i]
            cy    = cyvec[i]
            cz    = czvec[i]
            width = widthvec[i]

            scope = scope2selection(mesh,scope)  # convert scope to selection structure
            damp = false
            if cx != 0 || cy != 0 || cz != 0
                damp = true
            end

            nset    = scope.phnTags     # node tags in the set
            elems   = scope.elems       # elemet type => tags Dictionary
            setDim  = scope.dim         # dimension of the scoped geometry
            non     = length(nset)      # number of nodes of the ith esup

            pdimS = max(length(k),length(c)) # physical dimension of spring element

            if dim == 2
                pdim = 2
            elseif dim == 3
                pdim = 3
            end

            # initialization of a esup matrix
            if pdim == 3           # if physical dimension is 3 
                K0 = zeros(3,3)
                C0 = zeros(3,3)
            elseif pdim == 2       # if physical dimension is 2
                K0 = zeros(2,2)
                C0 = zeros(2,2)
            elseif pdim == 1       # if physical dimension is 1
                K0 = 0.0
                C0 = 0.0
            else
                error("dimension of the problem is $(problem.dim).")
            end

            for (eType, eTags) in elems # for each element type
                nnET    = mesh.elem[eType].nn    # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags  # element Tags of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix of selected elements
                nonpe   = mesh.elem[eType].non        # number of nodes per elements
                noe     = length(eTags);              # number of selected elements
                noqpe   = nonpe * pdim          # number of dofs per elements
                nov     = noqpe^2 * noe         # number of all individual stiffness values

                nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # properties of current element type:
                # rows: shape functions, colums: Gauss points
                HG  = mesh.elem[eType].HG[1]  # shape functions at Gauss points    
                HξG = mesh.elem[eType].HξG[1] # ∂h/∂ξ derivatives at Gauss points
                HηG = mesh.elem[eType].HηG[1] # ∂h/∂η derivatives at Gauss points
                HζG = mesh.elem[eType].HζG[1] # ∂h/∂ζ derivatives at Gauss points

                W   = mesh.elem[eType].GW[1]  # Gauss weigths      
                gpe = length(W)            # number of Gauss points per elements

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = xcrd[nn]      # elemental x coordinate matrix
                Y = ycrd[nn]      # elemental y coordinate matrix
                Z = zcrd[nn]      # elemental z coordinate matrix

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> Gauss points
                Xξ = X * HξG   # ∂x/∂ξ matrix
                Yξ = Y * HξG   # ∂y/∂ξ matrix
                Zξ = Z * HξG   # ∂z/∂ξ matrix
                Xη = X * HηG   # ∂x/∂η matrix
                Yη = Y * HηG   # ∂y/∂η matrix
                Zη = Z * HηG   # ∂z/∂η matrix
                Xζ = X * HζG   # ∂x/∂ζ matrix
                Yζ = Y * HζG   # ∂y/∂ζ matrix
                Zζ = Z * HζG   # ∂z/∂ζ matrix

                # Coordinates at Gauss points arranged in matrices:
                # rows -> elements, columns -> Gauss points
                XG = X * HG   # x coordinates at Gauss ponts of the elements
                YG = Y * HG   # y coordinates at Gauss ponts of the elements
                ZG = Z * HG   # z coordinates at Gauss ponts of the elements

                ridx1   = zeros(Int32, nov, 1)    # initialization of row index vector
                cidx1   = zeros(Int32, nov, 1)    # initialization of column index vector
                datak1   = zeros(Float64, nov, 1)  # initialization of data vector

                if damp
                    datac1   = zeros(Float64, nov, 1)  # initialization of data vector
                end

                # initialization of transpose of approximation matrix:
                A = zeros(Float64, pdim, nonpe * pdim)

                for e in 1:noe # for each selected element from current element type
                    Ke1 = zeros(Float64, pdim*nonpe, pdim*nonpe) # elemental nodal load vector ini
                    if damp
                        Ce1 = zeros(Float64, pdim*nonpe, pdim*nonpe)
                    end
                    for g in 1:gpe    # for each Gauss Point

                        # extracting Jacobi matrix elements:
                        xξ = Xξ[e,g]  # ∂x/∂ξ
                        yξ = Yξ[e,g]  # ∂y/∂ξ
                        zξ = Zξ[e,g]  # ∂z/∂ξ
                        xη = Xη[e,g]  # ∂x/∂η
                        yη = Yη[e,g]  # ∂y/∂η
                        zη = Zη[e,g]  # ∂z/∂η
                        xζ = Xζ[e,g]  # ∂x/∂ζ
                        yζ = Yζ[e,g]  # ∂y/∂ζ
                        zζ = Zζ[e,g]  # ∂z/∂ζ
                                    
                        # filling up aprroximation matrix:
                        for k ∈ 1:nonpe, l ∈ 1:pdim
                            AT[pdim*k-(pdim-l),l] = HG[k,g]
                        end   

                        # extracting coordinates & filling up esup matrices
                        x = XG[e,g] # x coordinate at Gauss point g of element e
                        K0[1,1] = isa(kx, Function) ? kx(x, y, z) : kx
                        if pdim > 1
                            y = YG[e,g] # y coordinate at Gauss point g of element e
                            K0[2,2] = isa(ky, Function) ? ky(x, y, z) : ky
                            if pdim > 2
                                z = ZG[e,g] # z coordinate at Gauss point g of element e
                                K0[3] = isa(kz, Function) ? kz(x, y, z) : kz
                            end
                        end

                        if damp
                            C0[1,1] = isa(cx, Function) ? cx(x, y, z) : cx
                            if pdim > 1
                                y = YG[e,g] # y coordinate at Gauss point g of element e
                                C0[2,2] = isa(cy, Function) ? cy(x, y, z) : cy
                                if pdim > 2
                                    z = ZG[e,g] # z coordinate at Gauss point g of element e
                                    C0[3] = isa(cz, Function) ? cz(x, y, z) : cz
                                end
                            end       
                        end                     

                        # transforming esup matrices
                        if dir == :gsys
                            K1 = K0
                            if damp
                                C1 = C0
                            end
                        end

                        # construction of Ja:
                        # volume support on 3D solid oe shell:
                        if dim == 3 && (setDim == 3 || (setDim == 2 && type == :shell))
                            Ja = xξ * (yη * zζ - yζ * zη) -
                                    yξ * (xη * zζ - xζ * zη) +
                                    zξ * (xη * yζ - xζ * yη)
                        elseif dim == 3 && setDim == 2  # surface support on 3D solid or shell
                            xy = xξ*yη - yξ*xη
                            yz = yξ*zη - zξ*yη
                            zx = zξ*xη - xξ*zη                        
                            Ja = √(xy^2 + yz^2 + zx^2)      
                        elseif dim == 3 && setDim == 1  # line support on 3D solid or shell
                            Ja = √(xξ^2+yξ^2+zξ^2)
                        elseif dim == 3 && setDim == 0  # point support on 3D solid or shell
                            Ja = 1
                        elseif setDim == 2 && type == :axisymX  # volume support on axisymX                
                            Ja = √(xξ * yη - xη * yξ) * 2π * y
                        elseif setDim == 2 && type == :axisymY  # volume support on axisymY       
                            Ja = √(xξ * yη - xη * yξ) * 2π * x   
                        elseif dim == 2 && setDim == 2          # volume support on 2D solid
                            Ja = √(xξ * yη - xη * yξ) * width    
                        elseif setDim == 1 && type == :axisymX  # surface support on axisymX    
                            Ja = √(xξ^2+yξ^2) * 2π * y  
                        elseif setDim == 1 && type == :axisymY  # surface support on axisymY
                            Ja = √(xξ^2+yξ^2) * 2π * x
                        elseif dim == 2 && setDim == 1          # surface support on 2D solid
                            Ja = √(xξ^2+yξ^2) * width
                        elseif setDim == 0 && type == :axisymX  # line support on axisymX 
                            Ja = 2π * y
                        elseif setDim == 0 && type == :axisymX  # line support on axisymY
                            Ja = 2π * x
                        elseif dim == 2 && setDim == 0          # line support on 2D solid   
                            Ja = width 
                        end  
                        
                        Ke1 += A' * K1 * A * Ja * W(g)
                        if damp
                            Ce1 += A' * C1 * A * Ja * W(g)
                        end
                    end

                    # filling up sparse matrix creator vectors:
                    for k1 in 1:noqpe      # for each dof per element
                        for k2 in 1:noqpe  # for each dof per element
                            ridxk1[ii] = nn2[e, k1]
                            cidxk1[ii] = nn2[e, k2]
                            datak1[ii] = Ke1[k1,k2]
                            if damp
                                datac1[ii] = Ce1[k1,k2]
                            end
                            ii += 1
                        end
                    end
                end
                append!(ridx,ridx1)
                append!(cidx,cidx1)
                append!(datak,datak1)
                if damp
                    append!(datac,datac1)
                end
            end
        end
    end
    sys.Ke = sparse(ridx, cidx, datak, dof, dof)  # creating the global stiffness matrix
    dropzeros!(sys.Ke)                           # drop zeros from stiffness matrix
    if damp
        sys.Ce = sparse(ridx, cidx, datac, dof, dof)  # creating the global stiffness matrix
        dropzeros!(sys.Ce)                           # drop zeros from stiffness matrix 
    end 

    # elastic dipslacement:

    for (esupKey, esupData) in pairs(esup)  # for each spring set
        scopevec  = esupData.scope
        dirvec    = esupData.dir     # vector of direction definitions in the current displacement set
        uxvec     = esupData.ux      # vector of displacements in x direction in the current displacement set
        uyvec     = esupData.uy      # vector of displacements in y direction in the current displacement set
        uzvec     = esupData.uz      # vector of displacements in z direction in the current displacement set
        rxvec     = esupData.rx      # vector of rotations about x axis in the current displacement set
        ryvec     = esupData.ry      # vector of rotations about y axis in the current displacement set
        rzvec     = esupData.rz      # vector of rotations about z axis in the current displacement set
        tdvec     = esupData.td      # vector of time dependency function in the current displacement set

        n = length(scopevec)         # number of defined esup in the current set  
        for i = 1:n
            scope = scopevec[i]
            dir   = dirvec[i]      # direction definition of the ith displacement
            ux    = uxvec[i]       # displacement in x direction
            uy    = uyvec[i]       # displacement in y direction
            uz    = uzvec[i]       # displacement in z direction
            rx    = rxvec[i]       # rotation about x axis
            ry    = ryvec[i]       # rotation about y axis
            rz    = rzvec[i]       # rotation about z axis
            td    = tdvec[i]       # time dependency function of the ith displacement  
            
            scope = scope2selection(mesh,scope)  # convert scope to selection structure
            nset  = scope.phnTags     # node tags in the set
            non   = length(nset)      # number of nodes of the ith esup            

            if dim == 2                 # if dimension is 2
                U = [ux,uy,rz]           # general displacement vector of a node
            elseif dim == 3             # if dimension is 3
                U = [ux,uy,uz,rx,ry,rz]  # general displacement vector of a node
            end

            Kset = Vector{Int32}(undef, 0)
            Kload = Vector{Float64}(undef, 0) 
            for j in eachindex(U)     # for each possible displacement coordinate
                if U[j] isa Number &&  U[j] != 0        # if the jth displacement coordinate is a number (but it is not zero)
                    append!(Kset, DoFs[nset,j])       # add the corresponding dof number to Kset dof vector    
                    append!(Kload,fill(U[j],non))                 # add the given displacement to Kload vevtor
                elseif U[j] isa Function        # if the jth displacement coordinate is a function
                    append!(Kset, DoFs[nset,j])       # add the corresponding dof number to Kset dof vector                     
                    xyz = mesh.xyz              # coordinate matrix
                    x   = xyz[nset,1]           # x coordinates
                    y   = xyz[nset,2]           # y coordinates
                    z   = xyz[nset,3]           # z coordinates
                    append!(Kload,U[j](x,y,z))  # add the numeric values (with substituted x,y,z coordinates) to Kload vector
                end
            end

            if td == :none && !isempty(Kset)  # if amplitude variation is not defined and Kset is not empty
                if isempty(bcs.fconst)
                    bcs.fconst = -Ke[:,Kset]*Kload
                else
                    bcs.fconst -= Ke[:,Kset]*Kload  # add kinematic load to constant nodal load vector
                end    
            elseif td isa Function && !isempty(Kset)  # if amplitude variation is a function and Kset is not empty
                if isempty(bcs.fbaseA)         # if global base load vector for functions is empty
                    bcs.fbaseA = -Ke[:,Kset]*Kload  # global base load matrix for functions = currnet base load vector for functions
                else                           # if global base load vector for functions is not empty
                    bcs.fbaseA = [bcs.fbaseA -Ke[:,Kset]*Kload] # place current base load vector for functions next to current global base load matrix for functions
                end
            elseif td isa Matrix && !isempty(Kset)  # if amplitude variation is a matrix and Kset is not empty
                if isempty(bcs.fbaseB)   # if global base load vector for numeric variation is empty
                    bcs.fbaseB = -Ke[:,Kset]*Kload # global base load matrix for numeric variation = currnet base load vector for numeric variation
                else                     # if global base load vector for numeric variation is not empty
                    bcs.fbaseB = [bcs.fbaseB -Ke[:,Kset]*Kload] # place current base load vector for numeric variation next to current global base load matrix for functions
                end
            end
        end
    end

return sys, bcs      
    
end

function initialDisplacement(mesh,idisp,DoFs)
    dof  = maximum(DoFs)    # number of degrees of freedom of the model
    dim = size(DoFs,2)/3+1  # dimension of the model
    Aset = collect(1:dof)
    q0 = zeros(dof,1)

    idisp = to_vector(idisp)      # if initial displacement is not a vector/dictionary, convert to a vector   
    
    for (idispKey, idispData) in pairs(idisp)  # for each displacement set
        scopevec  = idispData.scope   # vector of scoped geometry identifiers in the current displacement set
        dirvec    = idispData.dir     # vector of direction definitions in the current displacement set
        uxvec     = idispData.ux      # vector of displacements in x direction in the current displacement set
        uyvec     = idispData.uy      # vector of displacements in y direction in the current displacement set
        uzvec     = idispData.uz      # vector of displacements in z direction in the current displacement set
        rxvec     = idispData.rx      # vector of rotations about x axis in the current displacement set
        ryvec     = idispData.ry      # vector of rotations about y axis in the current displacement set
        rzvec     = idispData.rz      # vector of rotations about z axis in the current displacement set

        n = length(scopevec)         # number of defined displacements in the current displacement set   
        
        for i = 1:n                  # for each displacement in the current displacement set
            scope   = scopevec[i]    # scoped geometry of the ith displacement
            dir     = dirvec[i]      # direction definition of the ith displacement
            ux      = uxvec[i]       # displacement in x direction
            uy      = uyvec[i]       # displacement in y direction
            uz      = uzvec[i]       # displacement in z direction
            rx      = rxvec[i]       # rotation about x axis
            ry      = ryvec[i]       # rotation about y axis
            rz      = rzvec[i]       # rotation about z axis

            scope = scope2selection(mesh,scope)  # convert scope to selection structure  
            
            nset = scope.phnTags          # node ids belonging to the scoped geometry
            non  = length(nset)           # number of nodes belonging to the scoped geometry

            if dim == 2                  # if dimension is 2
                U = [ux,uy,rz]           # general displacement vector of a node
            elseif dim == 3              # if dimension is 3
                U = [ux,uy,uz,rx,ry,rz]  # general displacement vector of a node
            end

            for j in eachindex(U)     # for each possible displacement coordinate
                if U[j] isa Number          # if the jth displacement coordinate is a number (but it is not zero)
                    q0[DoFs[nset,j]] = fill(U[j],non)
                elseif U[j] isa Function        # if the jth displacement coordinate is a function
                    xyz = mesh.xyz              # coordinate matrix
                    x   = xyz[nset,1]           # x coordinates
                    y   = xyz[nset,2]           # y coordinates
                    z   = xyz[nset,3]           # z coordinates
                    q0[DoFs[nset,j]] = U[j](x,y,z)
                end
            end
        end
    end 
    return q0                      
end

function initialVelocity(mesh,ivelo,DoFs)
    dof  = maximum(DoFs)    # number of degrees of freedom of the model
    dim = size(DoFs,2)/3+1  # dimension of the model
    Aset = collect(1:dof)
    v0 = zeros(dof,1)

    ivelo = to_vector(ivelo)      # if initial displacement is not a vector/dictionary, convert to a vector   
    
    for (iveloKey, iveloData) in pairs(ivelo)  # for each displacement set
        scopevec  = iveloData.scope   # vector of scoped geometry identifiers in the current displacement set
        dirvec    = iveloData.dir     # vector of direction definitions in the current displacement set
        vxvec     = iveloData.vx      # vector of displacements in x direction in the current displacement set
        vyvec     = iveloData.vy      # vector of displacements in y direction in the current displacement set
        vzvec     = iveloData.vz      # vector of displacements in z direction in the current displacement set
        wxvec     = iveloData.wx      # vector of rotations about x axis in the current displacement set
        wyvec     = iveloData.wy      # vector of rotations about y axis in the current displacement set
        wzvec     = iveloData.wz      # vector of rotations about z axis in the current displacement set

        n = length(scopevec)         # number of defined displacements in the current displacement set   
        
        for i = 1:n                  # for each displacement in the current displacement set
            scope   = scopevec[i]    # scoped geometry of the ith displacement
            dir     = dirvec[i]      # direction definition of the ith displacement
            vx      = vxvec[i]       # displacement in x direction
            vy      = vyvec[i]       # displacement in y direction
            vz      = vzvec[i]       # displacement in z direction
            wx      = wxvec[i]       # rotation about x axis
            wy      = wyvec[i]       # rotation about y axis
            wz      = wzvec[i]       # rotation about z axis

            scope = scope2selection(mesh,scope)  # convert scope to selection structure  
            
            nset = scope.phnTags          # node ids belonging to the scoped geometry
            non  = length(nset)           # number of nodes belonging to the scoped geometry

            if dim == 2                  # if dimension is 2
                V = [vx,vy,vz]           # general displacement vector of a node
            elseif dim == 3              # if dimension is 3
                Val = [vx,vy,vz,wx,wy,wz]  # general displacement vector of a node
            end

            for j in eachindex(V)     # for each possible displacement coordinate
                if V[j] isa Number          # if the jth displacement coordinate is a number (but it is not zero)
                    v0[DoFs[nset,j]] = fill(V[j],non)
                elseif V[j] isa Function        # if the jth displacement coordinate is a function
                    xyz = mesh.xyz              # coordinate matrix
                    x   = xyz[nset,1]           # x coordinates
                    y   = xyz[nset,2]           # y coordinates
                    z   = xyz[nset,3]           # z coordinates
                    v0[DoFs[nset,j]] = V[j](x,y,z)
                end
            end
        end
    end 
    return v0                      
end

# A function creating stiffness matrix ---------------------------------------------------------------------------
function stiffnessMatrix(mesh,sec,dofTable)

    xyz = mesh.xyz           # coordinate matrix
    x   = xyz[:,1]           # x coordinates
    y   = xyz[:,2]           # y coordinates
    z   = xyz[:,3]           # z coordinates

    non = length(x)          # number of nodes
    dof = maximum(dofTable)  # number of degrees of freedom of the model

    ridx = Vector{Int64}()   # initializing row index vector
    cidx = Vector{Int64}()   # initializing col index vector
    data = Vector{Float64}() # initializing data vector

    for (secKey, secData) in pairs(sec)  # for each section

        scope   = secData.scope      # scope of the current section
        D       = secData.D          # material matrix
        dim     = secData.dim        # spatial dimension of the current section
        pdim    = secData.pdim       # physical dimension of the current section
        intid   = secData.intid      # integration rule identifier of the current section

        setData = scope2selection(mesh,scope)  # convert scope to selection structure      

        if secData isa Solid3D                   # 3D Solid section
            for (eType, eTags) in setData.elems  # for each element type
                nnET    = mesh.elem[eType].nn    # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags  # element tags of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix of current section with curent element type
                nonpe   = mesh.elem[eType].non   # number of nodes per elements
                noe     = length(eTags);         # number of elements
                noqpe   = nonpe * pdim           # number of dofs per elements
                nov     = noqpe^2 * noe          # number of all individual stiffness values

                nn2 = zeros(Int32, noe, nonpe*pdim)   # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # properties of current element type:
                # rows: shape functions, colums: Gauss points
                HξG = mesh.elem[eType].HξG[intid]  # ∂h/∂ξ derivatives at Gauss points   
                HηG = mesh.elem[eType].HηG[intid]  # ∂h/∂η derivatives at Gauss points   
                HζG = mesh.elem[eType].HζG[intid]  # ∂h/∂ζ derivatives at Gauss points   
                W   = mesh.elem[eType].GW[intid]   # Gauss weights

                gpe = length(W)             # number of Gauss points per elements

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix
                Z = z[nn]  # elemental z coordinate matrix

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> Gauss points
                Xξ = X * HξG # ∂x/∂ξ matrix
                Yξ = Y * HξG # ∂y/∂ξ matrix
                Zξ = Z * HξG # ∂z/∂ξ matrix
                Xη = X * HηG # ∂x/∂η matrix
                Yη = Y * HηG # ∂y/∂η matrix
                Zη = Z * HηG # ∂z/∂η matrix
                Xζ = X * HζG # ∂x/∂ζ matrix
                Yζ = Y * HζG # ∂y/∂ζ matrix
                Zζ = Z * HζG # ∂z/∂ζ matrix

                ridx1   = zeros(Int32, nov, 1)   # initialization of row index vector
                cidx1   = zeros(Int32, nov, 1)   # initialization of column index vector
                data1   = zeros(Float64, nov, 1) # initialization of data vector
                B       = zeros(Float64, 6, nonpe * pdim) # initialization of B matrix
                K1      = zeros(Float64, pdim*nonpe, pdim*nonpe)  # initialization of stiffness matrix for eth element

                DB = similar(B)             # Holds D * B (size: size(D,1) × size(B,2))
                BTDB = zeros(size(K1))      # Holds B' * (D * B) (size: pdim*nonpe × pdim*nonpe)

                i = 1
                for e in 1:noe                   # for each element
                    fill!(K1, 0.0)  # resetting stiffness matrix for the eth element

                    for g in 1:gpe               # for each Gauss point

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        xξ = Xξ[e,g]  # ∂x/∂ξ at the gth point of eth element
                        yξ = Yξ[e,g]  # ∂y/∂ξ at the gth point of eth element
                        zξ = Zξ[e,g]  # ∂z/∂ξ at the gth point of eth element
                        xη = Xη[e,g]  # ∂x/∂η at the gth point of eth element
                        yη = Yη[e,g]  # ∂y/∂η at the gth point of eth element
                        zη = Zη[e,g]  # ∂z/∂η at the gth point of eth element
                        xζ = Xζ[e,g]  # ∂x/∂ζ at the gth point of eth element
                        yζ = Yζ[e,g]  # ∂x/∂ζ at the gth point of eth element
                        zζ = Zζ[e,g]  # ∂x/∂ζ at the gth point of eth element

                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)

                        # assembling B matrix at the gth point of the eth element:
                        for k = 1:nonpe  # for each node of the eth element
                            hx = HξG[k,g] * ξx + HηG[k,g] * ηx + HζG[k,g] * ζx # ∂h/∂x 
                            hy = HξG[k,g] * ξy + HηG[k,g] * ηy + HζG[k,g] * ζy # ∂h/∂y
                            hz = HξG[k,g] * ξz + HηG[k,g] * ηz + HζG[k,g] * ζz # ∂h/∂z
                            B[1,3*k-2] = B[4,3*k-1] = B[6,3*k]   = hx  # filling up B matrix
                            B[2,3*k-1] = B[4,3*k-2] = B[5,3*k]   = hy  # filling up B matrix
                            B[3,3*k]   = B[5,3*k-1] = B[6,3*k-2] = hz  # filling up B matrix
                        end

                        # adding the stiffness part form the gth point to the eth stiffness matrix:
                         mul!(DB, D, B)                     # temp1 = D * B
                         mul!(BTDB, transpose(B), DB)     # temp2 = B' * (D * B)
                         @. K1 += BTDB * det * W[g]      # final update
                    end

                    # filling up sparse matrix creator vectors:
                    for k1 in 1:noqpe      # for each dof per element
                        for k2 in 1:noqpe  # for each dof per element
                            ridx1[i] = nn2[e, k1]
                            cidx1[i] = nn2[e, k2]
                            data1[i] = K1[k1,k2]
                            i += 1
                        end
                    end
                end
                append!(ridx,ridx1) # adding row indexer vector part from the eth element
                append!(cidx,cidx1) # adding column indexer vector part from the eth element
                append!(data,data1) # adding data vector part from the eth element
            end
        elseif secData isa Solid2D  # 2D Solid section
            secType = secData.type       # type of the 2D Solid section
            for (eType, eTags) in setData.elems # for each element type
                nnET    = mesh.elem[eType].nn   # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags # coonectivity matrix of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix
                nonpe   = mesh.elem[eType].non  # number of nodes per elements
                noe     = length(eTags)         # number of nodes
                noqpe   = nonpe * pdim          # number of dofs per elements
                nov     = noqpe^2 * noe         # number of all individual stiffness values

                nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # properties of current element type:
                # rows: shape functions, colums: Gauss points
                HG  = mesh.elem[eType].HG[intid]  # shape function values at Gauss points
                HξG = mesh.elem[eType].HξG[intid] # ∂h/∂ξ derivatives at Gauss points   
                HηG = mesh.elem[eType].HηG[intid] # ∂h/∂η derivatives at Gauss points   
                W   = mesh.elem[eType].GW[intid]  # Gauss weigths  
                
                gpe = length(W) # number of Gauss points per elements

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix             

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> Gauss points
                Xξ = X * HξG    # ∂x/∂ξ matrix
                Yξ = Y * HξG    # ∂y/∂ξ matrix
                Xη = X * HηG    # ∂x/∂η matrix
                Yη = Y * HηG    # ∂y/∂η matrix

                ridx1   = zeros(Int32, nov, 1)    # initialization of row index vector
                cidx1   = zeros(Int32, nov, 1)    # initialization of column index vector
                data1   = zeros(Float64, nov, 1)  # initialization of data vector
              
                # initialization of B matrices and integration constants based on mechanical model:
                if secType == :PlaneStrain  # planestrain model
                    B = zeros(Float64, 3, nonpe * pdim) # initialization of B matrix
                    C = 1
                    typeid = 1 # assigning typeid = 1 to this mechanical model
                elseif secType == :PlaneStress # planestress model
                    B = zeros(Float64, 3, nonpe * pdim) # initialization of B matrix
                    C = secData.width  # width of the current section
                    typeid = 2 # assigning typeid = 2 to this mechanical model
                elseif secType == :AxiSymmetricY  # axisymmetric model about axis y
                    B = zeros(Float64, 4, nonpe * pdim) # initialization of B matrix
                    R = X * HG  # radial coordinate matrix at the Gauss points
                    d = 1
                    typeid = 3  # assigning typeid = 3 to this mechanical model
                elseif secType == :AxiSymmetricX  # axisymmetric model about axis x
                    B = zeros(Float64, 4, nonpe * pdim) # initialization of B matrix
                    R = Y * HG  # radial coordinate matrix at the Gauss points
                    d = 0
                    typeid = 4  # assigning typeid = 4 to this mechanical model
                end
                K1 = zeros(Float64, pdim*nonpe, pdim*nonpe) # initialization of stiffness matrix for eth element   
                DB = similar(B)             # Holds D * B (size: size(D,1) × size(B,2))
                BTDB = zeros(size(K1))      # Holds B' * (D * B) (size: pdim*nonpe × pdim*nonpe)

                i = 1
                @inbounds for e in 1:noe  # for each element
                    # K1 = zeros(Float64, pdim*nonpe, pdim*nonpe) # initialization of stiffness matrix for eth element 
                    fill!(K1, 0.0)  # resetting stiffness matrix for the eth element  
                    @inbounds for g in 1:gpe     # for each Gauss point

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        xξ = Xξ[e,g]  # ∂x/∂ξ at the gth point of eth element
                        yξ = Yξ[e,g]  # ∂y/∂ξ at the gth point of eth element
                        xη = Xη[e,g]  # ∂x/∂η at the gth point of eth element
                        yη = Yη[e,g]  # ∂y/∂η at the gth point of eth element

                        # Jacobi matrix determinant at the gth point of the eth element
                        det = xξ * yη - xη * yξ

                        # elements of the Jacobi inverse at the gth point of the eth element:
                        ξx =  yη / det # ∂ξ/∂x at the gth point of eth element
                        ηx = -yξ / det # ∂η/∂x at the gth point of eth element
                        ξy = -xη / det # ∂ξ/∂y at the gth point of eth element
                        ηy =  xξ / det # ∂η/∂y at the gth point of eth element    

                        # assembling B matrix at the gth point of the eth element:
                        if typeid <= 2  # if mechanical model is planestrain or planestress
                            @inbounds for k = 1:nonpe  # for each node of the eth element
                                hx = HξG[k,g] * ξx + HηG[k,g] * ηx # ∂h/∂x 
                                hy = HξG[k,g] * ξy + HηG[k,g] * ηy # ∂h/∂y 
                                B[1,2*k-1] = B[3,2*k] = hx  # filling up B matrix with ∂h/∂x values
                                B[2,2*k] = B[3,2*k-1] = hy  # filling up B matrix with ∂h/∂y values
                            end
                        elseif typeid >=3  # if mechanical model is axisymmetric
                            r = R[e,g]     # radial coordinate of the gth point of the eth element
                            C = 2 * r * π  # integration constant (length of the corresponding circle)
                            @inbounds for k = 1:nonpe # for each node of the eth element
                                hx = HξG[k,g] * ξx + HηG[k,g] * ηx # ∂h/∂x 
                                hy = HξG[k,g] * ξy + HηG[k,g] * ηy # ∂h/∂y                                  
                                B[1,2*k-1] = B[4,2*k] = hx # filling up B matrix with ∂h/∂x values
                                B[2,2*k] = B[4,2*k-1] = hy # filling up B matrix with ∂h/∂y values
                                B[3,2*k-d] = HG[k,g]/r     # filling up B matrix with h/r values
                            end
                        end

                        # adding the stiffness part form the gth point to the eth stiffness matrix:
                        mul!(DB, D, B)                     
                        mul!(BTDB, transpose(B), DB)     
                        @. K1 += BTDB * C * det * W[g]     
                    end

                    # filling up sparse matrix creator vectors:
                    @inbounds for k1 in 1:noqpe # for each dof per element
                        @inbounds for k2 in 1:noqpe # for each dof per element
                            ridx1[i] = nn2[e, k1]
                            cidx1[i] = nn2[e, k2]
                            data1[i] = K1[k1,k2]
                            i += 1
                        end
                    end
                end
                append!(ridx,ridx1)
                append!(cidx,cidx1)
                append!(data,data1)
            end
        elseif secData isa Shell  # Shell section

            t    = secData.width  # width of the current shell section
            intm = secData.intm   # integration rule identifier for membrane+bending part
            ints = secData.ints   # integration rule identifier for shear part
            intt = secData.intt   # integration rule identifier for torsional part
            kt   = secData.kt     # torsional defomration factor
            G    = secData.mat.G  # shear modulus

            # preparing material matrices: 
            Dm = zeros(6,6) # initialization of material matrix for membrane+bending part
            Dm[1:3,1:3] = 2*D[1:3,1:3] 
            Dm[4:6,4:6] = (2/3)*D[1:3,1:3] # material matrix for membrane+bending part

            Ds = zeros(4,4) # initialization of material matrix for shear part
            Ds[1:2,1:2] = 2*D[4:5,4:5]
            Ds[3:4,3:4] = (2/3)*D[4:5,4:5] # material matrix for shear part

            for (eType, eTags) in setData.elems # for each element type of the current shell section
                nnET    = mesh.elem[eType].nn   # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags # coonectivity matrix of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix
                nonpe   = mesh.elem[eType].non  # number of nodes per elements
                noe     = length(eTags)         # number of elements
                noqpe   = nonpe * pdim          # number of dofs per elements
                nov     = noqpe^2 * noe         # number of all individual stiffness values  

                ridx1 = zeros(Int32, nov, 1)    # initialization of row index vector
                cidx1 = zeros(Int32, nov, 1)    # initialization of column index vector
                data1 = zeros(Float64, nov, 1)  # initialization of data vector

                nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # elemental coordinate matrices:
                # rows -> elements (e1,e2,e3...), columns -> nodes of the element (n1,n2,n3...)
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix  
                Z = z[nn]  # elemental z coordinate matrix  
                
                if intm > 0
                    # -> for membrane part of the stiffness matrix ----------------------------------------------------------
                    HGm  = mesh.elem[eType].HG[intm]   # shape function values at Gauss points for membrane part   
                    HξGm = mesh.elem[eType].HξG[intm]  # ∂h/∂ξ derivatives at Gauss points for membrane part   
                    HηGm = mesh.elem[eType].HηG[intm]  # ∂h/∂η derivatives at Gauss points for membrane part    
                    HζGm = mesh.elem[eType].HζG[intm]  # ∂h/∂ζ derivatives at Gauss points for membrane part    
                    Wm   = mesh.elem[eType].GW[intm]   # Gauss weights at Gauss points for membrane part
                    gpem = length(Wm)                  # number of Gauss points per elements for membrane part

                    # x,y,z derivatives with respect to ξ,η,ζ (Jacobi elements) at Gauss points for membrane part:
                    # rows -> elements (e1,e2,e3...), columns -> Gauss points (G1,G2,G3...)
                    XξGm = X * HξGm # ∂x/∂ξ matrix
                    YξGm = Y * HξGm # ∂y/∂ξ matrix
                    ZξGm = Z * HξGm # ∂z/∂ξ matrix
                    XηGm = X * HηGm # ∂x/∂η matrix
                    YηGm = Y * HηGm # ∂y/∂η matrix
                    ZηGm = Z * HηGm # ∂z/∂η matrix       
                else
                    gpem = 0         
                end 
                
                if ints>0
                    # -> for shear part of the stiffness matrix ----------------------------------------------------------
                    if ints == intm  # if the id for shear part is the same as the id for membrane part
                        HGs  = HGm   # -> reference the membrane part matrices
                        HξGs = HξGm  
                        HηGs = HηGm
                        HζGs = HζGm
                        Ws   = Wm 
                        gpes = gpem
                        XξGs = XξGm # ∂x/∂ξ matrix
                        YξGs = YξGm # ∂y/∂ξ matrix
                        ZξGs = ZξGm # ∂z/∂ξ matrix
                        XηGs = XηGm # ∂x/∂η matrix
                        YηGs = YηGm # ∂y/∂η matrix
                        ZηGs = ZηGm # ∂z/∂η matrix 
                    else # if the id for shear part is different from the id for membrane part
                        HGs  = mesh.elem[eType].HG[ints]  # shape function values at Gauss points for shear part 
                        HξGs = mesh.elem[eType].HξG[ints] # ∂h/∂ξ derivatives at Gauss points for shear part   
                        HηGs = mesh.elem[eType].HηG[ints] # ∂h/∂η derivatives at Gauss points for shear part   
                        HζGs = mesh.elem[eType].HζG[ints] # ∂h/∂ζ derivatives at Gauss points for shear part   
                        Ws   = mesh.elem[eType].GW[ints]  # Gauss weights at Gauss points for shear part
                        gpes = length(Ws)                 # number of Gauss points per elements for shear part
                        XξGs = X * HξGs # ∂x/∂ξ matrix
                        YξGs = Y * HξGs # ∂y/∂ξ matrix
                        ZξGs = Z * HξGs # ∂z/∂ξ matrix
                        XηGs = X * HηGs # ∂x/∂η matrix
                        YηGs = Y * HηGs # ∂y/∂η matrix
                        ZηGs = Z * HηGs # ∂z/∂η matrix                  
                    end
                else
                    gpes = 0
                end

                if intt > 0
                    # -> for torsional part of the stiffness matrix ----------------------------------------------------------
                    if intt == intm # if the id for torsional part is the same as the id for membrane part
                        HGt  = HGm  # -> reference the membrane part matrices
                        HξGt = HξGm 
                        HηGt = HηGm
                        HζGt = HζGm
                        Wt   = Wm 
                        gpet = gpem
                        XξGt = XξGm # ∂x/∂ξ matrix
                        YξGt = YξGm # ∂y/∂ξ matrix
                        ZξGt = ZξGm # ∂z/∂ξ matrix
                        XηGt = XηGm # ∂x/∂η matrix
                        YηGt = YηGm # ∂y/∂η matrix
                        ZηGt = ZηGm # ∂z/∂η matrix   
                    elseif intt == ints # if the id for torsional part is the same as the id for shear part
                        HGt  = HGs      # -> reference the shear part matrices
                        HξGt = HξGs 
                        HηGt = HηGs
                        HζGt = HζGs
                        Wt   = Ws 
                        gpet = gpes 
                        XξGt = XξGs # ∂x/∂ξ matrix
                        YξGt = YξGs # ∂y/∂ξ matrix
                        ZξGt = ZξGs # ∂z/∂ξ matrix
                        XηGt = XηGs # ∂x/∂η matrix
                        YηGt = YηGs # ∂y/∂η matrix
                        ZηGt = ZηGs # ∂z/∂η matrix                       
                    else  # if the id for torsional part is different from the id for membrane or shear part
                        HGt  = mesh.elem[eType].HG[intt]  # shape function values at Gauss points for torsional part 
                        HξGt = mesh.elem[eType].HξG[intt] # ∂h/∂ξ derivatives at Gauss points for torsional part   
                        HηGt = mesh.elem[eType].HηG[intt] # ∂h/∂η derivatives at Gauss points for torsional part  
                        HζGt = mesh.elem[eType].HζG[intt] # ∂h/∂ζ derivatives at Gauss points for torsional part
                        Wt   = mesh.elem[eType].GW[intt]  # Gauss weights at Gauss points for torsional part
                        gpet = length(Wt)                 # number of Gauss points per elements for torsional part
                        XξGt = X * HξGt # ∂x/∂ξ matrix
                        YξGt = Y * HξGt # ∂y/∂ξ matrix
                        ZξGt = Z * HξGt # ∂z/∂ξ matrix
                        XηGt = X * HηGt # ∂x/∂η matrix
                        YηGt = Y * HηGt # ∂y/∂η matrix
                        ZηGt = Z * HηGt # ∂z/∂η matrix                    
                    end
                else
                    gpet = 0 
                end 
                
                # initialization of matrices:
                Bm = zeros(Float64, 6, nonpe*pdim) # B matrix for membrane+bending part
                Bs = zeros(Float64, 4, nonpe*pdim) # B matrix for shear part
                Bt = zeros(Float64, 1, nonpe*pdim) # B matrix for torsional part   
                K1 = zeros(Float64, pdim*nonpe, pdim*nonpe) # initialization of stiffness matrix for eth element  
                DmBm = zeros(size(Bm))
                BmTDmBm = zeros(size(K1))
                DsBs = zeros(size(Bs))
                BsTDsBs = zeros(size(K1))
                BtTBt = zeros(size(K1))
                
                i = 1
                @inbounds for e in 1:noe  # for each element of the current element type from current section
                    fill!(K1, 0.0)  # resetting stiffness matrix for the eth element  
                
                    @inbounds for g in 1:gpem  # for each Gauss point of the membrane part

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        # -> tangent vectors at the gth point of eth element
                        xξ = XξGm[e,g]  # ∂x/∂ξ at the gth point of eth element (a1x)
                        yξ = YξGm[e,g]  # ∂y/∂ξ at the gth point of eth element (a1y)
                        zξ = ZξGm[e,g]  # ∂z/∂ξ at the gth point of eth element (a1z)
                        xη = XηGm[e,g]  # ∂x/∂η at the gth point of eth element (a2x)
                        yη = YηGm[e,g]  # ∂y/∂η at the gth point of eth element (a2y)
                        zη = ZηGm[e,g]  # ∂z/∂η at the gth point of eth element (a2z)

                        # generating local coordinate system at the gth point of the eth element:
                        e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)
                        xζ = e3x .*0.5 .* t # ∂x/∂η at the gth point of eth element
                        yζ = e3y .*0.5 .* t # ∂y/∂η at the gth point of eth element
                        zζ = e3z .*0.5 .* t # ∂z/∂η at the gth point of eth element

                        # Jacobi determinant and elements of inverse at the gth point of the eth element:
                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)
                                                    
                        @inbounds for n in 1:nonpe      # for each node of the eth element

                            # shape function and it's derivatives with respect to x,y,z:
                            hx = HξGm[n,g] * ξx + HηGm[n,g] * ηx + HζGm[n,g] * ζx # ∂h/∂x 
                            hy = HξGm[n,g] * ξy + HηGm[n,g] * ηy + HζGm[n,g] * ζy # ∂h/∂y
                            hz = HξGm[n,g] * ξz + HηGm[n,g] * ηz + HζGm[n,g] * ζz # ∂h/∂z  

                            # gradient operators:
                            Lh1 = e1x*hx + e1y*hy + e1z*hz
                            Lh2 = e2x*hx + e2y*hy + e2z*hz

                            # elements of the Bm matrix associated with the strain contribution
                            # of the inplane displacements at the nth node
                            B1m11 = Lh1*e1x
                            B1m12 = Lh1*e1y
                            B1m13 = Lh1*e1z
                            B1m21 = Lh2*e2x
                            B1m22 = Lh2*e2y
                            B1m23 = Lh2*e2z
                            B1m31 = Lh2*e1x + Lh1*e2x
                            B1m32 = Lh2*e1y + Lh1*e2y
                            B1m33 = Lh2*e1z + Lh1*e2z

                            # elements of the Bm matrix associated with the strain contribution
                            # of the rotations at the nth node (also includes curvature effect)
                            B3m11 = 0.5*t*(B1m13*e3y - B1m12*e3z)
                            B3m12 = 0.5*t*(B1m11*e3z - B1m13*e3x)
                            B3m13 = 0.5*t*(B1m12*e3x - B1m11*e3y)
                            B3m21 = 0.5*t*(B1m23*e3y - B1m22*e3z)
                            B3m22 = 0.5*t*(B1m21*e3z - B1m23*e3x)
                            B3m23 = 0.5*t*(B1m22*e3x - B1m21*e3y)
                            B3m31 = 0.5*t*(B1m33*e3y - B1m32*e3z)
                            B3m32 = 0.5*t*(B1m31*e3z - B1m33*e3x)
                            B3m33 = 0.5*t*(B1m32*e3x - B1m31*e3y)

                            # populating Bm matrix:
                            Bm[1,(n-1)*pdim+1] = B1m11
                            Bm[1,(n-1)*pdim+2] = B1m12
                            Bm[1,(n-1)*pdim+3] = B1m13
                            Bm[2,(n-1)*pdim+1] = B1m21
                            Bm[2,(n-1)*pdim+2] = B1m22
                            Bm[2,(n-1)*pdim+3] = B1m23  
                            Bm[3,(n-1)*pdim+1] = B1m31
                            Bm[3,(n-1)*pdim+2] = B1m32
                            Bm[3,(n-1)*pdim+3] = B1m33     
                            Bm[4,(n-1)*pdim+4] = B3m11
                            Bm[4,(n-1)*pdim+5] = B3m12
                            Bm[4,(n-1)*pdim+6] = B3m13
                            Bm[5,(n-1)*pdim+4] = B3m21
                            Bm[5,(n-1)*pdim+5] = B3m22
                            Bm[5,(n-1)*pdim+6] = B3m23
                            Bm[6,(n-1)*pdim+4] = B3m31
                            Bm[6,(n-1)*pdim+5] = B3m32
                            Bm[6,(n-1)*pdim+6] = B3m33  
                        end
                        # updating stiffness matrix with the membrane+bending part at the gth Gauss point
                        mul!(DmBm, Dm, Bm)                 
                        mul!(BmTDmBm, transpose(Bm), DmBm)    
                        @. K1 += BmTDmBm * det * Wm[g]
                    end

                    @inbounds for g in 1:gpes  # for each Gauss point of the shear part

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        # -> tangent vectors at the gth point of eth element
                        xξ = XξGs[e,g]  # ∂x/∂ξ at the gth point of eth element (a1x)
                        yξ = YξGs[e,g]  # ∂y/∂ξ at the gth point of eth element (a1y)
                        zξ = ZξGs[e,g]  # ∂z/∂ξ at the gth point of eth element (a1z)
                        xη = XηGs[e,g]  # ∂x/∂η at the gth point of eth element (a2x)
                        yη = YηGs[e,g]  # ∂y/∂η at the gth point of eth element (a2y)
                        zη = ZηGs[e,g]  # ∂z/∂η at the gth point of eth element (a2z)

                        # generating local coordinate system at the gth point of the eth element:
                        e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)
                        xζ = e3x .*0.5 .* t # ∂x/∂η at the gth point of eth element
                        yζ = e3y .*0.5 .* t # ∂y/∂η at the gth point of eth element
                        zζ = e3z .*0.5 .* t # ∂z/∂η at the gth point of eth element

                        # Jacobi determinant and elements of inverse at the gth point of the eth element:
                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)
                                                    
                        @inbounds for n in 1:nonpe      # for each node of the eth element

                            # shape function and it's derivatives with respect to x,y,z:
                            h  = HGs[n,g] # shape function value at the nth node
                            hx = HξGs[n,g] * ξx + HηGs[n,g] * ηx + HζGs[n,g] * ζx # ∂h/∂x 
                            hy = HξGs[n,g] * ξy + HηGs[n,g] * ηy + HζGs[n,g] * ζy # ∂h/∂y
                            hz = HξGs[n,g] * ξz + HηGs[n,g] * ηz + HζGs[n,g] * ζz # ∂h/∂z 

                            # gradient operators:
                            Lh1 = e1x*hx + e1y*hy + e1z*hz
                            Lh2 = e2x*hx + e2y*hy + e2z*hz
                            Lh3 = e3x*hx + e3y*hy + e3z*hz
                            Lζ1 = e1x*ζx + e1y*ζy + e1z*ζz
                            Lζ2 = e2x*ζx + e2y*ζy + e2z*ζz
                            Lζ3 = e3x*ζx + e3y*ζy + e3z*ζz

                            # preparing θᵀΦ product:
                            θTΦ11 = e1z*e3y - e1y*e3z
                            θTΦ12 = e1x*e3z - e1z*e3x
                            θTΦ13 = e1y*e3x - e1x*e3y
                            θTΦ21 = e2z*e3y - e2y*e3z
                            θTΦ22 = e2x*e3z - e2z*e3x
                            θTΦ23 = e2y*e3x - e2x*e3y
                            θTΦ31 = e3z*e3y - e3y*e3z
                            θTΦ32 = e3x*e3z - e3z*e3x
                            θTΦ33 = e3y*e3x - e3x*e3y
                            
                            # calculating elements of the Bs matrix:
                            B1s11 = Lh3*e1x + Lh1*e3x
                            B1s12 = Lh3*e1y + Lh1*e3y
                            B1s13 = Lh3*e1z + Lh1*e3z
                            B1s21 = Lh3*e2x + Lh2*e3x
                            B1s22 = Lh3*e2y + Lh2*e3y
                            B1s23 = Lh3*e2z + Lh2*e3z

                            B2s11 = 0.5*t*h*(Lζ3*θTΦ11 + Lζ1*θTΦ31)
                            B2s12 = 0.5*t*h*(Lζ3*θTΦ12 + Lζ1*θTΦ32)
                            B2s13 = 0.5*t*h*(Lζ3*θTΦ13 + Lζ1*θTΦ33)
                            B2s21 = 0.5*t*h*(Lζ3*θTΦ21 + Lζ2*θTΦ31)
                            B2s22 = 0.5*t*h*(Lζ3*θTΦ22 + Lζ2*θTΦ32)
                            B2s23 = 0.5*t*h*(Lζ3*θTΦ23 + Lζ2*θTΦ33)

                            B3s11 = 0.5*t*(Lh3*θTΦ11 + Lh1*θTΦ31)
                            B3s12 = 0.5*t*(Lh3*θTΦ12 + Lh1*θTΦ32)
                            B3s13 = 0.5*t*(Lh3*θTΦ13 + Lh1*θTΦ33)
                            B3s21 = 0.5*t*(Lh3*θTΦ21 + Lh2*θTΦ31)
                            B3s22 = 0.5*t*(Lh3*θTΦ22 + Lh2*θTΦ32)
                            B3s23 = 0.5*t*(Lh3*θTΦ23 + Lh2*θTΦ33)

                            # populating Bs matrix:
                            Bs[1,(n-1)*pdim+1] = B1s11
                            Bs[1,(n-1)*pdim+2] = B1s12
                            Bs[1,(n-1)*pdim+3] = B1s13
                            Bs[2,(n-1)*pdim+1] = B1s21
                            Bs[2,(n-1)*pdim+2] = B1s22
                            Bs[2,(n-1)*pdim+3] = B1s23
                            Bs[1,(n-1)*pdim+4] = B2s11
                            Bs[1,(n-1)*pdim+5] = B2s12
                            Bs[1,(n-1)*pdim+6] = B2s13
                            Bs[2,(n-1)*pdim+4] = B2s21
                            Bs[2,(n-1)*pdim+5] = B2s22
                            Bs[2,(n-1)*pdim+6] = B2s23  
                            Bs[3,(n-1)*pdim+4] = B3s11
                            Bs[3,(n-1)*pdim+5] = B3s12
                            Bs[3,(n-1)*pdim+6] = B3s13                          
                            Bs[4,(n-1)*pdim+4] = B3s21
                            Bs[4,(n-1)*pdim+5] = B3s22
                            Bs[4,(n-1)*pdim+6] = B3s23   
                        end
                        # updating stiffness matrix with the shear part at the gth Gauss point
                        mul!(DsBs, Ds, Bs)                 
                        mul!(BsTDsBs, transpose(Bs), DsBs)    
                        @. K1 += BsTDsBs * det * Ws[g]
                    end
                    
                    @inbounds for g in 1:gpet  # for each Gauss point of the torsional part

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        # -> tangent vectors at the gth point of eth element
                        xξ = XξGt[e,g]  # ∂x/∂ξ at the gth point of eth element (a1x)
                        yξ = YξGt[e,g]  # ∂y/∂ξ at the gth point of eth element (a1y)
                        zξ = ZξGt[e,g]  # ∂z/∂ξ at the gth point of eth element (a1z)
                        xη = XηGt[e,g]  # ∂x/∂η at the gth point of eth element (a2x)
                        yη = YηGt[e,g]  # ∂y/∂η at the gth point of eth element (a2y)
                        zη = ZηGt[e,g]  # ∂z/∂η at the gth point of eth element (a2z)

                        # generating local coordinate system at the gth point of the eth element:
                        e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)
                        xζ = e3x .*0.5 .* t # ∂x/∂η at the gth point of eth element
                        yζ = e3y .*0.5 .* t # ∂y/∂η at the gth point of eth element
                        zζ = e3z .*0.5 .* t # ∂z/∂η at the gth point of eth element

                        # Jacobi determinant and elements of inverse at the gth point of the eth element:
                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)
                                                    
                        @inbounds for n in 1:nonpe      # for each node of the eth element

                            # shape function and it's derivatives with respect to x,y,z:
                            h  = HGt[n,g] # shape function value at the nth node
                            hx = HξGt[n,g] * ξx + HηGt[n,g] * ηx + HζGt[n,g] * ζx # ∂h/∂x 
                            hy = HξGt[n,g] * ξy + HηGt[n,g] * ηy + HζGt[n,g] * ζy # ∂h/∂y
                            hz = HξGt[n,g] * ξz + HηGt[n,g] * ηz + HζGt[n,g] * ζz # ∂h/∂z                      

                            # gradient operators:
                            Lh1 = e1x*hx + e1y*hy + e1z*hz
                            Lh2 = e2x*hx + e2y*hy + e2z*hz

                            # calculating & populating Bt matrix:
                            Bt[(n-1)*pdim+1] = 0.5 * (Lh2*e1x - Lh1*e2x)
                            Bt[(n-1)*pdim+2] = 0.5 * (Lh2*e1y - Lh1*e2y)
                            Bt[(n-1)*pdim+3] = 0.5 * (Lh2*e1z - Lh1*e2z)
                            Bt[(n-1)*pdim+4] = h*e3x
                            Bt[(n-1)*pdim+5] = h*e3y
                            Bt[(n-1)*pdim+6] = h*e3z  
                        end
                        # updating stiffness matrix with the torsional part at the gth Gauss point
                        mul!(BtTBt, transpose(Bt), Bt)    
                        @. K1 += kt * G * t * BtTBt * det * Wt[g]
                    end

                    # filling up sparse matrix creator vectors:
                    @inbounds for k1 in 1:noqpe # for each dof per element
                        @inbounds for k2 in 1:noqpe # for each dof per element
                            ridx1[i] = nn2[e, k1]
                            cidx1[i] = nn2[e, k2]
                            data1[i] = K1[k1,k2]
                            i += 1
                        end
                    end 
                end
                append!(ridx,ridx1)
                append!(cidx,cidx1)
                append!(data,data1)
            end
        elseif secData isa Beam3D  # 3D Beam section
            # not yet implemented
        elseif secData isa Beam2D  # 2D Beam section  
            # not yet implemented
        elseif secData isa Beam1D  # 1D Beam section  
            # not yet implemented             
        elseif secData isa Truss3D # 3D Truss section
            # not yet implemented
        elseif secData isa Truss2D # 2D Truss section
            # not yet implemented
        elseif secData isa Truss1D # 1D Truss section
            # not yet implemented                    
        end
    end

    K = sparse(ridx, cidx, data, dof, dof)  # creating the global stiffness matrix
    dropzeros!(K)                           # drop zeros from stiffness matrix
    return 0.5*(K+K')   
end

function massMatrix(mesh,sec,DoFs;lumped=false)

    xyz = mesh.xyz           # coordinate matrix
    x   = xyz[:,1]           # x coordinates
    y   = xyz[:,2]           # y coordinates
    z   = xyz[:,3]           # z coordinates

    dof = maximum(DoFs)      # number of degrees of freedom of the model

    ridx = Vector{Int64}()   # initializing row index vector
    cidx = Vector{Int64}()   # initializing col index vector
    data = Vector{Float64}() # initializing data vector

    for (secKey, secData) in pairs(sec)  # for each section

        scope   = secData.scope      # scope of the current section
        ρ       = secData.mat.ρ      # mass density of the section   
        dim     = secData.dim        # spatial dimension of the current section
        pdim    = secData.pdim       # physical dimension of the current section
        intid   = secData.intid      # integration rule identifier of the current section

        setData = scope2selection(mesh,scope)  # convert scope to selection structure          

        if secData isa Solid3D                   # 3D Solid section
            for (eType, eTags) in setData.elems  # for each element type
                nnET    = mesh.elem[eType].nn    # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags  # element tags of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix of current section with curent element type
                nonpe   = mesh.elem[eType].non   # number of nodes per elements
                noe     = length(eTags);         # number of elements
                noqpe   = nonpe * pdim           # number of dofs per elements
                nov     = noqpe^2 * noe          # number of all individual stiffness values

                nn2 = zeros(Int32, noe, nonpe*pdim)   # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = DoFs[nn,i]
                end

                # properties of current element type:
                # rows: shape functions, colums: Gauss points
                HG  = mesh.elem[eType].HG[intid]  # ∂h/∂ξ derivatives at Gauss points  
                HξG = mesh.elem[eType].HξG[intid]  # ∂h/∂ξ derivatives at Gauss points   
                HηG = mesh.elem[eType].HηG[intid]  # ∂h/∂η derivatives at Gauss points   
                HζG = mesh.elem[eType].HζG[intid]  # ∂h/∂ζ derivatives at Gauss points   
                W   = mesh.elem[eType].GW[intid]   # Gauss weights

                gpe = length(W)             # number of Gauss points per elements

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix
                Z = z[nn]  # elemental z coordinate matrix

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> Gauss points
                Xξ = X * HξG # ∂x/∂ξ matrix
                Yξ = Y * HξG # ∂y/∂ξ matrix
                Zξ = Z * HξG # ∂z/∂ξ matrix
                Xη = X * HηG # ∂x/∂η matrix
                Yη = Y * HηG # ∂y/∂η matrix
                Zη = Z * HηG # ∂z/∂η matrix
                Xζ = X * HζG # ∂x/∂ζ matrix
                Yζ = Y * HζG # ∂y/∂ζ matrix
                Zζ = Z * HζG # ∂z/∂ζ matrix
                
                ridx1   = zeros(Int32, nov, 1)   # initialization of row index vector
                cidx1   = zeros(Int32, nov, 1)   # initialization of column index vector
                data1   = zeros(Float64, nov, 1) # initialization of data vector
                A       = zeros(Float64, 3, nonpe * pdim) # initialization of B matrix
                ATρA    = zeros(Float64, pdim*nonpe, pdim*nonpe)              
                M1      = zeros(Float64, pdim*nonpe, pdim*nonpe)  # initialization of stiffness matrix for eth element

                i = 1
                for e in 1:noe                   # for each element   
                    fill!(M1, 0.0)  # resetting stiffness matrix for the eth element
                    for g in 1:gpe               # for each Gauss point

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        xξ = Xξ[e,g]  # ∂x/∂ξ at the gth point of eth element
                        yξ = Yξ[e,g]  # ∂y/∂ξ at the gth point of eth element
                        zξ = Zξ[e,g]  # ∂z/∂ξ at the gth point of eth element
                        xη = Xη[e,g]  # ∂x/∂η at the gth point of eth element
                        yη = Yη[e,g]  # ∂y/∂η at the gth point of eth element
                        zη = Zη[e,g]  # ∂z/∂η at the gth point of eth element
                        xζ = Xζ[e,g]  # ∂x/∂ζ at the gth point of eth element
                        yζ = Yζ[e,g]  # ∂x/∂ζ at the gth point of eth element
                        zζ = Zζ[e,g]  # ∂x/∂ζ at the gth point of eth element

                        # Jacobi matrix determinant at the gth point of the eth element
                        det = xξ * (yη * zζ - yζ * zη) -
                              yξ * (xη * zζ - xζ * zη) +
                              zξ * (xη * yζ - xζ * yη)

                        for k = 1:nonpe
                            A[1, 3*k-2] = A[2, 3*k-1] = A[3, 3*k] = HG[k,g];
                        end

                        # adding the mass part form the gth point to the eth stiffness matrix:
                        mul!(ATρA, transpose(A), ρ*A)     # temp2 = B' * (D * B)
                        @. M1 += ATρA * det * W[g]        # final update
                   end

                    # filling up sparse matrix creator vectors:
                    for k1 in 1:noqpe      # for each dof per element
                        for k2 in 1:noqpe  # for each dof per element
                            ridx1[i] = nn2[e, k1]
                            cidx1[i] = nn2[e, k2]
                            data1[i] = M1[k1,k2]
                            i += 1
                        end
                    end 
                end
                append!(ridx,ridx1) # adding row indexer vector part from the eth element
                append!(cidx,cidx1) # adding column indexer vector part from the eth element
                append!(data,data1) # adding data vector part from the eth element
            end
        elseif secData isa Solid2D              # 2D Solid section
            secType = secData.type              # type of the 2D Solid section
            for (eType, eTags) in setData.elems # for each element type
                nnET    = mesh.elem[eType].nn   # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags # coonectivity matrix of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix
                nonpe   = mesh.elem[eType].non  # number of nodes per elements
                noe     = length(eTags)         # number of nodes
                noqpe   = nonpe * pdim          # number of dofs per elements
                nov     = noqpe^2 * noe         # number of all individual stiffness values

                nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = DoFs[nn,i]
                end
                
                # properties of current element type:
                # rows: shape functions, colums: Gauss points
                HG  = mesh.elem[eType].HG[intid]  # shape function values at Gauss points
                HξG = mesh.elem[eType].HξG[intid] # ∂h/∂ξ derivatives at Gauss points   
                HηG = mesh.elem[eType].HηG[intid] # ∂h/∂η derivatives at Gauss points   
                W   = mesh.elem[eType].GW[intid]  # Gauss weigths  
                
                gpe = length(W) # number of Gauss points per elements

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix             

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> Gauss points
                Xξ = X * HξG    # ∂x/∂ξ matrix
                Yξ = Y * HξG    # ∂y/∂ξ matrix
                Xη = X * HηG    # ∂x/∂η matrix
                Yη = Y * HηG    # ∂y/∂η matrix

                ridx1   = zeros(Int32, nov, 1)    # initialization of row index vector
                cidx1   = zeros(Int32, nov, 1)    # initialization of column index vector
                data1   = zeros(Float64, nov, 1)  # initialization of data vector 
                A       = zeros(Float64, 2, nonpe * pdim) # initialization of B matrix
                ATρA    = zeros(Float64, pdim*nonpe, pdim*nonpe)              
                M1      = zeros(Float64, pdim*nonpe, pdim*nonpe)  # initialization of stiffness matrix for eth element
                
                # initialization of B matrices and integration constants based on mechanical model:
                if secType == :PlaneStrain  # planestrain model
                    typeid = 1 # assigning typeid = 1 to this mechanical model
                    C = 1
                elseif secType == :PlaneStress # planestress model
                    typeid = 2 # assigning typeid = 2 to this mechanical model
                    C = secData.width  # width of the current section 
                elseif secType == :AxiSymmetricY  # axisymmetric model about axis y
                    typeid = 3  # assigning typeid = 3 to this mechanical model
                    R = X * HG  # radial coordinate matrix at the Gauss points
                elseif secType == :AxiSymmetricX  # axisymmetric model about axis x
                    typeid = 4  # assigning typeid = 4 to this mechanical model
                    R = Y * HG  # radial coordinate matrix at the Gauss points    
                end
                
                i = 1
                for e in 1:noe  # for each element
                    fill!(M1, 0.0)  # resetting stiffness matrix for the eth element  
                    for g in 1:gpe     # for each Gauss point

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        xξ = Xξ[e,g]  # ∂x/∂ξ at the gth point of eth element
                        yξ = Yξ[e,g]  # ∂y/∂ξ at the gth point of eth element
                        xη = Xη[e,g]  # ∂x/∂η at the gth point of eth element
                        yη = Yη[e,g]  # ∂y/∂η at the gth point of eth element

                        # Jacobi matrix determinant at the gth point of the eth element
                        det = xξ * yη - xη * yξ   
                        
                        for k = 1:nonpe
                            A[1, 2*k-1] = A[2, 2*k] = HG[k,g];
                        end   
                        
                        if typeid >= 3
                            r = R[e,g]
                            C = 2*r*π
                        end

                        # adding the mass part form the gth point to the eth stiffness matrix:
                        mul!(ATρA, transpose(A), ρ*A)     # temp2 = B' * (D * B)
                        @. M1 += ATρA * C * det * W[g]        # final update
                   end   
                   
                    # filling up sparse matrix creator vectors:
                    for k1 in 1:noqpe      # for each dof per element
                        for k2 in 1:noqpe  # for each dof per element
                            ridx1[i] = nn2[e, k1]
                            cidx1[i] = nn2[e, k2]
                            data1[i] = M1[k1,k2]
                            i += 1
                        end
                    end 
                end
                append!(ridx,ridx1) # adding row indexer vector part from the eth element
                append!(cidx,cidx1) # adding column indexer vector part from the eth element
                append!(data,data1) # adding data vector part from the eth element
            end                   
        elseif secData isa Shell  # Shell section

            t     = secData.width  # width of the current shell section
            intid = secData.intid   
            
            for (eType, eTags) in setData.elems # for each element type of the current shell section
                nnET    = mesh.elem[eType].nn   # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags # coonectivity matrix of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix
                nonpe   = mesh.elem[eType].non  # number of nodes per elements
                noe     = length(eTags)         # number of elements
                noqpe   = nonpe * pdim          # number of dofs per elements
                nov     = noqpe^2 * noe         # number of all individual stiffness values  

                ridx1 = zeros(Int32, nov, 1)    # initialization of row index vector
                cidx1 = zeros(Int32, nov, 1)    # initialization of column index vector
                data1 = zeros(Float64, nov, 1)  # initialization of data vector
                A     = zeros(Float64, 3, nonpe * pdim) # initialization of B matrix
                ATA  = zeros(Float64, pdim*nonpe, pdim*nonpe)              
                M1    = zeros(Float64, pdim*nonpe, pdim*nonpe)  # initialization of stiffness matrix for eth element

                nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = DoFs[nn,i]
                end

                # elemental coordinate matrices:
                # rows -> elements (e1,e2,e3...), columns -> nodes of the element (n1,n2,n3...)
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix  
                Z = z[nn]  # elemental z coordinate matrix
                
                # -> for membrane part of the stiffness matrix ----------------------------------------------------------
                HG  = mesh.elem[eType].HG[intid]   # shape function values at Gauss points for membrane part   
                HξG = mesh.elem[eType].HξG[intid]  # ∂h/∂ξ derivatives at Gauss points for membrane part   
                HηG = mesh.elem[eType].HηG[intid]  # ∂h/∂η derivatives at Gauss points for membrane part    
                HζG = mesh.elem[eType].HζG[intid]  # ∂h/∂ζ derivatives at Gauss points for membrane part    
                W   = mesh.elem[eType].GW[intid]   # Gauss weights at Gauss points for membrane part
                gpe = length(W)                  # number of Gauss points per elements for membrane part

                # x,y,z derivatives with respect to ξ,η,ζ (Jacobi elements) at Gauss points for membrane part:
                # rows -> elements (e1,e2,e3...), columns -> Gauss points (G1,G2,G3...)
                XξG = X * HξG # ∂x/∂ξ matrix
                YξG = Y * HξG # ∂y/∂ξ matrix
                ZξG = Z * HξG # ∂z/∂ξ matrix
                XηG = X * HηG # ∂x/∂η matrix
                YηG = Y * HηG # ∂y/∂η matrix
                ZηG = Z * HηG # ∂z/∂η matrix  

                HN  = mesh.elem[eType].HN
                HξN = mesh.elem[eType].HξN  # ∂h/∂ξ derivatives at Gauss points for membrane part   
                HηN = mesh.elem[eType].HηN  # ∂h/∂η derivatives at Gauss points for membrane part      

                # x,y,z derivatives with respect to ξ,η,ζ (Jacobi elements) at Gauss points for membrane part:
                # rows -> elements (e1,e2,e3...), columns -> Gauss points (G1,G2,G3...)
                XξN = X * HξN # ∂x/∂ξ matrix
                YξN = Y * HξN # ∂y/∂ξ matrix
                ZξN = Z * HξN # ∂z/∂ξ matrix
                XηN = X * HηN # ∂x/∂η matrix
                YηN = Y * HηN # ∂y/∂η matrix
                ZηN = Z * HηN # ∂z/∂η matrix  

                i = 1
                for e in 1:noe  # for each element of the current element type from current section
                    fill!(M1, 0.0)  # resetting stiffness matrix for the eth element 
                    if lumped 
                        m1 = zeros(nonpe)
                        g1 = 0  
                        g2 = 0 
                    end
                    for g in 1:gpe  # for each Gauss point of the membrane part

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        # -> tangent vectors at the gth point of eth element
                        xξ = XξG[e,g]  # ∂x/∂ξ at the gth point of eth element (a1x)
                        yξ = YξG[e,g]  # ∂y/∂ξ at the gth point of eth element (a1y)
                        zξ = ZξG[e,g]  # ∂z/∂ξ at the gth point of eth element (a1z)
                        xη = XηG[e,g]  # ∂x/∂η at the gth point of eth element (a2x)
                        yη = YηG[e,g]  # ∂y/∂η at the gth point of eth element (a2y)
                        zη = ZηG[e,g]  # ∂z/∂η at the gth point of eth element (a2z)

                        # generating local coordinate system at the gth point of the eth element:
                        e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)
                        xζ = e3x .*0.5 .* t # ∂x/∂η at the gth point of eth element
                        yζ = e3y .*0.5 .* t # ∂y/∂η at the gth point of eth element
                        zζ = e3z .*0.5 .* t # ∂z/∂η at the gth point of eth element   
                        
                        # Jacobi matrix determinant at the gth point of the eth element
                        det = xξ * (yη * zζ - yζ * zη) -
                              yξ * (xη * zζ - xζ * zη) +
                              zξ * (xη * yζ - xζ * yη)
          
                        if !lumped
                            for k = 1:nonpe
                                h = HG[k,g]
                                A[1, 6*k-5] = A[2, 6*k-4] = A[3, 6*k-3] = sqrt(2) * h;
                                A[1, 6*k-1] =  h*t/2*e3z * (1/3)
                                A[1, 6*k  ] = -h*t/2*e3y * (1/3)
                                A[2, 6*k-2] = -h*t/2*e3z * (1/3)
                                A[2, 6*k  ] =  h*t/2*e3x * (1/3)
                                A[3, 6*k-2] =  h*t/2*e3y * (1/3)
                                A[3, 6*k-1] = -h*t/2*e3x * (1/3)
                            end

                            # adding the mass part form the gth point to the eth stiffness matrix:
                            mul!(ATA, transpose(A), A)     # temp2 = B' * (D * B)
                            @. M1 += ATA * ρ * det * W[g]        # final update
                        else
                            g2 += ρ * det * W[g] 
                            for k = 1:nonpe
                                h = HG[k,g]   
                                g1 += h^2 * ρ * det * W[g] 
                                m1[k] += h^2 *2 * ρ * det * W[g] 
                            end
                        end
                    end

                    if lumped
                        m1 = m1./g1*g2
                        for k = 1:nonpe
                            xξ = XξN[e,k]  # ∂x/∂ξ at the gth point of eth element (a1x)
                            yξ = YξN[e,k]  # ∂y/∂ξ at the gth point of eth element (a1y)
                            zξ = ZξN[e,k]  # ∂z/∂ξ at the gth point of eth element (a1z)
                            xη = XηN[e,k]  # ∂x/∂η at the gth point of eth element (a2x)
                            yη = YηN[e,k]  # ∂y/∂η at the gth point of eth element (a2y)
                            zη = ZηN[e,k]  # ∂z/∂η at the gth point of eth element (a2z)  
                            
                            e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)

                            I1 = (m1[k]/4)*t^2
                            p = 6*(k-1)
                            M1[p+1,p+1] = M1[p+2,p+2] = M1[p+3,p+3] = m1[k]
                            M1[p+4,p+4] = (e2x^2+e1x^2)*I1
                            M1[p+5,p+5] = (e2y^2+e1y^2)*I1
                            M1[p+6,p+6] = (e2z^2+e1z^2)*I1
                            M1[p+4,p+5] = M1[p+5,p+4] = (e2x*e2y+e1y*e1x)*I1
                            M1[p+4,p+6] = M1[p+6,p+4] = (e2x*e2z+e1z*e1x)*I1
                            M1[p+5,p+6] = M1[p+6,p+5] = (e2y*e2z+e1z*e1y)*I1
                        end
                    end
                           
                    # filling up sparse matrix creator vectors:
                    for k1 in 1:noqpe      # for each dof per element
                        for k2 in 1:noqpe  # for each dof per element
                            ridx1[i] = nn2[e, k1]
                            cidx1[i] = nn2[e, k2]
                            data1[i] = M1[k1,k2]
                            i += 1
                        end
                    end 
                end                   
                append!(ridx,ridx1) # adding row indexer vector part from the eth element
                append!(cidx,cidx1) # adding column indexer vector part from the eth element
                append!(data,data1) # adding data vector part from the eth element
            end
        elseif secData isa Beam3D  # 3D Beam section
            # not yet implemented
        elseif secData isa Beam2D  # 2D Beam section  
            # not yet implemented
        elseif secData isa Beam1D  # 1D Beam section  
            # not yet implemented             
        elseif secData isa Truss3D # 3D Truss section
            # not yet implemented
        elseif secData isa Truss2D # 2D Truss section
            # not yet implemented
        elseif secData isa Truss1D # 1D Truss section
            # not yet implemented                    
        end
    end

    M = sparse(ridx, cidx, data, dof, dof)  # creating the global stiffness matrix
    dropzeros!(M)                           # drop zeros from stiffness matrix
    return 0.5*(M+M')   
end

function dampingMatrix(K,M,bcs,dampdef=Damping();wmax=nothing,Q=nothing,wn=nothing)

    # merging boundary conditions
    if bcs isa Vector || bcs isa Dict
        bcs = mergeBoundaryConditions(bcs)
    end

    Aset = collect(1:size(K,1))      # Aset vector (vector of all dof ids)
    Nset = setdiff(Aset,[bcs.Xset;bcs.Kset]) # vector of free dofs
    Ndim = length(Nset)

    Knn = manual_sparse_submatrix(K, Nset)
    Mnn = manual_sparse_submatrix(M, Nset)

    type = dampdef.type

    if type == :none
        Cnn = nothing
    elseif type == :Raylight
        α = dampdef.alpha
        β = dampdef.beta
        ξ = dampdef.xi

        if β === nothing
            if wmax == nothing
                # wmax calculation!
            end
            β = [2ξ[i]/(wmax)^(2i-1) for i in eachindex(ξ)]
        else
            β = to_vector(β)
        end

        Cnn = α * Mnn + β[1] * Knn
        if length(β) > 1
            if Ndim == nnz(Mnn)
                invMnn = spdiagm(1 ./ diag(Mnn))
            else
                invMnn = inv(Mnn)
            end
            MK = copy(Knn)
            iMK = invMnn * Knn
            for i in 2:length(β)
                MK *= iMK
                Cnn += β[i] * MK
            end
        end
        dropzeros!(Cnn)
    elseif type == :modal
        zeta = dampdef.zeta
        maxmode = dampdef.maxmode

        if maxmode == :last
            n = Ndim
        else
            n = max(maxmode,Ndim)
        end

        if Q === nothing || wn === nothing
            fmin = 1.01
            mi = 1000000
            ωmin² = (2π * fmin)^2 # converting the minimum frequency (Hz) into the square of the angular frequency
            ωn², Φraw = eigs(Knn, Mnn, nev=n, which=:LR, sigma=ωmin², maxiter=mi)
            ωn = sqrt.(abs.(real(ωn²))) # converting eigenvalues to angular frequencies
            Φ = real(Φraw) # taking real part of mode shapes
            norms = sqrt.(sum((Φ' * Mnn) .* Φ', dims=2))
            Q = Φ ./ norms'
        end

        if zeta isa Function
            ζ = zeta(ωn)
        else
            ζ = to_vector(zeta)
        end

        σ = ζ .* ωn
        Cm = spdiagm(0 => 2 .* σ )
       #  Cnn = sparse(Q' \ (Cm * (Q \ I)))
        Cnn = sparse(Q * Cm * Q')
        dropzeros!(Cnn)
    end
    return Cnn
end

function lingen(P, dura, dt)
    x = unique([0:dt:dura, dura])

    if P[1, 1] > 0
        P = vcat([0 0], P)
    end

    if P[end, 1] < dura
        P = vcat(P, [dura P[end, 2]])
    end

    n = length(x)
    N = size(P, 1)

    # Create a NearestNeighbor search tree
    tree = NearestNeighbors.KDTree(P[:, 1])  # Search on the first column of P
    idx = knn(tree, x, 1).indices  # Find the index of the closest value for each x
    
    X = x[idx]
    f = zeros(n)

    for i in 1:N-1
        xA = X[i]
        xB = X[i+1]
        yA = P[i, 2]
        yB = P[i+1, 2]
        rng = (x .>= xA) .& (x .<= xB)
        xC = x[rng]
        f[rng] .= (xC .- xA) ./ (xB .- xA) .* (yB .- yA) .+ yA
    end

    return f
end

function updateStiffness!(K,sys)

    # merging boundary conditions
    if sys isa Vector || sys isa Dict
        sys = mergeSystemMatrices(sys)
    end

    if !isempty(sys.Ke) # if the stiffness matrix from elastic support is not empty
        K +=  sys.Ke  # adding stiffness from elastic support
    end

    if !isempty(sys.Ks) # if the stiffness matrix from springs is not empty
        K += sys.Ks  # adding stiffness from springs
    end

    if !isempty(sys.G) # if the mpc matrix is not empty
        Gdim = size(sys.G,1) # number of MPC
        K = [K sys.G';sys.G zeros(Gdim,Gdim)]  # extending stiffness matrix with mpc matrix
    end

end

function updateStiffness(K,sys)

    # merging boundary conditions
    if sys isa Vector || sys isa Dict
        sys = mergeSystemMatrices(sys)
    end

    if !isempty(sys.Ke) # if the stiffness matrix from elastic support is not empty
        K += sys.Ke  # adding stiffness from elastic support
    end

    if !isempty(sys.Ks) # if the stiffness matrix from springs is not empty
        K += sys.Ks  # adding stiffness from springs
    end

    if !isempty(sys.G) # if the mpc matrix is not empty
        Gdim = size(sys.G,1) # number of MPCs
        K = [K sys.G';sys.G zeros(Gdim,Gdim)]  # extending stiffness matrix with mpc matrix
    end

    return K

end

# A function solving Kq = f directly (with f arrenged in a matrix) -------------------------------------------------
function statSolDirect(K,bcs;dt=:none,tmax=:none,deco=:none)

    # merging boundary conditions
    if bcs isa Vector || bcs isa Dict
        bcs = mergeBoundaryConditions(bcs)
    end

    Aset = collect(1:size(K,1))      # Aset vector (vector of all dof ids)
    Adim = length(Aset)              # number of all dofs
    Nset = setdiff(Aset,[bcs.Xset;bcs.Kset]) # vector of free dofs

    Knn = manual_sparse_submatrix(K, Nset)
    #  Knn = K[Nset,Nset]               # condensed stiffness matrix

    if deco == :Chol     # if Cholesky decomposition is required
        Knn = cholesky(Knn, check=false)  # decompose stiffness matrix
    end

    s1 = size(bcs.fconst,1);
    s2 = size(bcs.qconst,1);
    s3 = size(bcs.fbaseA,1);
    s4 = size(bcs.qbaseA,1);
    s5 = size(bcs.qbaseB,1);
    s6 = size(bcs.qbaseB,1);

    if dt != :none && tmax != :none  # if time domain is defined
        t = unique([collect(0:dt:tmax);tmax])  # creating time domain vector
        tNum = length(t)             # number of time steps
        f = zeros(Adim,tNum)         # initialization of nodal load matrix
        q = zeros(Adim,tNum)         # initialization of nodal displacenet matrix
        if s1 > 0     # if constant nodal load vector is not empty
             f[1:s1,:] .+= bcs.fconst        # inserting constant nodal load vector to all columns of nodal load matrix
        end
        if s2 > 0      # if constand nodal displacement vector is not empty 
            f[1:s2,:] .+= bcs.qconst    # inserting constant nodal displacement vector to all columns of nodal displacement matrix
        end
        if s3 > 0     # if base load matrix for variation with function is not empty
            fvarA = hcat([f.(t') for f in bcs.fvarA]...) # inserting the time doamin to the functions
            f[1:s3,:] .+= bcs.fbaseA * fvarA  # adding the variatonal part with functions to the global load matrix
        end
        if s4 > 0     # if base displacement matrix for variation with function is not empty
            qvarA = hcat(bcs.qvarA.(t)...)  # inserting the time doamin to the functions
            q[1:s4,:] .+= bcs.qbaseA * qvarA; # adding the variatonal part with functions to the global displacement matrix
        end        
        if s5 > 0      # if base load matrix for numeric variation is not empty
            f[1:s5,:] .+= bcs.fbaseB * bcs.fvarB # adding the variatonal part to the global load matrix
        end
        if s6 > 0      # if base displacement matrix for numeric variation is not empty
            q[1:s6,:] .+= bcs.qbaseB * bcs.qvarB # adding the variatonal part to the global displacement matrix
        end
    else
        t = [0.0]  # if time domain is not defined
        f = zeros(Adim,1) # initializing global nodal load vector
        q = zeros(Adim,1) # initializing global nodal displacement vector
        if s1 > 0 # if constand nodal load vector is not empty 
            f[1:s1] .+= bcs.fconst # adding constant global load vector to global nodal load vector
        end
        if s2 > 0 # if constand nodal displacement vector is not empty 
            q[1:s2] .+= bcs.qconst  # adding constant global displacement vector to global nodal displacement vector
        end
    end

    fn = f[Nset,:]  # condensed nodal load vector (matrix)

    qn = Knn\fn     # calculatng condensed nodal displacement vector (matrix)

    q[Nset,:] = qn  # insert condensed nodal displacement vector (matrix) to global nodal displacement vector (matrix)
   # v = Matrix{Float64}(undef,0,0)
   # a = Matrix{Float64}(undef,0,0)
  #  ls = [q[:,end] zeros(Adim,2)]  # last state (displacement, velocity, acceleration)

  #  return Solution(t,q,v,a,ls)
  return q, t
end

# A function calculating modal properties:
function modSol(K,M,bcs;n=6,fmin=1.01,mi=1000000,norm=:mass)

    # merging boundary conditions
    if bcs isa Vector || bcs isa Dict
        bcs = mergeBoundaryConditions(bcs)
    end

    Aset = collect(1:size(K,1))      # Aset vector (vector of all dof ids)
    Adim = length(Aset)              # number of all dofs
    Nset = setdiff(Aset,[bcs.Xset;bcs.Kset]) # vector of free dofs

    Knn = manual_sparse_submatrix(K, Nset)
    Mnn = manual_sparse_submatrix(M, Nset)

    ωmin² = (2π * fmin)^2 # converting the minimum frequency (Hz) into the square of the angular frequency
    ω², Φraw = eigs(Knn, Mnn, nev=n, which=:LR, sigma=ωmin², maxiter=mi)
    ω = sqrt.(abs.(real(ω²))) # converting eigenvalues to angular frequencies
    f = ω / 2π  # converting angular frequencies to frequencies in Hz
    Φ = real(Φraw) # taking real part of mode shapes

    # mode shape mormalization:
    if norm == :mass
        norms = sqrt.(sum((Φ' * Mnn) .* Φ', dims=2))
        Ψ = Φ ./ norms'
    elseif norm == :unit
        norms = sqrt.(sum(abs2, Φ; dims=1))  # 2-norm of each column
        Ψ = Φ ./ norms
    end

    Q = zeros(Adim,n) # initializing global nodal displacement vector
    Q[Nset,:] = Ψ  # insert condensed nodal displacement vector (matrix) to global nodal displacement vector (matrix)

    return Q, f
end

# A function calculating stress field --------------------------------------------------------------------------
function postCalc(mesh,sec,q,dofTable,fields=[:Sn];zeta=:top)

    Se = ElementNodeData() # initialization of an empty elemental (unaveraged) stress field
    Sn = NodeData()   # initialization of an empty nodal (averaged) stress field
    Ee = ElementNodeData() # initialization of an empty elemental (unaveraged) stress field
    En = NodeData()   # initialization of an empty nodal (averaged) stress field
    Ue = ElementNodeData() # initialization of an empty elemental (unaveraged) stress field
    Un = NodeData()   # initialization of an empty nodal (averaged) stress field

    sec = to_vector(sec)
    fields = to_vector(fields)

    xyz = mesh.xyz  # coordinate matrix
    x   = xyz[:,1]  # x coordinates
    y   = xyz[:,2]  # y coordinates
    z   = xyz[:,3]  # z coordinates

    non   = length(x)         # number of nodes
    nTags = collect(1:non)    # node tags
    dof   = maximum(dofTable) # number of degrees of freedom of the model
    epn   = zeros(Int, non)   # initialization of element per node vector

    tNum = size(q,2)          # number of timesteps

    c = 6  # number of independent stress components
    if :Se ∈ fields  # if elemental stress field is required
        Se.data = Array{Float64}(undef, 6, tNum, 0)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
        Se.name = :Se
    end
    if :Sn ∈ fields  # if nodal stress field is required
        Sn.data = zeros(6, tNum, non)  # initialization of nodal stress 3Dmatrix (6 x tNum x non)
        Sn.name = :Sn
    end
    if :Ee ∈ fields  # if elemental stress field is required
        Ee.data = Array{Float64}(undef, 6, tNum, 0)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
        Ee.name = :Ee
    end
    if :En ∈ fields  # if nodal stress field is required
        En.data = zeros(6, tNum, non)  # initialization of nodal stress 3Dmatrix (6 x tNum x non)
        En.name = :En
    end
    if :Ue ∈ fields  # if elemental stress field is required
        Ue.data = Array{Float64}(undef, 1, tNum, 0)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
        Ue.name = :Ue
    end
    if :Un ∈ fields  # if nodal stress field is required
        Un.data = zeros(1, tNum, non)  # initialization of nodal stress 3Dmatrix (6 x tNum x non)
        Un.name = :Un
    end

    for (secKey, secData) in pairs(sec)  # for each section

        scope = secData.scope      # Section name of the current section
        D     = secData.D          # material matrix
        dim   = secData.dim        # spatial dimension of the current section
        pdim  = secData.pdim       # physical dimension of the current section

        #setData = set(mesh,secName)  # retrieve set data from section name   
        
        setData = scope2selection(mesh,scope)  # convert scope to selection structure
        setDim  = setData.dim        # dimension of the set
        
        if secData isa Solid3D # 3D Solid section

            for (eType, eTags) in setData.elems  # for each element type of the current section
                nnET    = mesh.elem[eType].nn    # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags  # element tags of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix of current section with current element type
                nonpe   = mesh.elem[eType].non   # number of nodes per elements
                noe     = length(eTags)          # number of elements
                noqpe   = nonpe * pdim           # number of dofs per elements
                nov     = noqpe^2 * noe          # number of individual tensorfield values for whole model

              #  if ftype == :elemental || ftype == :both # if elemetal stress field is required
              #      SeET = zeros(6, tNum, noe * nonpe) # initialization of 3D tensorfield matrix for current element type (6 x tNum x noe*nonpe)
              #  end     
                if :Ee ∈ fields  # if elemental stress field is required
                    EeET = zeros(6, tNum, noe * nonpe)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end
                if :Se ∈ fields  # if elemental stress field is required
                    SeET = zeros(6, tNum, noe * nonpe)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end                
                if :Ue ∈ fields  # if elemental stress field is required
                    UeET = zeros(1, tNum, noe * nonpe)   # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end                

                nn2 = zeros(Int32, noe, nonpe*pdim) # initialization of dof connectivity matrix of current section with current element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # extracting shape function derivatives at nodes:
                HξN = mesh.elem[eType].HξN # ∂h/∂ξ at the nodes
                HηN = mesh.elem[eType].HηN # ∂h/∂η at the nodes 
                HζN = mesh.elem[eType].HζN # ∂h/∂ζ at the nodes

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate
                Y = y[nn]  # elemental y coordinate
                Z = z[nn]  # elemental z coordinate

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> nodes
                Xξ = X * HξN  # ∂x/∂ξ matrix at the nodes
                Yξ = Y * HξN  # ∂y/∂ξ matrix at the nodes
                Zξ = Z * HξN  # ∂z/∂ξ matrix at the nodes
                Xη = X * HηN  # ∂x/∂η matrix at the nodes
                Yη = Y * HηN  # ∂y/∂η matrix at the nodes
                Zη = Z * HηN  # ∂z/∂η matrix at the nodes
                Xζ = X * HζN  # ∂x/∂ζ matrix at the nodes
                Yζ = Y * HζN  # ∂y/∂ζ matrix at the nodes
                Zζ = Z * HζN  # ∂z/∂ζ matrix at the nodes

                B = zeros(Float64, c, nonpe * pdim) # initialization of B matrix

                i = 1
                for e = 1:noe     # for each element
                    q1 = q[nn2[e,:],:]  # extracting nodal displacements associated with the eth element
                    for n = 1:nonpe   # for each node of the eth element

                        # extracting Jacobi matrix elements at the nth node of the eth element:
                        xξ = Xξ[e,n] # ∂x/∂ξ at the nth node of eth element
                        yξ = Yξ[e,n] # ∂y/∂ξ at the nth node of eth element
                        zξ = Zξ[e,n] # ∂z/∂ξ at the nth node of eth element
                        xη = Xη[e,n] # ∂x/∂η at the nth node of eth element
                        yη = Yη[e,n] # ∂y/∂η at the nth node of eth element
                        zη = Zη[e,n] # ∂z/∂η at the nth node of eth element
                        xζ = Xζ[e,n] # ∂x/∂ζ at the nth node of eth element
                        yζ = Yζ[e,n] # ∂y/∂ζ at the nth node of eth element
                        zζ = Zζ[e,n] # ∂z/∂ζ at the nth node of eth element

                        # Jacobi matrix determinant at the nth node of the eth element
                        det = xξ * (yη * zζ - yζ * zη) -
                              yξ * (xη * zζ - xζ * zη) +
                              zξ * (xη * yζ - xζ * yη)

                        # elements of the Jacobi inverse at the nth node of the eth element:
                        ξx = (yη * zζ - yζ * zη) / det # ∂ξ/∂x at the nth node of eth element  
                        ηx = (yζ * zξ - yξ * zζ) / det # ∂η/∂x at the nth node of eth element  
                        ζx = (yξ * zη - yη * zξ) / det # ∂ζ/∂x at the nth node of eth element  
                        ξy = (xζ * zη - xη * zζ) / det # ∂ξ/∂y at the nth node of eth element  
                        ηy = (xξ * zζ - xζ * zξ) / det # ∂η/∂y at the nth node of eth element  
                        ζy = (xη * zξ - xξ * zη) / det # ∂ζ/∂y at the nth node of eth element  
                        ξz = (xη * yζ - xζ * yη) / det # ∂ξ/∂z at the nth node of eth element  
                        ηz = (xζ * yξ - xξ * yζ) / det # ∂η/∂z at the nth node of eth element  
                        ζz = (xξ * yη - xη * yξ) / det # ∂ζ/∂z at the nth node of eth element  

                        # assembling B matrix at the nth node of the eth element:
                        for k = 1:nonpe
                            hx = HξN[k,n] * ξx + HηN[k,n] * ηx + HζN[k,n] * ζx  # ∂h/∂x
                            hy = HξN[k,n] * ξy + HηN[k,n] * ηy + HζN[k,n] * ζy  # ∂h/∂y
                            hz = HξN[k,n] * ξz + HηN[k,n] * ηz + HζN[k,n] * ζz  # ∂h/∂z
                            B[1,3*k-2] = B[4,3*k-1] = B[6,3*k]   = hx  # filling up B matrix
                            B[2,3*k-1] = B[4,3*k-2] = B[5,3*k]   = hy  # filling up B matrix
                            B[3,3*k]   = B[5,3*k-1] = B[6,3*k-2] = hz  # filling up B matrix
                        end
                        E1 = B * q1
                        S1 = D * E1  # stress matrix at the nth node of the eth element (6 x tNum)

                        if :Ee ∈ fields  # if elemental stress field is required
                            EeET[:,:,(e-1)*nonpe+n] = E1    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :En ∈ fields  # if nodal stress field is required
                            En.data[:,:,nn[e,n]] += E1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end 
                        if :Se ∈ fields  # if elemental stress field is required
                            SeET[:,:,(e-1)*nonpe+n] = S1    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :Sn ∈ fields  # if nodal stress field is required
                            Sn.data[:,:,nn[e,n]] += S1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end                         
                        if :Ue ∈ fields  # if elemental stress field is required
                            UeET[1,:,(e-1)*nonpe+n] = sum(0.5*S1.*E1,dims=1)   # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :Un ∈ fields  # if nodal stress field is required
                            Un.data[1,:,nn[e,n]] += sum(0.5*S1.*E1,dims=1) # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end  


                    #=    if ftype == :elemental || ftype == :both # if unavareged stress is required
                            SeET[:,:,(e-1)*nonpe+n] = S1    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if ftype == :nodal || ftype == :both  # if nodal stress is required
                            Snd.tBlock[:,:,nn[e,n]] += S1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end =#
                    end
                end

                if :Ee ∈ fields  # if elemental stress field is required
                    Ee.data = cat(Ee.data, EeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Ee.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ee.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end
                if :Se ∈ fields  # if elemental stress field is required
                    Se.data = cat(Se.data, SeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Se.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end
                if :Ue ∈ fields  # if elemental stress field is required
                    Ue.data = cat(Ue.data, UeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Ue.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ue.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end

                if :En ∈ fields || :Sn ∈ fields || :Un ∈ fields # if elemental stress field is required
                    for x in nn
                        epn[x] += 1  # updating element per node counter
                    end
                end

               # if ftype == :elemental || ftype == :both  # if elemental stress is required
               #     Se.tBlock = cat(Se.tBlock, SeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
               #                                                  # of the global 3D elemental stress matrix       
               #     append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
               #     append!(Se.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
               # end

               # if ftype == :nodal || ftype == :both  # if nodal stress is required
               #     for x in nn
               #         epn[x] += 1  # updating element per node counter
               #     end
               # end
            end

            # S.type = :Solid3D

        elseif secData isa Solid2D # 2D Solid section
            secType = secData.type       # type of the 2D Solid section
            ν = secData.mat.ν

            for (eType, eTags) in setData.elems  # for each element type of the current section
                nnET    = mesh.elem[eType].nn    # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags  # element tags of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix of current section with current element type
                nonpe   = mesh.elem[eType].non   # number of nodes per elements
                noe     = length(eTags)          # number of elements
                noqpe   = nonpe * pdim           # number of dofs per elements
                nov     = noqpe^2 * noe          # number of individual stress values for whole model 
  
                if :Ee ∈ fields  # if elemental stress field is required
                    EeET = zeros(6, tNum, noe * nonpe)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end
                if :Se ∈ fields  # if elemental stress field is required
                    SeET = zeros(6, tNum, noe * nonpe)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end                
                if :Ue ∈ fields  # if elemental stress field is required
                    UeET = zeros(1, tNum, noe * nonpe)   # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end                   
                

               # if ftype == :elemental || ftype == :both # if unavareged stress filed is required
               #     Suavget = zeros(c, tNum, noe * nonpe) # initialization of 3D stress matrix for current element type (6 x tNum x noe*nonpe)
               # end        

                nn2 = zeros(Int32, noe, nonpe*pdim) # initialization of dof connectivity matrix of current section with current element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # extracting shape functions and derivatives at nodes:
                HN  = mesh.elem[eType].HN  # shape function values at the nodes
                HξN = mesh.elem[eType].HξN # ∂h/∂ξ at the nodes
                HηN = mesh.elem[eType].HηN # ∂h/∂η at the nodes 
                
                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate
                Y = y[nn]  # elemental y coordinate

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> nodes
                Xξ = X * HξN  # ∂x/∂ξ matrix at the nodes
                Yξ = Y * HξN  # ∂y/∂ξ matrix at the nodes
                Xη = X * HηN  # ∂x/∂η matrix at the nodes
                Yη = Y * HηN  # ∂y/∂η matrix at the nodes

                # initialization of B matrices based on mechanical model:
                if secType == :PlaneStrain  # planestrain model
                    r = [1,2,4,3]
                    B = zeros(Float64, 3, nonpe * pdim) # initialization of B matrix
                    typeid = 1 # assigning typeid = 1 to this mechanical model
                elseif secType == :PlaneStress # planestress model
                    r = [1,2,4]
                    B = zeros(Float64, 3, nonpe * pdim) # initialization of B matrix
                    typeid = 2 # assigning typeid = 2 to this mechanical model
                elseif secType == :AxiSymmetricY  # axisymmetric model about axis y
                    r = [1,2,3,4]
                    B = zeros(Float64, 4, nonpe * pdim) # initialization of B matrix
                    R = X * HN  # radial coordinate matrix at the nodes
                    d = 1
                    typeid = 3  # assigning typeid = 3 to this mechanical model
                elseif secType == :AxiSymmetricX  # axisymmetric model about axis x
                    r = [1,2,3,4]
                    B = zeros(Float64, 4, nonpe * pdim) # initialization of B matrix
                    R = Y * HN  # radial coordinate matrix at the nodes
                    d = 0
                    typeid = 4  # assigning typeid = 4 to this mechanical model
                end

                i = 1
                for e = 1:noe     # for each element
                    q1 = q[nn2[e,:],:]  # extracting nodal displacements associated with the eth element    
                    for n = 1:nonpe   # for each node of the eth element

                        # extracting Jacobi matrix elements at the nth node of the eth element:
                        xξ = Xξ[e,n] # ∂x/∂ξ at the nth node of eth element
                        yξ = Yξ[e,n] # ∂y/∂ξ at the nth node of eth element
                        xη = Xη[e,n] # ∂x/∂η at the nth node of eth element
                        yη = Yη[e,n] # ∂y/∂η at the nth node of eth element

                        # Jacobi matrix determinant at the nth node of the eth element
                        det = xξ * yη - xη * yξ     
                        
                        # elements of the Jacobi inverse at the nth node of the eth element:
                        ξx =  yη / det # ∂ξ/∂x at the gth point of eth element
                        ηx = -yξ / det # ∂η/∂x at the gth point of eth element
                        ξy = -xη / det # ∂ξ/∂y at the gth point of eth element
                        ηy =  xξ / det # ∂η/∂y at the gth point of eth element 

                        # assembling B matrix at the nth node of the eth element:
                        if typeid <= 2  # if mechanical model is planestrain or planestress
                            for k = 1:nonpe  # for each node of the eth element
                                hx = HξN[k,n] * ξx + HηN[k,n] * ηx # ∂h/∂x 
                                hy = HξN[k,n] * ξy + HηN[k,n] * ηy # ∂h/∂y 
                                B[1,2*k-1] = B[3,2*k] = hx  # filling up B matrix with ∂h/∂x values
                                B[2,2*k] = B[3,2*k-1] = hy  # filling up B matrix with ∂h/∂y values
                            end
                        elseif typeid >=3  # if mechanical model is axisymmetric
                            r = R[e,n]     # radial coordinate of the nth node of the eth element
                            C = 2 * r * π  # integration constant (length of the corresponding circle)
                            for k = 1:nonpe # for each node of the eth element
                                hx = HξN[k,n] * ξx + HηN[k,n] * ηx # ∂h/∂x 
                                hy = HξN[k,n] * ξy + HηN[k,n] * ηy # ∂h/∂y                                  
                                B[1,2*k-1] = B[4,2*k] = hx # filling up B matrix with ∂h/∂x values
                                B[2,2*k] = B[4,2*k-1] = hy # filling up B matrix with ∂h/∂y values
                                B[3,2*k-d] = HG[k,n]/r     # filling up B matrix with h/r values
                            end
                        end

                        E1 = B * q1
                        S1 = D * E1  # stress matrix at the nth node of the eth element (6 x tNum)

                        if :Ue ∈ fields  # if elemental stress field is required
                            UeET[1,:,(e-1)*nonpe+n] = sum(0.5*S1.*E1,dims=1)    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :Un ∈ fields  # if nodal stress field is required
                            Un.data[1,:,nn[e,n]] += sum(0.5*S1.*E1,dims=1) # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end  

                        if secType == :PlaneStrain
                            S1 = [S1;ν*(S1[1:1,:]+S1[2:2,:])]
                        end
                        if secType == :PlaneStrain
                            E1 = [E1;-(ν/(1-ν))*(SE1[1:1,:]+E1[2:2,:])]
                        end

                        if :Ee ∈ fields  # if elemental stress field is required
                            EeET[r,:,(e-1)*nonpe+n] = E1    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :En ∈ fields  # if nodal stress field is required
                            En.data[r,:,nn[e,n]] += E1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end 
                        if :Se ∈ fields  # if elemental stress field is required
                            SeET[r,:,(e-1)*nonpe+n] = S1    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :Sn ∈ fields  # if nodal stress field is required
                            Sn.data[r,:,nn[e,n]] += S1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end                         

                       #= if secType == :PlaneStrain
                            S1 = [S1;ν*(S1[1:1,:]+S1[2:2,:])]
                        end

                        if ftype == :elemental || ftype == :both # if unavareged stress is required
                            Suavget[r,:,(e-1)*nonpe+n] = S1    # inserting S1 to the next layer of the 3D elemental stress matrix for current element type                    
                        end
                        if ftype == :nodal || ftype == :both  # if nodal stress is required
                            Sn.tBlock[r,:,nn[e,n]] += S1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end =#
                    end
                end

                if :Ee ∈ fields  # if elemental stress field is required
                    Ee.data = cat(Ee.data, EeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Ee.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ee.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end
                if :Se ∈ fields  # if elemental stress field is required
                    Se.data = cat(Se.data, SeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Se.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end
                if :Ue ∈ fields  # if elemental stress field is required
                    Ue.data = cat(Ue.data, UeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Ue.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ue.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end

                if :En ∈ fields || :Sn ∈ fields || :Un ∈ fields # if elemental stress field is required
                    for x in nn
                        epn[x] += 1  # updating element per node counter
                    end
                end
                
              #=  if ftype == :elemental || ftype == :both  # if elemental stress is required
                    Se.tBlock = cat(Se.tBlock, Suavget; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                        # of the global 3D elemental stress matrix       
                    append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Se.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end

                if ftype == :nodal || ftype == :both  # if nodal stress is required
                    for x in nn
                        epn[x] += 1  # updating element per node counter
                    end
                end =#
            end
        elseif secData isa Shell  # Shell section
            t    = secData.width  # width of the current shell section
            intm = secData.intm   # integration rule identifier for membrane+bending part
            ints = secData.ints   # integration rule identifier for shear part
            intt = secData.intt   # integration rule identifier for torsional part
            kt   = secData.kt     # torsional defomration factor
            G    = secData.mat.G  # shear modulus

            DD = zeros(6,6)
            DD[1:5,1:5] = D
            DD[6,6] = 1

            # preparing material matrices: 
            # Dm = zeros(6,6) # initialization of material matrix for membrane+bending part
            # Dm[1:3,1:3] = 2*D[1:3,1:3] 
            # Dm[4:6,4:6] = (2/3)*D[1:3,1:3] # material matrix for membrane+bending part

            # Ds = zeros(4,4) # initialization of material matrix for shear part
            # Ds[1:2,1:2] = 2*D[4:5,4:5]
            # Ds[3:4,3:4] = (2/3)*D[4:5,4:5] # material matrix for shear part

            ζ = zeta ==  :bot ? -1 : zeta == :mid ? 0 : zeta == :top ? 1 : zeta
            
            for (eType, eTags) in setData.elems # for each element type of the current shell section
                nnET    = mesh.elem[eType].nn   # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags # coonectivity matrix of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix
                nonpe   = mesh.elem[eType].non  # number of nodes per elements
                noe     = length(eTags)         # number of elements
                noqpe   = nonpe * pdim          # number of dofs per elements
                nov     = noqpe^2 * noe         # number of all individual stiffness values  

                if :Ee ∈ fields  # if elemental stress field is required
                    EeET = zeros(6, tNum, noe * nonpe)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end
                if :Se ∈ fields  # if elemental stress field is required
                    SeET = zeros(6, tNum, noe * nonpe)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end                
                if :Ue ∈ fields  # if elemental stress field is required
                    UeET = zeros(1, tNum, noe * nonpe)   # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
                end  

               # if ftype == :elemental || ftype == :both # if elemetal stress field is required
               #     SetfET = zeros(c, tNum, noe * nonpe) # initialization of 3D tensorfield matrix for current element type (6 x tNum x noe*nonpe)
               # end   

                nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # extracting shape function derivatives at nodes:
                HN  = mesh.elem[eType].HN  # shape function values at the nodesr
                HξN = mesh.elem[eType].HξN # ∂h/∂ξ at the nodes
                HηN = mesh.elem[eType].HηN # ∂h/∂η at the nodes 
                HζN = mesh.elem[eType].HζN # ∂h/∂ζ at the nodes

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate
                Y = y[nn]  # elemental y coordinate
                Z = z[nn]  # elemental z coordinate
                
                # x,y,z derivatives with respect to ξ,η,ζ (Jacobi elements) at Nodes:
                # rows -> elements (e1,e2,e3...), columns -> nodes of the elements
                XξN = X * HξN # ∂x/∂ξ matrix
                YξN = Y * HξN # ∂y/∂ξ matrix
                ZξN = Z * HξN # ∂z/∂ξ matrix
                XηN = X * HηN # ∂x/∂η matrix
                YηN = Y * HηN # ∂y/∂η matrix
                ZηN = Z * HηN # ∂z/∂η matrix 
                
                B = zeros(Float64, c, nonpe * pdim) # initialization of B matrix

                i = 1
                for e = 1:noe     # for each element
                    n1 = nn2[e,:]
                    q1 = q[n1,:]  # extracting nodal displacements associated with the eth element
                    for n = 1:nonpe   # for each node of the eth element

                        # extracting Jacobi matrix elements at the nth node of the eth element:
                        # -> tangent vectors at the gth point of eth element
                        xξ = XξN[e,n]  # ∂x/∂ξ at the nth node of eth element (a1x)
                        yξ = YξN[e,n]  # ∂y/∂ξ at the nth node of eth element (a1y)
                        zξ = ZξN[e,n]  # ∂z/∂ξ at the nth node of eth element (a1z)
                        xη = XηN[e,n]  # ∂x/∂η at the nth node of eth element (a2x)
                        yη = YηN[e,n]  # ∂y/∂η at the nth node of eth element (a2y)
                        zη = ZηN[e,n]  # ∂z/∂η at the nth node of eth element (a2z)

                        # generating local coordinate system at the nth node of the eth element:
                        e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)
                        xζ = e3x .*0.5 .* t # ∂x/∂η at the nth node of eth element
                        yζ = e3y .*0.5 .* t # ∂y/∂η at the nth node of eth element
                        zζ = e3z .*0.5 .* t # ∂z/∂η at the nth node of eth element  
                                                
                        # Jacobi determinant and elements of inverse at the gth point of the eth element:
                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)  
                        
                        for k in 1:nonpe      # for each node of the eth element

                            # shape function and it's derivatives with respect to x,y,z:
                            h  = HN[k,n] # shape function value at the nth node
                            hx = HξN[k,n] * ξx + HηN[k,n] * ηx + HζN[k,n] * ζx # ∂h/∂x 
                            hy = HξN[k,n] * ξy + HηN[k,n] * ηy + HζN[k,n] * ζy # ∂h/∂y
                            hz = HξN[k,n] * ξz + HηN[k,n] * ηz + HζN[k,n] * ζz # ∂h/∂z  

                            # gradient operators:
                            Lh1 = e1x*hx + e1y*hy + e1z*hz
                            Lh2 = e2x*hx + e2y*hy + e2z*hz
                            Lh3 = e3x*hx + e3y*hy + e3z*hz
                            Lζ1 = e1x*ζx + e1y*ζy + e1z*ζz
                            Lζ2 = e2x*ζx + e2y*ζy + e2z*ζz
                            Lζ3 = e3x*ζx + e3y*ζy + e3z*ζz

                            # preparing θᵀΦ product:
                            θTΦ11 = e1z*e3y - e1y*e3z
                            θTΦ12 = e1x*e3z - e1z*e3x
                            θTΦ13 = e1y*e3x - e1x*e3y
                            θTΦ21 = e2z*e3y - e2y*e3z
                            θTΦ22 = e2x*e3z - e2z*e3x
                            θTΦ23 = e2y*e3x - e2x*e3y
                            θTΦ31 = e3z*e3y - e3y*e3z
                            θTΦ32 = e3x*e3z - e3z*e3x
                            θTΦ33 = e3y*e3x - e3x*e3y

                            # elements of the Bm matrix associated with the strain contribution
                            # of the inplane displacements at the nth node
                            B1m11 = Lh1*e1x
                            B1m12 = Lh1*e1y
                            B1m13 = Lh1*e1z
                            B1m21 = Lh2*e2x
                            B1m22 = Lh2*e2y
                            B1m23 = Lh2*e2z
                            B1m31 = Lh2*e1x + Lh1*e2x
                            B1m32 = Lh2*e1y + Lh1*e2y
                            B1m33 = Lh2*e1z + Lh1*e2z  

                            B2m11 = 0.5*t*h*Lζ1*θTΦ11
                            B2m12 = 0.5*t*h*Lζ1*θTΦ12
                            B2m13 = 0.5*t*h*Lζ1*θTΦ13
                            B2m21 = 0.5*t*h*Lζ2*θTΦ21
                            B2m22 = 0.5*t*h*Lζ2*θTΦ22
                            B2m23 = 0.5*t*h*Lζ2*θTΦ23
                            B2m31 = 0.5*t*h*(Lζ2*θTΦ11 + Lζ1*θTΦ21)
                            B2m32 = 0.5*t*h*(Lζ2*θTΦ12 + Lζ1*θTΦ22)
                            B2m33 = 0.5*t*h*(Lζ2*θTΦ13 + Lζ1*θTΦ23)

                            # elements of the Bm matrix associated with the strain contribution
                            # of the rotations at the nth node (also includes curvature effect)
                           # B3m11 = 0.5*t*(B1m13*e3y - B1m12*e3z)
                           # B3m12 = 0.5*t*(B1m11*e3z - B1m13*e3x)
                           # B3m13 = 0.5*t*(B1m12*e3x - B1m11*e3y)
                           # B3m21 = 0.5*t*(B1m23*e3y - B1m22*e3z)
                           # B3m22 = 0.5*t*(B1m21*e3z - B1m23*e3x)
                           # B3m23 = 0.5*t*(B1m22*e3x - B1m21*e3y)
                           # B3m31 = 0.5*t*(B1m33*e3y - B1m32*e3z)
                           # B3m32 = 0.5*t*(B1m31*e3z - B1m33*e3x)
                           # B3m33 = 0.5*t*(B1m32*e3x - B1m31*e3y)   

                            B3m11 = 0.5*t*Lh1*θTΦ11
                            B3m12 = 0.5*t*Lh1*θTΦ12
                            B3m13 = 0.5*t*Lh1*θTΦ13
                            B3m21 = 0.5*t*Lh2*θTΦ21
                            B3m22 = 0.5*t*Lh2*θTΦ22
                            B3m23 = 0.5*t*Lh2*θTΦ23
                            B3m31 = 0.5*t*(Lh2*θTΦ11 + Lh1*θTΦ21)
                            B3m32 = 0.5*t*(Lh2*θTΦ12 + Lh1*θTΦ22)
                            B3m33 = 0.5*t*(Lh2*θTΦ13 + Lh1*θTΦ23)   

                            # calculating elements of the Bs matrix:
                            B1s11 = Lh3*e1x + Lh1*e3x
                            B1s12 = Lh3*e1y + Lh1*e3y
                            B1s13 = Lh3*e1z + Lh1*e3z
                            B1s21 = Lh3*e2x + Lh2*e3x
                            B1s22 = Lh3*e2y + Lh2*e3y
                            B1s23 = Lh3*e2z + Lh2*e3z

                            B2s11 = 0.5*t*h*(Lζ3*θTΦ11 + Lζ1*θTΦ31)
                            B2s12 = 0.5*t*h*(Lζ3*θTΦ12 + Lζ1*θTΦ32)
                            B2s13 = 0.5*t*h*(Lζ3*θTΦ13 + Lζ1*θTΦ33)
                            B2s21 = 0.5*t*h*(Lζ3*θTΦ21 + Lζ2*θTΦ31)
                            B2s22 = 0.5*t*h*(Lζ3*θTΦ22 + Lζ2*θTΦ32)
                            B2s23 = 0.5*t*h*(Lζ3*θTΦ23 + Lζ2*θTΦ33)

                            B3s11 = 0.5*t*(Lh3*θTΦ11 + Lh1*θTΦ31)
                            B3s12 = 0.5*t*(Lh3*θTΦ12 + Lh1*θTΦ32)
                            B3s13 = 0.5*t*(Lh3*θTΦ13 + Lh1*θTΦ33)
                            B3s21 = 0.5*t*(Lh3*θTΦ21 + Lh2*θTΦ31)
                            B3s22 = 0.5*t*(Lh3*θTΦ22 + Lh2*θTΦ32)
                            B3s23 = 0.5*t*(Lh3*θTΦ23 + Lh2*θTΦ33)

                            B[1,(k-1)*pdim+1] = B1m11
                            B[1,(k-1)*pdim+2] = B1m12
                            B[1,(k-1)*pdim+3] = B1m13
                            B[2,(k-1)*pdim+1] = B1m21
                            B[2,(k-1)*pdim+2] = B1m22
                            B[2,(k-1)*pdim+3] = B1m23
                            B[3,(k-1)*pdim+1] = B1m31
                            B[3,(k-1)*pdim+2] = B1m32
                            B[3,(k-1)*pdim+3] = B1m33

                         #   B[4,(k-1)*pdim+1] = B1s11
                         #   B[4,(k-1)*pdim+2] = B1s12
                         #   B[4,(k-1)*pdim+3] = B1s13
                         #   B[5,(k-1)*pdim+1] = B1s21
                         #   B[5,(k-1)*pdim+2] = B1s22
                         #   B[5,(k-1)*pdim+3] = B1s23

                            B[1,(k-1)*pdim+4] = B2m11 + ζ*B3m11
                            B[1,(k-1)*pdim+5] = B2m12 + ζ*B3m12
                            B[1,(k-1)*pdim+6] = B2m13 + ζ*B3m13
                            B[2,(k-1)*pdim+4] = B2m21 + ζ*B3m21
                            B[2,(k-1)*pdim+5] = B2m22 + ζ*B3m22
                            B[2,(k-1)*pdim+6] = B2m23 + ζ*B3m23  
                            B[3,(k-1)*pdim+4] = B2m31 + ζ*B3m31
                            B[3,(k-1)*pdim+5] = B2m32 + ζ*B3m32
                            B[3,(k-1)*pdim+6] = B2m33 + ζ*B3m33 

                          #  B[4,(k-1)*pdim+4] = B2s11 + ζ*B3s11
                          #  B[4,(k-1)*pdim+5] = B2s12 + ζ*B3s12
                          #  B[4,(k-1)*pdim+6] = B2s13 + ζ*B3s13
                          #  B[5,(k-1)*pdim+4] = B2s21 + ζ*B3s21
                          #  B[5,(k-1)*pdim+5] = B2s22 + ζ*B3s22
                          #  B[5,(k-1)*pdim+6] = B2s23 + ζ*B3s23  

                          #  B[6,(k-1)*pdim+1] = 0.5 * (Lh2*e1x - Lh1*e2x)
                          #  B[6,(k-1)*pdim+2] = 0.5 * (Lh2*e1y - Lh1*e2y)
                          #  B[6,(k-1)*pdim+3] = 0.5 * (Lh2*e1z - Lh1*e2z)
                          #  B[6,(k-1)*pdim+4] = h*e3x
                          #  B[6,(k-1)*pdim+5] = h*e3y
                          #  B[6,(k-1)*pdim+6] = h*e3z 
                        end
                      #  S1 = DD * B * q1  # stress matrix at the nth node of the eth element (6 x tNum)
                        
                        E1 = B * q1
                        S1 = DD * E1  # stress matrix at the nth node of the eth element (6 x tNum)

                        if :Ee ∈ fields  # if elemental stress field is required
                            EeET[[1,2,4,5,6,3],:,(e-1)*nonpe+n] = E1    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :En ∈ fields  # if nodal stress field is required
                            En.data[[1,2,4,5,6,3],:,nn[e,n]] += E1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end 
                        if :Se ∈ fields  # if elemental stress field is required
                            SeET[[1,2,4,5,6,3],:,(e-1)*nonpe+n] = S1    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :Sn ∈ fields  # if nodal stress field is required
                            Sn.data[[1,2,4,5,6,3],:,nn[e,n]] += S1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end                         
                        if :Ue ∈ fields  # if elemental stress field is required
                            UeET[1,:,(e-1)*nonpe+n] = sum(0.5*S1.*E1,dims=1)    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                        end
                        if :Un ∈ fields  # if nodal stress field is required
                            Un.data[1,:,nn[e,n]] += sum(0.5*S1.*E1,dims=1) # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                        end  

                    #   if ftype == :elemental || ftype == :both # if unavareged stress is required
                     #       SetfET[[1,2,4,5,6,3],:,(e-1)*nonpe+n] = S1    # inserting S1 to the next layer of the 3D unaveraged stress matrix for current element type                    
                     #   end
                     #   if ftype == :nodal || ftype == :both  # if nodal stress is required
                     #       Sn.tBlock[[1,2,4,5,6,3],:,nn[e,n]] += S1 # adding S1 to the corresponting layer (global node number of the nth node) of the 3D averaged stress matrix
                     #   end 
                    end
                end

                if :Ee ∈ fields  # if elemental stress field is required
                    Ee.data = cat(Ee.data, EeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Ee.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ee.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end
                if :Se ∈ fields  # if elemental stress field is required
                    Se.data = cat(Se.data, SeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Se.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end
                if :Ue ∈ fields  # if elemental stress field is required
                    Ue.data = cat(Ue.data, UeET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
                                                          # of the global 3D elemental stress matrix                           
                    append!(Ue.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ue.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
                end

                if :En ∈ fields || :Sn ∈ fields || :Un ∈ fields # if elemental stress field is required
                    for x in nn
                        epn[x] += 1  # updating element per node counter
                    end
                end


               # if ftype == :elemental || ftype == :both  # if elemental stress is required
               #     Se.tBlock = cat(Se.tBlock, SetfET; dims=3)  # placing the 3D unavareged stress matrix of current element type behind the last layer
               #                                                  # of the global 3D elemental stress matrix       
               #     append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
               #     append!(Se.eTags,eTags)  # adding element tags of current element type from current section to global element tags vector
               # end

               # if ftype == :nodal || ftype == :both  # if nodal stress is required
               #     for x in nn
               #         epn[x] += 1  # updating element per node counter
               #     end
               # end
            end 
        elseif secData isa Beam3D  # 3D Beam section
            # not yet implemented
        elseif secData isa Beam2D  # 2D Beam section  
            # not yet implemented
        elseif secData isa Beam1D  # 1D Beam section  
            # not yet implemented             
        elseif secData isa Truss3D # 3D Truss section
            # not yet implemented
        elseif secData isa Truss2D # 2D Truss section
            # not yet implemented
        elseif secData isa Truss1D # 1D Truss section
            # not yet implemented                    
        end
    end

    if :En ∈ fields || :Sn ∈ fields || :Un ∈ fields
        epn[epn .== 0] .= 1  # replacing epn = 0 values to 1 for division
        if :En ∈ fields
            En.data ./= reshape(epn, 1, 1, :)  # dividing stress values of 3D nodal stress matrix with number of nodes per elements
        end
        if :Sn ∈ fields
            Sn.data ./= reshape(epn, 1, 1, :)  # dividing stress values of 3D nodal stress matrix with number of nodes per elements
        end
        if :Un ∈ fields
            Un.data ./= reshape(epn, 1, 1, :)  # dividing stress values of 3D nodal stress matrix with number of nodes per elements
        end
    end


     #if ftype == :nodal || ftype == :both   # if nodal stress is required
     #   epn[epn .== 0] .= 1  # replacing epn = 0 values to 1 for division
     #   Sn.tBlock ./= reshape(epn, 1, 1, :)  # dividing stress values of 3D nodal stress matrix with number of nodes per elements
    #end
    return Sn,Se,En,Ee,Un,Ue
end

# A function calculating strain energy ---------------------------------------------------------------------------
function postGauss(mesh,sec,q,dofTable,fields;zeta=:top)

    Se = ElementNodeData() # initialization of an empty elemental (unaveraged) stress field
    Sn = NodeData()   # initialization of an empty nodal (averaged) stress field
    Ee = ElementNodeData() # initialization of an empty elemental (unaveraged) stress field
    En = NodeData()   # initialization of an empty nodal (averaged) stress field
    Ue = ElementData() # initialization of an empty elemental (unaveraged) stress field

    sec = to_vector(sec)
    fields = to_vector(fields)

    xyz = mesh.xyz           # coordinate matrix
    x   = xyz[:,1]           # x coordinates
    y   = xyz[:,2]           # y coordinates
    z   = xyz[:,3]           # z coordinates

    non = length(x)          # number of nodes
    dof = maximum(dofTable)  # number of degrees of freedom of the model

    tNum = size(q,2)  # number of timesteps
    epn   = zeros(Int, non)   # initialization of element per node vector

    c = 6  # number of independent stress components
    if :Se ∈ fields  # if elemental stress field is required
        Se.data = Array{Float64}(undef, 6, tNum, 0)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
        Se.name = :Sge
    end
    if :Sn ∈ fields  # if nodal stress field is required
        Sn.data = zeros(6, tNum, non)  # initialization of nodal stress 3Dmatrix (6 x tNum x non)
        Sn.name = :Sgn
    end
    if :Ee ∈ fields  # if elemental stress field is required
        Ee.data = Array{Float64}(undef, 6, tNum, 0)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
        Ee.name = :Ege
    end
    if :En ∈ fields  # if nodal stress field is required
        En.data = zeros(6, tNum, non)  # initialization of nodal stress 3Dmatrix (6 x tNum x non)
        En.name = :Egn
    end
    if :Ue ∈ fields  # if elemental stress field is required
        Ue.data = Array{Float64}(undef, 1, tNum, 0)  # initialization of unavareged stress 3Dmatrix (6 x tNum x 0)
        Ue.name = :Uge
    end

    for (secKey, secData) in pairs(sec)  # for each section

        scope   = secData.scope      # scope of the current section
        D       = secData.D          # material matrix
        dim     = secData.dim        # spatial dimension of the current section
        pdim    = secData.pdim       # physical dimension of the current section
        intid   = secData.intid      # integration rule identifier of the current section
        ν       = secData.mat.ν

        setData = scope2selection(mesh,scope)  # convert scope to selection structure      

        if secData isa Solid3D                   # 3D Solid section
            for (eType, eTags) in setData.elems  # for each element type
                nnET    = mesh.elem[eType].nn    # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags  # element tags of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix of current section with curent element type
                nonpe   = mesh.elem[eType].non   # number of nodes per elements
                noe     = length(eTags);         # number of elements
                noqpe   = nonpe * pdim           # number of dofs per elements
                nov     = noqpe^2 * noe          # number of all individual stiffness values

                nn2 = zeros(Int32, noe, nonpe*pdim)   # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end  

                # properties of current element type:
                # rows: shape functions, colums: Gauss points
                HG  = mesh.elem[eType].HG[intid]
                HξG = mesh.elem[eType].HξG[intid]  # ∂h/∂ξ derivatives at Gauss points   
                HηG = mesh.elem[eType].HηG[intid]  # ∂h/∂η derivatives at Gauss points   
                HζG = mesh.elem[eType].HζG[intid]  # ∂h/∂ζ derivatives at Gauss points   
                W   = mesh.elem[eType].GW[intid]   # Gauss weights

                iHG = pinv(HG)

                gpe = length(W)             # number of Gauss points per elements

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix
                Z = z[nn]  # elemental z coordinate matrix

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> Gauss points
                Xξ = X * HξG # ∂x/∂ξ matrix
                Yξ = Y * HξG # ∂y/∂ξ matrix
                Zξ = Z * HξG # ∂z/∂ξ matrix
                Xη = X * HηG # ∂x/∂η matrix
                Yη = Y * HηG # ∂y/∂η matrix
                Zη = Z * HηG # ∂z/∂η matrix
                Xζ = X * HζG # ∂x/∂ζ matrix
                Yζ = Y * HζG # ∂y/∂ζ matrix
                Zζ = Z * HζG # ∂z/∂ζ matrix

                B   = zeros(Float64, 6, nonpe * pdim) # initialization of B matrix 
                ϵ   = zeros(Float64, 6, tNum)
                σ   = zeros(Float64, 6, tNum)

                if :En ∈ fields  || :Ee ∈ fields           
                    E1g = zeros(Float64,6,tNum,gpe)
                    E1n = zeros(Float64,6,tNum,nonpe)
                end
                if :Sn ∈ fields  || :Se ∈ fields  
                    S1g = zeros(Float64,6,tNum,gpe)
                    S1n = zeros(Float64,6,tNum,nonpe)
                end
                if :Ue ∈ fields
                    UET = zeros(Float64,1,tNum,noe)
                    U1  = zeros(Float64,1,tNum)
                end

                i = 1
                for e in 1:noe                   # for each element
                    if :En ∈ fields  || :Ee ∈ fields 
                        fill!(E1g, 0.0)
                    end
                    if :Sn ∈ fields  || :Se ∈ fields  
                        fill!(S1g, 0.0)
                    end
                    if :Ue ∈ fields
                        fill!(U1, 0.0) 
                    end

                    n1 = nn2[e,:]
                    q1 = q[n1,:]  # extracting nodal displacements associated with the eth element
                                   
                    for g in 1:gpe               # for each Gauss point

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        xξ = Xξ[e,g]  # ∂x/∂ξ at the gth point of eth element
                        yξ = Yξ[e,g]  # ∂y/∂ξ at the gth point of eth element
                        zξ = Zξ[e,g]  # ∂z/∂ξ at the gth point of eth element
                        xη = Xη[e,g]  # ∂x/∂η at the gth point of eth element
                        yη = Yη[e,g]  # ∂y/∂η at the gth point of eth element
                        zη = Zη[e,g]  # ∂z/∂η at the gth point of eth element
                        xζ = Xζ[e,g]  # ∂x/∂ζ at the gth point of eth element
                        yζ = Yζ[e,g]  # ∂x/∂ζ at the gth point of eth element
                        zζ = Zζ[e,g]  # ∂x/∂ζ at the gth point of eth element

                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)

                        # assembling B matrix at the gth point of the eth element:
                        for k = 1:nonpe  # for each node of the eth element
                            hx = HξG[k,g] * ξx + HηG[k,g] * ηx + HζG[k,g] * ζx # ∂h/∂x 
                            hy = HξG[k,g] * ξy + HηG[k,g] * ηy + HζG[k,g] * ζy # ∂h/∂y
                            hz = HξG[k,g] * ξz + HηG[k,g] * ηz + HζG[k,g] * ζz # ∂h/∂z
                            B[1,3*k-2] = B[4,3*k-1] = B[6,3*k]   = hx  # filling up B matrix
                            B[2,3*k-1] = B[4,3*k-2] = B[5,3*k]   = hy  # filling up B matrix
                            B[3,3*k]   = B[5,3*k-1] = B[6,3*k-2] = hz  # filling up B matrix
                        end

                        mul!(ϵ, B, q1)              
                        mul!(σ, D, ϵ) 
                        if :En ∈ fields  || :Ee ∈ fields 
                            E1g[:,:,g] = ϵ
                        end
                        if :Sn ∈ fields  || :Se ∈ fields  
                            S1g[:,:,g] = σ
                        end
                        if :Ue ∈ fields
                            U1 += 0.5*sum(σ.*ϵ,dims=1) * det * W[g]  
                        end
                    end
                    if :Ue ∈ fields
                        UET[1,:,e] = U1
                    end

                    for t = 1:tNum
                        if :En ∈ fields  || :Ee ∈ fields 
                            E1n[:,t,:] = view(E1g,:,i,:) * iHG
                        end
                        if :Sn ∈ fields  || :Se ∈ fields  
                            S1n[:,t,:] = view(S1g,:,i,:) * iHG
                        end
                    end

                    if :Ee ∈ fields 
                        Ee.data = cat(Ee.data, E1n; dims=3)
                    end
                    if :Se ∈ fields 
                        Se.data = cat(Se.data, S1n; dims=3)  
                    end
                    if :En ∈ fields                                                              
                        En.data[:,:,nn[e,:]] += E1n
                    end
                    if :Sn ∈ fields  
                        Sn.data[:,:,nn[e,:]] += S1n   
                    end      
                end

                if :Ee ∈ fields 
                    append!(Ee.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ee.eTags,eTags)  # adding element     
                end
                if :Se ∈ fields         
                    append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Se.eTags,eTags)  # adding element 
                end
                if :Ue ∈ fields 
                    Ue.data = cat(Ue.data, UET; dims=3) 
                    append!(Ue.eTags,eTags)  # adding element 
                end

                if :En ∈ fields || :Sn ∈ fields || # if nodal field is required
                    for x in nn
                        epn[x] += 1  # updating element per node counter
                    end
                end
            end
        elseif secData isa Solid2D  # 2D Solid section
            secType = secData.type       # type of the 2D Solid section
            for (eType, eTags) in setData.elems # for each element type
                nnET    = mesh.elem[eType].nn   # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags # coonectivity matrix of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix
                nonpe   = mesh.elem[eType].non  # number of nodes per elements
                noe     = length(eTags)         # number of nodes
                noqpe   = nonpe * pdim          # number of dofs per elements
                nov     = noqpe^2 * noe         # number of all individual stiffness values

                nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # properties of current element type:
                # rows: shape functions, colums: Gauss points
                HG  = mesh.elem[eType].HG[intid]  # shape function values at Gauss points
                HξG = mesh.elem[eType].HξG[intid] # ∂h/∂ξ derivatives at Gauss points   
                HηG = mesh.elem[eType].HηG[intid] # ∂h/∂η derivatives at Gauss points   
                W   = mesh.elem[eType].GW[intid]  # Gauss weigths  

                iHG = pinv(HG)
                
                gpe = length(W) # number of Gauss points per elements

                # elemental coordinate matrices:
                # rows -> elements, colums -> element nodes
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix             

                # Jacobi matrix elements arranged in matrices:
                # rows -> elements, colums -> Gauss points
                Xξ = X * HξG    # ∂x/∂ξ matrix
                Yξ = Y * HξG    # ∂y/∂ξ matrix
                Xη = X * HηG    # ∂x/∂η matrix
                Yη = Y * HηG    # ∂y/∂η matrix
              
                # initialization of B matrices and integration constants based on mechanical model:
                if secType == :PlaneStrain  # planestrain model
                    er = 3                
                    sr = 4  
                    re = [1,2,4]              
                    rs = [1,2,4,3]
                    B = zeros(Float64, 3, nonpe * pdim) # initialization of B matrix
                    ϵ   = zeros(Float64, 3, tNum)
                    σ   = zeros(Float64, 3, tNum)
                    C = 1
                    typeid = 1 # assigning typeid = 1 to this mechanical model
                elseif secType == :PlaneStress # planestress model
                    er = 4
                    sr = 3
                    re = [1,2,4,3]
                    rs = [1,2,4]
                    B = zeros(Float64, 3, nonpe * pdim) # initialization of B matrix
                    ϵ   = zeros(Float64, 3, tNum)
                    σ   = zeros(Float64, 3, tNum)
                    C = secData.width  # width of the current section
                    typeid = 2 # assigning typeid = 2 to this mechanical model
                elseif secType == :AxiSymmetricY  # axisymmetric model about axis y
                    er = es = 4
                    re = rs = [1,2,3,4]
                    B = zeros(Float64, 4, nonpe * pdim) # initialization of B matrix
                    ϵ   = zeros(Float64, 4, tNum)
                    σ   = zeros(Float64, 4, tNum)
                    R = X * HG  # radial coordinate matrix at the Gauss points
                    d = 1
                    typeid = 3  # assigning typeid = 3 to this mechanical model
                elseif secType == :AxiSymmetricX  # axisymmetric model about axis x
                    er = es = 4
                    re = rs = [1,2,3,4]
                    B = zeros(Float64, 4, nonpe * pdim) # initialization of B matrix
                    ϵ   = zeros(Float64, 4, tNum)
                    σ   = zeros(Float64, 4, tNum)
                    R = Y * HG  # radial coordinate matrix at the Gauss points
                    d = 0
                    typeid = 4  # assigning typeid = 4 to this mechanical model
                end

                if :En ∈ fields  || :Ee ∈ fields           
                    E1g = zeros(Float64,er,tNum,gpe)
                    E1n = zeros(Float64,6,tNum,nonpe)
                end
                if :Sn ∈ fields  || :Se ∈ fields  
                    S1g = zeros(Float64,sr,tNum,gpe)
                    S1n = zeros(Float64,6,tNum,nonpe)
                end
                if :Ue ∈ fields
                    UET = zeros(Float64,1,tNum,noe)
                    U1  = zeros(Float64,1,tNum)
                end

                i = 1
                for e in 1:noe  # for each element
                    if :En ∈ fields  || :Ee ∈ fields 
                        fill!(E1g, 0.0)
                    end
                    if :Sn ∈ fields  || :Se ∈ fields  
                        fill!(S1g, 0.0)
                    end
                    if :Ue ∈ fields
                        fill!(U1, 0.0) 
                    end

                    n1 = nn2[e,:]
                    q1 = q[n1,:]  # extracting nodal displacements associated with the eth element

                    for g in 1:gpe     # for each Gauss point

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        xξ = Xξ[e,g]  # ∂x/∂ξ at the gth point of eth element
                        yξ = Yξ[e,g]  # ∂y/∂ξ at the gth point of eth element
                        xη = Xη[e,g]  # ∂x/∂η at the gth point of eth element
                        yη = Yη[e,g]  # ∂y/∂η at the gth point of eth element

                        # Jacobi matrix determinant at the gth point of the eth element
                        det = xξ * yη - xη * yξ

                        # elements of the Jacobi inverse at the gth point of the eth element:
                        ξx =  yη / det # ∂ξ/∂x at the gth point of eth element
                        ηx = -yξ / det # ∂η/∂x at the gth point of eth element
                        ξy = -xη / det # ∂ξ/∂y at the gth point of eth element
                        ηy =  xξ / det # ∂η/∂y at the gth point of eth element    

                        # assembling B matrix at the gth point of the eth element:
                        if typeid <= 2  # if mechanical model is planestrain or planestress
                            @inbounds for k = 1:nonpe  # for each node of the eth element
                                hx = HξG[k,g] * ξx + HηG[k,g] * ηx # ∂h/∂x 
                                hy = HξG[k,g] * ξy + HηG[k,g] * ηy # ∂h/∂y 
                                B[1,2*k-1] = B[3,2*k] = hx  # filling up B matrix with ∂h/∂x values
                                B[2,2*k] = B[3,2*k-1] = hy  # filling up B matrix with ∂h/∂y values
                            end
                        elseif typeid >=3  # if mechanical model is axisymmetric
                            r = R[e,g]     # radial coordinate of the gth point of the eth element
                            C = 2 * r * π  # integration constant (length of the corresponding circle)
                            @inbounds for k = 1:nonpe # for each node of the eth element
                                hx = HξG[k,g] * ξx + HηG[k,g] * ηx # ∂h/∂x 
                                hy = HξG[k,g] * ξy + HηG[k,g] * ηy # ∂h/∂y                                  
                                B[1,2*k-1] = B[4,2*k] = hx # filling up B matrix with ∂h/∂x values
                                B[2,2*k] = B[4,2*k-1] = hy # filling up B matrix with ∂h/∂y values
                                B[3,2*k-d] = HG[k,g]/r     # filling up B matrix with h/r values
                            end
                        end

                        # adding the stiffness part form the gth point to the eth stiffness matrix:

                        mul!(ϵ, B, q1)              
                        mul!(σ, D, ϵ) 
                        if :Ue ∈ fields
                            U1 += 0.5*sum(σ.*ϵ,dims=1) * det * W[g]  
                        end

                        if secType == :PlaneStrain
                            σe = [σ;ν*(σ[1:1,:]+σ[2:2,:])]
                        else
                            σe = σ
                        end
                        if secType == :PlaneStress
                            ϵe = [ϵ;-(ν/(1-ν))*(ϵ[1:1,:]+ϵ[2:2,:])]
                        else
                            ϵe = ϵ
                        end

                        if :En ∈ fields  || :Ee ∈ fields 
                            E1g[:,:,g] = ϵe
                        end
                        if :Sn ∈ fields  || :Se ∈ fields  
                            S1g[:,:,g] = σe
                        end                         
                    end

                    if :Ue ∈ fields
                        UET[1,:,e] = U1
                    end

                    for t = 1:tNum
                        if :En ∈ fields  || :Ee ∈ fields 
                            E1n[re,t,:] = view(E1g,:,t,:) * iHG
                        end
                        if :Sn ∈ fields  || :Se ∈ fields  
                            S1n[rs,t,:] = view(S1g,:,t,:) * iHG
                        end
                    end

                    if :Ee ∈ fields 
                        Ee.data = cat(Ee.data, E1n; dims=3)
                    end
                    if :Se ∈ fields 
                        Se.data = cat(Se.data, S1n; dims=3)  
                    end
                    if :En ∈ fields                                                              
                        En.data[:,:,nn[e,:]] += E1n
                    end
                    if :Sn ∈ fields  
                        Sn.data[:,:,nn[e,:]] += S1n   
                    end      
                end

                if :Ee ∈ fields 
                    append!(Ee.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ee.eTags,eTags)  # adding element     
                end
                if :Se ∈ fields         
                    append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Se.eTags,eTags)  # adding element 
                end
                if :Ue ∈ fields 
                    Ue.data = cat(Ue.data, UET; dims=3) 
                    append!(Ue.eTags,eTags)  # adding element 
                end

                if :En ∈ fields || :Sn ∈ fields # if nodal field is required
                    for x in nn
                        epn[x] += 1  # updating element per node counter
                    end

                end
            end
        elseif secData isa Shell  # Shell section

            t    = secData.width  # width of the current shell section
            intm = secData.intm   # integration rule identifier for membrane+bending part
            ints = secData.ints   # integration rule identifier for shear part
            intt = secData.intt   # integration rule identifier for torsional part
            kt   = secData.kt     # torsional defomration factor
            G    = secData.mat.G  # shear modulus

            # preparing material matrices: 

            # preparing material matrices: 
            Dm = zeros(6,6) # initialization of material matrix for membrane+bending part
            Dm[1:3,1:3] = 2*D[1:3,1:3] 
            Dm[4:6,4:6] = (2/3)*D[1:3,1:3] # material matrix for membrane+bending part

            Ds = zeros(4,4) # initialization of material matrix for shear part
            Ds[1:2,1:2] = 2*D[4:5,4:5]
            Ds[3:4,3:4] = (2/3)*D[4:5,4:5] # material matrix for shear part

            Dm1 = D[1:3,1:3] 
            Ds1 = D[4:5,4:5]

            ζ = zeta ==  :bot ? -1 : zeta == :mid ? 0 : zeta == :top ? 1 : zeta

            for (eType, eTags) in setData.elems # for each element type of the current shell section
                nnET    = mesh.elem[eType].nn   # coonectivity matrix of current element type
                eTagsET = mesh.elem[eType].tags # coonectivity matrix of current element type
                nn      = nnET[indexin(eTags, eTagsET), :]  # connectivity matrix
                nonpe   = mesh.elem[eType].non  # number of nodes per elements
                noe     = length(eTags)         # number of elements
                noqpe   = nonpe * pdim          # number of dofs per elements
                nov     = noqpe^2 * noe         # number of all individual stiffness values  

                nn2 = zeros(Int32, noe, nonpe*pdim)  # initialization of dof connectivity matrix of current section with curent element type
                for i = 1:pdim
                    nn2[:,i:pdim:pdim*nonpe-(pdim-i)] = dofTable[nn,i]
                end

                # elemental coordinate matrices:
                # rows -> elements (e1,e2,e3...), columns -> nodes of the element (n1,n2,n3...)
                X = x[nn]  # elemental x coordinate matrix
                Y = y[nn]  # elemental y coordinate matrix  
                Z = z[nn]  # elemental z coordinate matrix  
                
                if intm > 0
                    # -> for membrane part of the stiffness matrix ----------------------------------------------------------
                    HGm  = mesh.elem[eType].HG[intm]   # shape function values at Gauss points for membrane part   
                    HξGm = mesh.elem[eType].HξG[intm]  # ∂h/∂ξ derivatives at Gauss points for membrane part   
                    HηGm = mesh.elem[eType].HηG[intm]  # ∂h/∂η derivatives at Gauss points for membrane part    
                    HζGm = mesh.elem[eType].HζG[intm]  # ∂h/∂ζ derivatives at Gauss points for membrane part    
                    Wm   = mesh.elem[eType].GW[intm]   # Gauss weights at Gauss points for membrane part
                    gpem = length(Wm)                  # number of Gauss points per elements for membrane part
                    iHGm = pinv(HGm)

                    # x,y,z derivatives with respect to ξ,η,ζ (Jacobi elements) at Gauss points for membrane part:
                    # rows -> elements (e1,e2,e3...), columns -> Gauss points (G1,G2,G3...)
                    XξGm = X * HξGm # ∂x/∂ξ matrix
                    YξGm = Y * HξGm # ∂y/∂ξ matrix
                    ZξGm = Z * HξGm # ∂z/∂ξ matrix
                    XηGm = X * HηGm # ∂x/∂η matrix
                    YηGm = Y * HηGm # ∂y/∂η matrix
                    ZηGm = Z * HηGm # ∂z/∂η matrix    
                    
                    ϵm   = zeros(Float64, 6, tNum)
                    σm   = zeros(Float64, 6, tNum)
                    ϵmζ  = zeros(Float64, 3, tNum)
                    σmζ  = zeros(Float64, 3, tNum) 

                    if :En ∈ fields  || :Ee ∈ fields           
                        E1gm = zeros(Float64,3,tNum,gpem)
                    end
                    if :Sn ∈ fields  || :Se ∈ fields  
                        S1gm = zeros(Float64,3,tNum,gpem)
                    end

                else
                    gpem = 0         
                end 
                
                if ints>0
                    # -> for shear part of the stiffness matrix ----------------------------------------------------------
                    if ints == intm  # if the id for shear part is the same as the id for membrane part
                        HGs  = HGm   # -> reference the membrane part matrices
                        HξGs = HξGm  
                        HηGs = HηGm
                        HζGs = HζGm
                        Ws   = Wm 
                        gpes = gpem
                        iHGs = iHGm
                        XξGs = XξGm # ∂x/∂ξ matrix
                        YξGs = YξGm # ∂y/∂ξ matrix
                        ZξGs = ZξGm # ∂z/∂ξ matrix
                        XηGs = XηGm # ∂x/∂η matrix
                        YηGs = YηGm # ∂y/∂η matrix
                        ZηGs = ZηGm # ∂z/∂η matrix 
                    else # if the id for shear part is different from the id for membrane part
                        HGs  = mesh.elem[eType].HG[ints]  # shape function values at Gauss points for shear part 
                        HξGs = mesh.elem[eType].HξG[ints] # ∂h/∂ξ derivatives at Gauss points for shear part   
                        HηGs = mesh.elem[eType].HηG[ints] # ∂h/∂η derivatives at Gauss points for shear part   
                        HζGs = mesh.elem[eType].HζG[ints] # ∂h/∂ζ derivatives at Gauss points for shear part   
                        Ws   = mesh.elem[eType].GW[ints]  # Gauss weights at Gauss points for shear part
                        gpes = length(Ws)                 # number of Gauss points per elements for shear part
                        iHGs = pinv(HGs)
                        XξGs = X * HξGs # ∂x/∂ξ matrix
                        YξGs = Y * HξGs # ∂y/∂ξ matrix
                        ZξGs = Z * HξGs # ∂z/∂ξ matrix
                        XηGs = X * HηGs # ∂x/∂η matrix
                        YηGs = Y * HηGs # ∂y/∂η matrix
                        ZηGs = Z * HηGs # ∂z/∂η matrix                  
                    end
                    ϵs   = zeros(Float64, 4, tNum)
                    σs   = zeros(Float64, 4, tNum)
                    ϵsζ  = zeros(Float64, 2, tNum)
                    σsζ  = zeros(Float64, 2, tNum) 
                    if :En ∈ fields  || :Ee ∈ fields           
                        E1gs = zeros(Float64,2,tNum,gpes)
                    end
                    if :Sn ∈ fields  || :Se ∈ fields  
                        S1gs = zeros(Float64,2,tNum,gpes)
                    end
                else
                    gpes = 0
                end

                if intt > 0
                    # -> for torsional part of the stiffness matrix ----------------------------------------------------------
                    if intt == intm # if the id for torsional part is the same as the id for membrane part
                        HGt  = HGm  # -> reference the membrane part matrices
                        HξGt = HξGm 
                        HηGt = HηGm
                        HζGt = HζGm
                        Wt   = Wm 
                        gpet = gpem
                        iHGt = iHGm
                        XξGt = XξGm # ∂x/∂ξ matrix
                        YξGt = YξGm # ∂y/∂ξ matrix
                        ZξGt = ZξGm # ∂z/∂ξ matrix
                        XηGt = XηGm # ∂x/∂η matrix
                        YηGt = YηGm # ∂y/∂η matrix
                        ZηGt = ZηGm # ∂z/∂η matrix   
                    elseif intt == ints # if the id for torsional part is the same as the id for shear part
                        HGt  = HGs      # -> reference the shear part matrices
                        HξGt = HξGs 
                        HηGt = HηGs
                        HζGt = HζGs
                        Wt   = Ws 
                        gpet = gpes 
                        iHGt = iHGs
                        XξGt = XξGs # ∂x/∂ξ matrix
                        YξGt = YξGs # ∂y/∂ξ matrix
                        ZξGt = ZξGs # ∂z/∂ξ matrix
                        XηGt = XηGs # ∂x/∂η matrix
                        YηGt = YηGs # ∂y/∂η matrix
                        ZηGt = ZηGs # ∂z/∂η matrix                       
                    else  # if the id for torsional part is different from the id for membrane or shear part
                        HGt  = mesh.elem[eType].HG[intt]  # shape function values at Gauss points for torsional part 
                        HξGt = mesh.elem[eType].HξG[intt] # ∂h/∂ξ derivatives at Gauss points for torsional part   
                        HηGt = mesh.elem[eType].HηG[intt] # ∂h/∂η derivatives at Gauss points for torsional part  
                        HζGt = mesh.elem[eType].HζG[intt] # ∂h/∂ζ derivatives at Gauss points for torsional part
                        Wt   = mesh.elem[eType].GW[intt]  # Gauss weights at Gauss points for torsional part
                        gpet = length(Wt)                 # number of Gauss points per elements for torsional part
                        iHGt = pinv(HGt)
                        XξGt = X * HξGt # ∂x/∂ξ matrix
                        YξGt = Y * HξGt # ∂y/∂ξ matrix
                        ZξGt = Z * HξGt # ∂z/∂ξ matrix
                        XηGt = X * HηGt # ∂x/∂η matrix
                        YηGt = Y * HηGt # ∂y/∂η matrix
                        ZηGt = Z * HηGt # ∂z/∂η matrix                    
                    end
                    ϵt   = zeros(Float64, 1, tNum)
                    σt   = zeros(Float64, 1, tNum)
                    if :En ∈ fields  || :Ee ∈ fields           
                        E1gt = zeros(Float64,1,tNum,gpet)
                    end
                    if :Sn ∈ fields  || :Se ∈ fields  
                        S1gt = zeros(Float64,1,tNum,gpet)
                    end
                else
                    gpet = 0 
                end 
                
                # initialization of matrices:
                Bm = zeros(Float64, 6, nonpe*pdim) # B matrix for membrane+bending part
                Bs = zeros(Float64, 4, nonpe*pdim) # B matrix for shear part
                Bt = zeros(Float64, 1, nonpe*pdim) # B matrix for torsional part  
                Bmζ = zeros(Float64, 3, nonpe*pdim) # B matrix for membrane+bending part 
                Bsζ = zeros(Float64, 2, nonpe*pdim) # B matrix for shear part
                K1 = zeros(Float64, pdim*nonpe, pdim*nonpe) # initialization of stiffness matrix for eth element  
                DmBm = zeros(size(Bm))
                BmTDmBm = zeros(size(K1))
                DsBs = zeros(size(Bs))
                BsTDsBs = zeros(size(K1))
                BtTBt = zeros(size(K1))
                if :En ∈ fields  || :Ee ∈ fields           
                    E1n = zeros(Float64,6,tNum,nonpe)
                end
                if :Sn ∈ fields  || :Se ∈ fields  
                    S1n = zeros(Float64,6,tNum,nonpe)
                end
                if :Ue ∈ fields
                    UET = zeros(Float64,1,tNum,noe)
                    U1  = zeros(Float64,1,tNum)
                end
                
                i = 1
                @inbounds for e in 1:noe  # for each element of the current element type from current section
                    fill!(K1, 0.0)  # resetting stiffness matrix for the eth element  
                    if :En ∈ fields  || :Ee ∈ fields 
                        fill!(E1gm, 0.0)
                        fill!(E1gs, 0.0)
                        fill!(E1gt, 0.0)
                    end
                    if :Sn ∈ fields  || :Se ∈ fields  
                        fill!(S1gm, 0.0)
                        fill!(S1gs, 0.0)
                        fill!(S1gt, 0.0)
                    end
                    if :Ue ∈ fields
                        fill!(U1, 0.0) 
                    end
                    n1 = nn2[e,:]
                    q1 = q[n1,:]  # extracting nodal displacements associated with the eth element
                
                    @inbounds for g in 1:gpem  # for each Gauss point of the membrane part

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        # -> tangent vectors at the gth point of eth element
                        xξ = XξGm[e,g]  # ∂x/∂ξ at the gth point of eth element (a1x)
                        yξ = YξGm[e,g]  # ∂y/∂ξ at the gth point of eth element (a1y)
                        zξ = ZξGm[e,g]  # ∂z/∂ξ at the gth point of eth element (a1z)
                        xη = XηGm[e,g]  # ∂x/∂η at the gth point of eth element (a2x)
                        yη = YηGm[e,g]  # ∂y/∂η at the gth point of eth element (a2y)
                        zη = ZηGm[e,g]  # ∂z/∂η at the gth point of eth element (a2z)

                        # generating local coordinate system at the gth point of the eth element:
                        e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)
                        xζ = e3x .*0.5 .* t # ∂x/∂η at the gth point of eth element
                        yζ = e3y .*0.5 .* t # ∂y/∂η at the gth point of eth element
                        zζ = e3z .*0.5 .* t # ∂z/∂η at the gth point of eth element

                        # Jacobi determinant and elements of inverse at the gth point of the eth element:
                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)
                                                    
                        @inbounds for n in 1:nonpe      # for each node of the eth element

                            # shape function and it's derivatives with respect to x,y,z:
                            h  = HGm[n,g] # shape function value at the nth node
                            hx = HξGm[n,g] * ξx + HηGm[n,g] * ηx + HζGm[n,g] * ζx # ∂h/∂x 
                            hy = HξGm[n,g] * ξy + HηGm[n,g] * ηy + HζGm[n,g] * ζy # ∂h/∂y
                            hz = HξGm[n,g] * ξz + HηGm[n,g] * ηz + HζGm[n,g] * ζz # ∂h/∂z  

                            # gradient operators:
                            Lh1 = e1x*hx + e1y*hy + e1z*hz
                            Lh2 = e2x*hx + e2y*hy + e2z*hz
                            Lh3 = e3x*hx + e3y*hy + e3z*hz
                            Lζ1 = e1x*ζx + e1y*ζy + e1z*ζz
                            Lζ2 = e2x*ζx + e2y*ζy + e2z*ζz
                            Lζ3 = e3x*ζx + e3y*ζy + e3z*ζz

                            # preparing θᵀΦ product:
                            θTΦ11 = e1z*e3y - e1y*e3z
                            θTΦ12 = e1x*e3z - e1z*e3x
                            θTΦ13 = e1y*e3x - e1x*e3y
                            θTΦ21 = e2z*e3y - e2y*e3z
                            θTΦ22 = e2x*e3z - e2z*e3x
                            θTΦ23 = e2y*e3x - e2x*e3y
                            θTΦ31 = e3z*e3y - e3y*e3z
                            θTΦ32 = e3x*e3z - e3z*e3x
                            θTΦ33 = e3y*e3x - e3x*e3y

                            # elements of the Bm matrix associated with the strain contribution
                            # of the inplane displacements at the nth node
                            B1m11 = Lh1*e1x
                            B1m12 = Lh1*e1y
                            B1m13 = Lh1*e1z
                            B1m21 = Lh2*e2x
                            B1m22 = Lh2*e2y
                            B1m23 = Lh2*e2z
                            B1m31 = Lh2*e1x + Lh1*e2x
                            B1m32 = Lh2*e1y + Lh1*e2y
                            B1m33 = Lh2*e1z + Lh1*e2z

                            B2m11 = 0.5*t*h*Lζ1*θTΦ11
                            B2m12 = 0.5*t*h*Lζ1*θTΦ12
                            B2m13 = 0.5*t*h*Lζ1*θTΦ13
                            B2m21 = 0.5*t*h*Lζ2*θTΦ21
                            B2m22 = 0.5*t*h*Lζ2*θTΦ22
                            B2m23 = 0.5*t*h*Lζ2*θTΦ23
                            B2m31 = 0.5*t*h*(Lζ2*θTΦ11 + Lζ1*θTΦ21)
                            B2m32 = 0.5*t*h*(Lζ2*θTΦ12 + Lζ1*θTΦ22)
                            B2m33 = 0.5*t*h*(Lζ2*θTΦ13 + Lζ1*θTΦ23)

                            # elements of the Bm matrix associated with the strain contribution
                            # of the rotations at the nth node (also includes curvature effect)
                            B3m11 = 0.5*t*(B1m13*e3y - B1m12*e3z)
                            B3m12 = 0.5*t*(B1m11*e3z - B1m13*e3x)
                            B3m13 = 0.5*t*(B1m12*e3x - B1m11*e3y)
                            B3m21 = 0.5*t*(B1m23*e3y - B1m22*e3z)
                            B3m22 = 0.5*t*(B1m21*e3z - B1m23*e3x)
                            B3m23 = 0.5*t*(B1m22*e3x - B1m21*e3y)
                            B3m31 = 0.5*t*(B1m33*e3y - B1m32*e3z)
                            B3m32 = 0.5*t*(B1m31*e3z - B1m33*e3x)
                            B3m33 = 0.5*t*(B1m32*e3x - B1m31*e3y)

                            # populating Bm matrix:
                            Bm[1,(n-1)*pdim+1] = B1m11
                            Bm[1,(n-1)*pdim+2] = B1m12
                            Bm[1,(n-1)*pdim+3] = B1m13
                            Bm[2,(n-1)*pdim+1] = B1m21
                            Bm[2,(n-1)*pdim+2] = B1m22
                            Bm[2,(n-1)*pdim+3] = B1m23  
                            Bm[3,(n-1)*pdim+1] = B1m31
                            Bm[3,(n-1)*pdim+2] = B1m32
                            Bm[3,(n-1)*pdim+3] = B1m33     
                            Bm[4,(n-1)*pdim+4] = B3m11
                            Bm[4,(n-1)*pdim+5] = B3m12
                            Bm[4,(n-1)*pdim+6] = B3m13
                            Bm[5,(n-1)*pdim+4] = B3m21
                            Bm[5,(n-1)*pdim+5] = B3m22
                            Bm[5,(n-1)*pdim+6] = B3m23
                            Bm[6,(n-1)*pdim+4] = B3m31
                            Bm[6,(n-1)*pdim+5] = B3m32
                            Bm[6,(n-1)*pdim+6] = B3m33  

                            # populating Bm matrix:
                            Bmζ[1,(n-1)*pdim+1] = B1m11
                            Bmζ[1,(n-1)*pdim+2] = B1m12
                            Bmζ[1,(n-1)*pdim+3] = B1m13
                            Bmζ[2,(n-1)*pdim+1] = B1m21
                            Bmζ[2,(n-1)*pdim+2] = B1m22
                            Bmζ[2,(n-1)*pdim+3] = B1m23  
                            Bmζ[3,(n-1)*pdim+1] = B1m31
                            Bmζ[3,(n-1)*pdim+2] = B1m32
                            Bmζ[3,(n-1)*pdim+3] = B1m33     
                            Bmζ[1,(n-1)*pdim+4] = ζ*B3m11
                            Bmζ[1,(n-1)*pdim+5] = ζ*B3m12
                            Bmζ[1,(n-1)*pdim+6] = ζ*B3m13
                            Bmζ[2,(n-1)*pdim+4] = ζ*B3m21
                            Bmζ[2,(n-1)*pdim+5] = ζ*B3m22
                            Bmζ[2,(n-1)*pdim+6] = ζ*B3m23
                            Bmζ[3,(n-1)*pdim+4] = ζ*B3m31
                            Bmζ[3,(n-1)*pdim+5] = ζ*B3m32
                            Bmζ[3,(n-1)*pdim+6] = ζ*B3m33                              
                        end
                        # updating stiffness matrix with the membrane+bending part at the gth Gauss point
                        mul!(ϵm, Bm, q1)              
                        mul!(σm, Dm, ϵm) 
                        if :Ue ∈ fields
                            U1 += 0.5*sum(σm.*ϵm,dims=1) * det * Wm[g]  
                        end

                        mul!(ϵmζ, Bmζ, q1)              
                        mul!(σmζ, Dm1, ϵmζ) 
                        if :En ∈ fields  || :Ee ∈ fields 
                            E1gm[:,:,g] = ϵmζ
                        end
                        if :Sn ∈ fields  || :Se ∈ fields  
                            S1gm[:,:,g] = σmζ
                        end                       
                    end

                    @inbounds for g in 1:gpes  # for each Gauss point of the shear part

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        # -> tangent vectors at the gth point of eth element
                        xξ = XξGs[e,g]  # ∂x/∂ξ at the gth point of eth element (a1x)
                        yξ = YξGs[e,g]  # ∂y/∂ξ at the gth point of eth element (a1y)
                        zξ = ZξGs[e,g]  # ∂z/∂ξ at the gth point of eth element (a1z)
                        xη = XηGs[e,g]  # ∂x/∂η at the gth point of eth element (a2x)
                        yη = YηGs[e,g]  # ∂y/∂η at the gth point of eth element (a2y)
                        zη = ZηGs[e,g]  # ∂z/∂η at the gth point of eth element (a2z)

                        # generating local coordinate system at the gth point of the eth element:
                        e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)
                        xζ = e3x .*0.5 .* t # ∂x/∂η at the gth point of eth element
                        yζ = e3y .*0.5 .* t # ∂y/∂η at the gth point of eth element
                        zζ = e3z .*0.5 .* t # ∂z/∂η at the gth point of eth element

                        # Jacobi determinant and elements of inverse at the gth point of the eth element:
                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)
                                                    
                        @inbounds for n in 1:nonpe      # for each node of the eth element

                            # shape function and it's derivatives with respect to x,y,z:
                            h  = HGs[n,g] # shape function value at the nth node
                            hx = HξGs[n,g] * ξx + HηGs[n,g] * ηx + HζGs[n,g] * ζx # ∂h/∂x 
                            hy = HξGs[n,g] * ξy + HηGs[n,g] * ηy + HζGs[n,g] * ζy # ∂h/∂y
                            hz = HξGs[n,g] * ξz + HηGs[n,g] * ηz + HζGs[n,g] * ζz # ∂h/∂z 

                            # gradient operators:
                            Lh1 = e1x*hx + e1y*hy + e1z*hz
                            Lh2 = e2x*hx + e2y*hy + e2z*hz
                            Lh3 = e3x*hx + e3y*hy + e3z*hz
                            Lζ1 = e1x*ζx + e1y*ζy + e1z*ζz
                            Lζ2 = e2x*ζx + e2y*ζy + e2z*ζz
                            Lζ3 = e3x*ζx + e3y*ζy + e3z*ζz

                            # preparing θᵀΦ product:
                            θTΦ11 = e1z*e3y - e1y*e3z
                            θTΦ12 = e1x*e3z - e1z*e3x
                            θTΦ13 = e1y*e3x - e1x*e3y
                            θTΦ21 = e2z*e3y - e2y*e3z
                            θTΦ22 = e2x*e3z - e2z*e3x
                            θTΦ23 = e2y*e3x - e2x*e3y
                            θTΦ31 = e3z*e3y - e3y*e3z
                            θTΦ32 = e3x*e3z - e3z*e3x
                            θTΦ33 = e3y*e3x - e3x*e3y
                            
                            # calculating elements of the Bs matrix:
                            B1s11 = Lh3*e1x + Lh1*e3x
                            B1s12 = Lh3*e1y + Lh1*e3y
                            B1s13 = Lh3*e1z + Lh1*e3z
                            B1s21 = Lh3*e2x + Lh2*e3x
                            B1s22 = Lh3*e2y + Lh2*e3y
                            B1s23 = Lh3*e2z + Lh2*e3z

                            B2s11 = 0.5*t*h*(Lζ3*θTΦ11 + Lζ1*θTΦ31)
                            B2s12 = 0.5*t*h*(Lζ3*θTΦ12 + Lζ1*θTΦ32)
                            B2s13 = 0.5*t*h*(Lζ3*θTΦ13 + Lζ1*θTΦ33)
                            B2s21 = 0.5*t*h*(Lζ3*θTΦ21 + Lζ2*θTΦ31)
                            B2s22 = 0.5*t*h*(Lζ3*θTΦ22 + Lζ2*θTΦ32)
                            B2s23 = 0.5*t*h*(Lζ3*θTΦ23 + Lζ2*θTΦ33)

                            B3s11 = 0.5*t*(Lh3*θTΦ11 + Lh1*θTΦ31)
                            B3s12 = 0.5*t*(Lh3*θTΦ12 + Lh1*θTΦ32)
                            B3s13 = 0.5*t*(Lh3*θTΦ13 + Lh1*θTΦ33)
                            B3s21 = 0.5*t*(Lh3*θTΦ21 + Lh2*θTΦ31)
                            B3s22 = 0.5*t*(Lh3*θTΦ22 + Lh2*θTΦ32)
                            B3s23 = 0.5*t*(Lh3*θTΦ23 + Lh2*θTΦ33)

                            # populating Bs matrix:
                            Bs[1,(n-1)*pdim+1] = B1s11
                            Bs[1,(n-1)*pdim+2] = B1s12
                            Bs[1,(n-1)*pdim+3] = B1s13
                            Bs[2,(n-1)*pdim+1] = B1s21
                            Bs[2,(n-1)*pdim+2] = B1s22
                            Bs[2,(n-1)*pdim+3] = B1s23
                            Bs[1,(n-1)*pdim+4] = B2s11
                            Bs[1,(n-1)*pdim+5] = B2s12
                            Bs[1,(n-1)*pdim+6] = B2s13
                            Bs[2,(n-1)*pdim+4] = B2s21
                            Bs[2,(n-1)*pdim+5] = B2s22
                            Bs[2,(n-1)*pdim+6] = B2s23  
                            Bs[3,(n-1)*pdim+4] = B3s11
                            Bs[3,(n-1)*pdim+5] = B3s12
                            Bs[3,(n-1)*pdim+6] = B3s13                          
                            Bs[4,(n-1)*pdim+4] = B3s21
                            Bs[4,(n-1)*pdim+5] = B3s22
                            Bs[4,(n-1)*pdim+6] = B3s23   

                            Bsζ[1,(n-1)*pdim+1] = B1s11
                            Bsζ[1,(n-1)*pdim+2] = B1s12
                            Bsζ[1,(n-1)*pdim+3] = B1s13
                            Bsζ[2,(n-1)*pdim+1] = B1s21
                            Bsζ[2,(n-1)*pdim+2] = B1s22
                            Bsζ[2,(n-1)*pdim+3] = B1s23
                            Bsζ[1,(n-1)*pdim+4] = B2s11 + ζ*B3s11
                            Bsζ[1,(n-1)*pdim+5] = B2s12 + ζ*B3s12
                            Bsζ[1,(n-1)*pdim+6] = B2s13 + ζ*B3s13  
                            Bsζ[2,(n-1)*pdim+4] = B2s21 + ζ*B3s21
                            Bsζ[2,(n-1)*pdim+5] = B2s22 + ζ*B3s22
                            Bsζ[2,(n-1)*pdim+6] = B2s23 + ζ*B3s23                              
                        end
                        # updating stiffness matrix with the shear part at the gth Gauss point
                        mul!(ϵs, Bs, q1)              
                        mul!(σs, Ds, ϵs) 
                        if :Ue ∈ fields
                            U1 += 0.5*sum(σs.*ϵs,dims=1) * det * Ws[g]  
                        end

                        mul!(ϵsζ, Bsζ, q1)              
                        mul!(σsζ, Ds1, ϵsζ) 
                        if :En ∈ fields  || :Ee ∈ fields 
                            E1gs[:,:,g] = ϵsζ
                        end
                        if :Sn ∈ fields  || :Se ∈ fields  
                            S1gs[:,:,g] = σsζ
                        end   
                    end
                    
                    @inbounds for g in 1:gpet  # for each Gauss point of the torsional part

                        # extracting Jacobi matrix elements at the gth point of the eth element:
                        # -> tangent vectors at the gth point of eth element
                        xξ = XξGt[e,g]  # ∂x/∂ξ at the gth point of eth element (a1x)
                        yξ = YξGt[e,g]  # ∂y/∂ξ at the gth point of eth element (a1y)
                        zξ = ZξGt[e,g]  # ∂z/∂ξ at the gth point of eth element (a1z)
                        xη = XηGt[e,g]  # ∂x/∂η at the gth point of eth element (a2x)
                        yη = YηGt[e,g]  # ∂y/∂η at the gth point of eth element (a2y)
                        zη = ZηGt[e,g]  # ∂z/∂η at the gth point of eth element (a2z)

                        # generating local coordinate system at the gth point of the eth element:
                        e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z = csysgen(xξ, yξ, zξ, xη, yη, zη)
                        xζ = e3x .*0.5 .* t # ∂x/∂η at the gth point of eth element
                        yζ = e3y .*0.5 .* t # ∂y/∂η at the gth point of eth element
                        zζ = e3z .*0.5 .* t # ∂z/∂η at the gth point of eth element

                        # Jacobi determinant and elements of inverse at the gth point of the eth element:
                        det,ξx,ηx,ζx,ξy,ηy,ζy,ξz,ηz,ζz = invgen3(xξ,yξ,zξ,xη,yη,zη,xζ,yζ,zζ)
                                                    
                        @inbounds for n in 1:nonpe      # for each node of the eth element

                            # shape function and it's derivatives with respect to x,y,z:
                            h  = HGt[n,g] # shape function value at the nth node
                            hx = HξGt[n,g] * ξx + HηGt[n,g] * ηx + HζGt[n,g] * ζx # ∂h/∂x 
                            hy = HξGt[n,g] * ξy + HηGt[n,g] * ηy + HζGt[n,g] * ζy # ∂h/∂y
                            hz = HξGt[n,g] * ξz + HηGt[n,g] * ηz + HζGt[n,g] * ζz # ∂h/∂z                      

                            # gradient operators:
                            Lh1 = e1x*hx + e1y*hy + e1z*hz
                            Lh2 = e2x*hx + e2y*hy + e2z*hz

                            # calculating & populating Bt matrix:
                            Bt[(n-1)*pdim+1] = 0.5 * (Lh2*e1x - Lh1*e2x)
                            Bt[(n-1)*pdim+2] = 0.5 * (Lh2*e1y - Lh1*e2y)
                            Bt[(n-1)*pdim+3] = 0.5 * (Lh2*e1z - Lh1*e2z)
                            Bt[(n-1)*pdim+4] = h*e3x
                            Bt[(n-1)*pdim+5] = h*e3y
                            Bt[(n-1)*pdim+6] = h*e3z  
                        end
                        # updating stiffness matrix with the torsional part at the gth Gauss point
                        mul!(ϵt, Bt, q1)              
                        mul!(σt, (kt * G * t), ϵt) 
                        if :Ue ∈ fields
                            U1 += 0.5*sum(σt.*ϵt,dims=1) * det * Wt[g]  
                        end
                        if :En ∈ fields  || :Ee ∈ fields 
                            E1gt[:,:,g] = ϵt
                        end
                        if :Sn ∈ fields  || :Se ∈ fields  
                            S1gt[:,:,g] = σt
                        end 
                    end

                    if :Ue ∈ fields
                        UET[1,:,e] = U1
                    end

                    for t = 1:tNum
                        if :En ∈ fields  || :Ee ∈ fields 
                            E1n[[1,2,4],t,:] = view(E1gm,:,t,:) * iHGm
                            E1n[[5,6],t,:]   = view(E1gs,:,t,:) * iHGs
                            E1n[3,t,:]       = view(E1gt,:,t,:) * iHGt
                        end
                        if :Sn ∈ fields  || :Se ∈ fields  
                            S1n[[1,2,4],t,:] = view(S1gm,:,t,:) * iHGm
                          #  S1n[[5,6],t,:]   = view(S1gs,:,t,:) * iHGs
                          #  S1n[3,t,:]       = view(S1gt,:,t,:) * iHGt
                        end
                    end

                    if :Ee ∈ fields 
                        Ee.data = cat(Ee.data, E1n; dims=3)
                    end
                    if :Se ∈ fields 
                        Se.data = cat(Se.data, S1n; dims=3)  
                    end
                    if :En ∈ fields                                                              
                        En.data[:,:,nn[e,:]] += E1n
                    end
                    if :Sn ∈ fields  
                        Sn.data[:,:,nn[e,:]] += S1n   
                    end 
                end

                if :Ee ∈ fields 
                    append!(Ee.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Ee.eTags,eTags)  # adding element     
                end
                if :Se ∈ fields         
                    append!(Se.nonpe,fill(nonpe,noe))      # adding number of nodes per element numbers of current element type from current section to global vector
                    append!(Se.eTags,eTags)  # adding element 
                end
                if :Ue ∈ fields 
                    Ue.data = cat(Ue.data, UET; dims=3) 
                    append!(Ue.eTags,eTags)  # adding element 
                end

                if :En ∈ fields || :Sn ∈ fields # if nodal field is required
                    for x in nn
                        epn[x] += 1  # updating element per node counter
                    end
                end
            end
        elseif secData isa Beam3D  # 3D Beam section
            # not yet implemented
        elseif secData isa Beam2D  # 2D Beam section  
            # not yet implemented
        elseif secData isa Beam1D  # 1D Beam section  
            # not yet implemented             
        elseif secData isa Truss3D # 3D Truss section
            # not yet implemented
        elseif secData isa Truss2D # 2D Truss section
            # not yet implemented
        elseif secData isa Truss1D # 1D Truss section
            # not yet implemented                    
        end
    end

    if :En ∈ fields || :Sn ∈ fields
        epn[epn .== 0] .= 1  # replacing epn = 0 values to 1 for division
        if :En ∈ fields
            En.data ./= reshape(epn, 1, 1, :)  # dividing stress values of 3D nodal stress matrix with number of nodes per elements
        end
        if :Sn ∈ fields
            Sn.data ./= reshape(epn, 1, 1, :)  # dividing stress values of 3D nodal stress matrix with number of nodes per elements
        end
    end

    return Sn,Se,En,Ee,Ue
end

function addView(mesh,DoFs,field,t=[0.0];comp=:all,name=:default,scf=100,show=1,depict=:disp)
    modelName = mesh.name
    non       = mesh.non  # number of nodes
    nTags     = collect(1:non)
    tNum      = length(t)
    
    # exclude orphan nodes:
    ActiveNodes= findall(i -> any(!iszero, DoFs[i, :]), 1:size(DoFs, 1))
    nTags = nTags[ActiveNodes]
    DoFs = DoFs[ActiveNodes,:]
    noan = length(nTags)  # number of active nodes

    comp = to_vector(comp)
    name = to_vector(name)
    pdim = size(DoFs,2)

    if field isa ElementNodeData
        r = size(field.data,1)  # number of field components
        eTags = field.eTags
        noe = length(eTags)
        nonpe = field.nonpe
        fdata = field.data
        ncount = [0;cumsum(nonpe)]
        if :p1 ∈ comp || :p2 ∈ comp || :p3 ∈ comp || :p ∈ comp || :Mohr ∈ comp || :Coulomb ∈ comp
            l = size(field.data,3)
            pdata = Array{Float64, 3}(undef, 3, tNum, l)
            for n in 1:l
                for t in 1:tNum
                    # Extract components
                    sxx, syy, szz, sxy, sxz, syz = fdata[:, t, n]

                    # Reconstruct full tensor
                    S = [sxx  sxy  sxz;
                         sxy  syy  syz;
                         sxz  syz  szz]

                    # Compute eigenvalues (principal stresses)
                    pdata[:, t, n] .= eigvals(S) |> sort
                end
            end
        end
        for c in eachindex(comp)
            comp1 = comp[c]
            if name[1] == :default
                name1 = string(field.name, "_", comp1)
            else
                name1 = name[c]
            end
            viewTag   = gmsh.view.add(name1)
            if comp1 == :xx || comp1 == :x
                fdata1 = fdata[1:1,:,:]
            elseif (comp1 == :yy || comp1 == :y) && r > 1
                fdata1 = fdata[2:2,:,:]
            elseif (comp1 == :zz || comp1 == :z) && r > 2
                fdata1 = fdata[3:3,:,:]
            elseif (comp1 == :xy || comp1 == :yx) && r > 3
                fdata1 = fdata[4:4,:,:]
            elseif (comp1 == :xz || comp1 == :zx) && r > 4
                fdata1 = fdata[5:5,:,:]      
            elseif (comp1 == :yz || comp1 == :zy) && r > 5
                fdata1 = fdata[6:6,:,:]  
            elseif (comp1 == :Mises || comp1 == :vonMises || comp1 == :HMH) && r == 6
                fdata1 = sqrt.(0.5*((fdata[1:1,:,:]-fdata[2:2,:,:]).^2+
                                    (fdata[1:1,:,:]-fdata[3:3,:,:]).^2+
                                    (fdata[2:2,:,:]-fdata[3:3,:,:]).^2+
                                    6*((fdata[4:4,:,:]).^2+(fdata[5:5,:,:]).^2+(fdata[6:6,:,:]).^2)))
            elseif comp1 == :p1
                fdata1 = pdata[1:1,:,:]
            elseif comp1 == :p2
                fdata1 = pdata[2:2,:,:]
            elseif comp1 == :p3
                fdata1 = pdata[3:3,:,:]                
            elseif comp1 == :p
                fdata1 = pdata
            elseif comp1 == :Mohr
                fdata1 = pdata[3:3,:,:]-pdata[1:1,:,:]
            elseif comp1 == :Coulomb
                fdata1 = max.(abs.(pdata[1:1,:,:]),abs.(pdata[3:3,:,:]))
            else
                fdata1 = fdata
            end 

            r1 = size(fdata1,1)
            c1 = r1==6 ? 9 : r1==3 ? 3 : r1==1 ? 1 : error("not supported")
            data  = [zeros(nonpe[i]*c1) for i in 1:noe]
            for i = 1:tNum
                for j = 1:noe
                    h = ncount[j]
                    for k = 0:nonpe[j]-1   
                        if r1 == 6
                            data[j][k*9+1] = fdata1[1,i,h+k+1]   # Txx
                            data[j][k*9+5] = fdata1[2,i,h+k+1]   # Tyy
                            data[j][k*9+9] = fdata1[3,i,h+k+1]   # Tzz
                            data[j][k*9+2] = data[j][k*9+4] = fdata1[4,i,h+k+1]   # Txy = Tyx
                            data[j][k*9+6] = data[j][k*9+8] = fdata1[5,i,h+k+1]   # Tyz = Tzy
                            data[j][k*9+3] = data[j][k*9+7] = fdata1[6,i,h+k+1]   # Txz = Tzx   
                        elseif r1 == 3
                            data[j][k*3+1] = fdata1[1,i,h+k+1]   # Vx
                            data[j][k*3+2] = fdata1[2,i,h+k+1]   # vy
                            data[j][k*3+3] = fdata1[3,i,h+k+1]   # Vz  
                        else
                            data[j][k+1]   = fdata1[1,i,h+k+1]   
                        end
                    end
                end
                gmsh.view.addModelData(viewTag, i-1, modelName, "ElementNodeData", eTags, data, t[i], c1)  
            end
            gmsh.view.option.setNumber(viewTag, "IntervalsType", 3)
            gmsh.view.option.setNumber(viewTag, "Visible", show)   
        end
    elseif field isa NodeData 
        r = size(field.data,1)  # number of field components
        fdata = field.data
        if :p1 ∈ comp || :p2 ∈ comp || :p3 ∈ comp || :p ∈ comp || :Mohr ∈ comp || :Coulomb ∈ comp
            l = size(field.data,3)
            pdata = Array{Float64, 3}(undef, 3, tNum, l)
            for n in 1:l
                for t in 1:tNum
                    # Extract components
                    sxx, syy, szz, sxy, sxz, syz = fdata[:, t, n]

                    # Reconstruct full tensor
                    S = [sxx  sxy  sxz;
                        sxy  syy  syz;
                        sxz  syz  szz]

                    # Compute eigenvalues (principal stresses)
                    pdata[:, t, n] .= eigvals(S) |> sort
                end
            end
        end
        for c in eachindex(comp)
            comp1 = comp[c]
            if name[1] == :default
                name1 = string(field.name, "_", comp1)
            else
                name1 = name[c]
            end
            viewTag   = gmsh.view.add(name1)
            if comp1 == :xx || comp1 == :x
                fdata1 = fdata[1:1,:,:]
            elseif (comp1 == :yy || comp1 == :y) && r > 1
                fdata1 = fdata[2:2,:,:]
            elseif (comp1 == :zz || comp1 == :z) && r > 2
                fdata1 = fdata[3:3,:,:]
            elseif (comp1 == :xy || comp1 == :yx) && r > 3
                fdata1 = fdata[4:4,:,:]
            elseif (comp1 == :xz || comp1 == :zx) && r > 4
                fdata1 = fdata[5:5,:,:]      
            elseif (comp1 == :yz || comp1 == :zy) && r > 5
                fdata1 = fdata[6:6,:,:]  
            elseif (comp1 == :Mises || comp1 == :vonMises || comp1 == :HMH) && r == 6
                fdata1 = sqrt.(0.5*((fdata[1:1,:,:]-fdata[2:2,:,:]).^2+
                                    (fdata[1:1,:,:]-fdata[3:3,:,:]).^2+
                                    (fdata[2:2,:,:]-fdata[3:3,:,:]).^2+
                                    6*((fdata[4:4,:,:]).^2+(fdata[5:5,:,:]).^2+(fdata[6:6,:,:]).^2)))
            elseif comp1 == :p1
                fdata1 = pdata[1:1,:,:]
            elseif comp1 == :p2
                fdata1 = pdata[2:2,:,:]
            elseif comp1 == :p3
                fdata1 = pdata[3:3,:,:]               
            elseif comp1 == :p
                fdata1 = pdata
            elseif comp1 == :Mohr
                fdata1 = pdata[3:3,:,:]-pdata[1:1,:,:]
            elseif comp1 == :Coulomb
                fdata1 = max.(abs.(pdata[1:1,:,:]),abs.(pdata[3:3,:,:]))
            else
                fdata1 = fdata
            end

            r1 = size(fdata1,1)
            c1 = r1==6 ? 9 : r1==3 ? 3 : r1==1 ? 1 : error("not supported")

            data = zeros(c1*non)
            for i = 1:tNum
                for j = 0:non-1  
                    if r1 == 6
                        data[j*9+1] = fdata1[1,i,j+1]   # Txx    
                        data[j*9+5] = fdata1[2,i,j+1]   # Tyy
                        data[j*9+9] = fdata1[3,i,j+1]   # Tzz   
                        data[j*9+2] = data[j*9+4] = fdata1[4,i,j+1]   # Txy = Tyx
                        data[j*9+6] = data[j*9+8] = fdata1[5,i,j+1]   # Tyz = Tzy
                        data[j*9+3] = data[j*9+7] = fdata1[6,i,j+1]   # Txz = Tzx  
                    elseif r1 == 3
                        data[j*3+1] = fdata1[1,i,j+1]   # Vx   
                        data[j*3+2] = fdata1[1,i,j+2]   # Vy
                        data[j*3+3] = fdata1[1,i,j+3]   # Vz
                    else
                        data[j+1] = fdata1[1,i,j+1]
                    end
                end
                gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],c1)
            end 
            gmsh.view.option.setNumber(viewTag, "IntervalsType", 3)
            gmsh.view.option.setNumber(viewTag, "Visible", show)   
        end
    elseif field isa ElementData 
        fdata1 = field.data
        eTags = field.eTags
        noe = length(eTags)
        for c in eachindex(comp)
            comp1 = comp[c]
            if name[1] == :default
                name1 = string(field.name, "_", comp1)
            else
                name1 = name[c]
            end
            viewTag   = gmsh.view.add("energy")
            data  = [[0.0] for i in 1:noe]
           # data = zeros(noe)
            for i = 1:tNum
                for j = 0:noe-1    
                    data[j+1] = [fdata1[1,i,j+1]]
                end  
                gmsh.view.addModelData(viewTag, i-1, modelName, "ElementData", eTags, data, t[i], 1)  
            end
            gmsh.view.option.setNumber(viewTag, "IntervalsType", 3)
            gmsh.view.option.setNumber(viewTag, "Visible", show)   
        end
    else
        if pdim == 3
            for c in eachindex(comp)
                comp1 = comp[c]
                if name[1] == :default
                    name1 = string("U_",comp1)
                else
                    name1 = name[c]
                end
                viewTag   = gmsh.view.add(name1)                    
                if comp1 == :all
                    data = zeros(3*noan)
                    for i = 1:tNum
                        for j = 0:noan-1
                            data[3*j+1] = field[DoFs[j+1,1],i]
                            data[3*j+2] = field[DoFs[j+1,2],i]
                        end
                        gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],3)
                    end
                elseif comp1 == :abs
                    data = zeros(noan)
                    for i = 1:tNum
                        data = sqrt(field[DoFs[:,1],i].^2+field[DoFs[:,2],i].^2)
                        gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
                    end
                elseif comp1 == :x
                    data = field[DoFs[:,1]]
                    gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
                elseif comp1 == :y
                    data = field[DoFs[:,2]]
                    gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
                else                  
                    data = field[DoFs[:,comp1]]
                    gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
                end
                
                vt = depict == :disp ? 5 : 2
                gmsh.view.option.setNumber(viewTag, "IntervalsType", 3)
                gmsh.view.option.setNumber(viewTag, "VectorType", vt)
                gmsh.view.option.setNumber(viewTag, "DisplacementFactor",scf)
                gmsh.view.option.setNumber(viewTag, "Visible", show)
            end
        elseif pdim == 6
            for c in eachindex(comp)
                comp1 = comp[c]
                if name[1] == :default
                    name1 = string("U_",comp1)
                else
                    name1 = name[c]
                end
                viewTag   = gmsh.view.add(name1)                  
                if comp1 == :all
                    data = zeros(3*noan)
                    for i = 1:tNum
                        for j = 0:noan-1
                            data[3*j+1] = field[DoFs[j+1,1],i]
                            data[3*j+2] = field[DoFs[j+1,2],i]
                            data[3*j+3] = field[DoFs[j+1,3],i]
                        end
                        gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],3)    
                    end
                elseif comp1 == :abs
                    data = zeros(nc*noan)
                    for i = 1:tNum
                        data = sqrt(field[DoFs[:,1]].^2+field[DoFs[:,2]].^2+field[DoFs[:,3]].^2)
                        gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
                    end    
                elseif comp1 == :x
                    data = field[DoFs[:,1]]
                    gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
                elseif comp1 == :y
                    data = field[DoFs[:,2]]
                    gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)  
                elseif comp1 == :z
                    data = field[DoFs[:,3]]
                    gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)                                                    
                else
                    data = field[DoFs[:,comp1]]
                    gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
                end
                vt = depict == :disp ? 5 : 2
                gmsh.view.option.setNumber(viewTag, "IntervalsType", 3)
                gmsh.view.option.setNumber(viewTag, "VectorType", vt)
                gmsh.view.option.setNumber(viewTag, "DisplacementFactor",scf)
                gmsh.view.option.setNumber(viewTag, "Visible", show)
            end
        end

    end
end

# A function adding views for Vector filed in GMSH
function addVectorResult(mesh,vf,DoFs,name,t=[0.0];comp=:all,scf=100,show=1,depict=:disp)
    modelName = mesh.name
    non       = mesh.non  # number of nodes
    nTags     = collect(1:non)
    tNum      = length(t)
    viewTag   = gmsh.view.add(name)

    # exclude orphan nodes:
    ActiveNodes= findall(i -> any(!iszero, DoFs[i, :]), 1:size(DoFs, 1))
    nTags = nTags[ActiveNodes]
    DoFs = DoFs[ActiveNodes,:]
    noan = length(nTags)  # number of active nodes

    pdim = size(DoFs,2)
    if pdim == 3
        if comp == :all
            data = zeros(3*noan)
            for i = 1:tNum
                for j = 0:noan-1
                    data[3*j+1] = vf[DoFs[j+1,1]]
                    data[3*j+2] = vf[DoFs[j+1,2]]
                end
                gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],3)
            end
        elseif comp == :abs
            data = zeros(noan)
            for i = 1:tNum
                data = sqrt(vf[DoFs[:,1]].^2+vf[DoFs[:,2]].^2)
                gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
            end
        else
            data = vf[DoFs[:,comp]]
            gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
        end
    elseif pdim == 6
        if comp == :all
            data = zeros(3*noan)
            for i = 1:tNum
                for j = 0:noan-1
                    data[3*j+1] = vf[DoFs[j+1,1]]
                    data[3*j+2] = vf[DoFs[j+1,2]]
                    data[3*j+3] = vf[DoFs[j+1,3]]
                end
                gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],3)    
            end
        elseif comp == :abs
            data = zeros(nc*noan)
            for i = 1:tNum
                data = sqrt(vf[DoFs[:,1]].^2+vf[DoFs[:,2]].^2+vf[DoFs[:,3]].^2)
                gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
            end                
        else
            data = vf[DoFs[:,comp]]
            gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],1)
        end
    end 

    vt = depict == :disp ? 5 : 2
    gmsh.view.option.setNumber(viewTag, "IntervalsType", 3)
    gmsh.view.option.setNumber(viewTag, "VectorType", vt)
    gmsh.view.option.setNumber(viewTag, "DisplacementFactor",scf)
    gmsh.view.option.setNumber(viewTag, "Visible", show)
   # gmsh.write("modified_model.msh")
end

# A function adding views for Tensor filed in GMSH
function addTensorResult(mesh,T,name,t=[0.0];comp=:all,show=0)
    modelName = mesh.name
    non       = mesh.non  # number of nodes
    nTags     = collect(1:non)
    tNum      = length(t)
    viewTag   = gmsh.view.add(name)  

    if T isa ElementalTensorField
        if comp == :all
            eTags = T.eTags
            noe = length(eTags)
            nonpe = T.nonpe
            ncount = [0;cumsum(nonpe)]
            data  = [zeros(nonpe[i]*9) for i in 1:noe]

            for i = 1:tNum
                for j = 1:noe
                    h = ncount[j]
                    for k = 0:nonpe[j]-1
                        data[j][k*9+1] = T.tBlock[1,i,h+k+1]   # Txx
                        data[j][k*9+5] = T.tBlock[2,i,h+k+1]   # Tyy
                        data[j][k*9+9] = T.tBlock[3,i,h+k+1]   # Tzz
                        data[j][k*9+2] = data[j][k*9+4] = T.tBlock[4,i,h+k+1]   # Txy = Tyx
                        data[j][k*9+6] = data[j][k*9+8] = T.tBlock[5,i,h+k+1]   # Tyz = Tzy
                        data[j][k*9+3] = data[j][k*9+7] = T.tBlock[6,i,h+k+1]   # Txz = Tzx   
                    end
                end
                gmsh.view.addModelData(viewTag, i-1, modelName, "ElementNodeData", eTags, data, t[i], 9)
            end
        elseif comp == :vonMises
            # ...
        end
        # gmsh.view.option.setNumber(viewTag, "NormalRaise", 0.2)    
    elseif T isa NodalTensorField  
        if comp == :all
            data = zeros(9*non)
            for i = 1:tNum
                for j = 0:non-1  
                    data[j*9+1] = T.tBlock[1,i,j+1]   # Txx    
                    data[j*9+5] = T.tBlock[2,i,j+1]   # Tyy
                    data[j*9+9] = T.tBlock[3,i,j+1]   # Tzz   
                    data[j*9+2] = data[j*9+4] = T.tBlock[4,i,j+1]   # Txy = Tyx
                    data[j*9+6] = data[j*9+8] = T.tBlock[5,i,j+1]   # Tyz = Tzy
                    data[j*9+3] = data[j*9+7] = T.tBlock[6,i,j+1]   # Txz = Tzx  
                end
                gmsh.view.addHomogeneousModelData(viewTag,i-1,modelName,"NodeData",nTags,data,t[i],9)
            end
        elseif comp == :vonMises
            # ...
        end   
    end
    gmsh.view.option.setNumber(viewTag, "IntervalsType", 3)
    gmsh.view.option.setNumber(viewTag, "Visible", show)    
end           
            
# A function starting GMSH gui & finalizes it
function startGMSH()
    gmsh.fltk.run()
    gmsh.finalize()
end

# A function to convert a general scope to a selection structure
function scope2selection(mesh,scope)
    if scope isa Selector
        if scope.object == :n
            scope = nodeSelelection(mesh,scope.seldef,scope.selpar,scope.presel,scope.name,scope.id,scope.tol)
        else
            scope =  elementSelelection(mesh,scope.object,scope.seldef,scope.selpar,scope.presel,scope.name,scope.id,scope.tol)
        end
    elseif scope isa Symbol
        if scope == :n
            scope =  nodeSelelection(mesh)
        else
            scope =  elementSelelection(mesh,scope)
        end
    elseif scope isa Tuple && scope[1] isa Symbol && scope[2] isa Symbol && scope[1] != :n
            scope =  elementSelelection(mesh,scope[1],scope[2],scope[3])
    elseif  scope isa Tuple && scope[1] isa Symbol   
            if length(scope) == 2    
                scope = nodeSelelection(mesh,scope[1],scope[2])
            elseif length(scope) == 3
                scope = nodeSelelection(mesh,scope[2],scope[3])
            end
    elseif !isa(scope,Selection)
        scope = mesh.set[mesh.findset[scope]] 
    end

    return scope
end

# A function to find a set base on an setID
function set(mesh,setID)
    return y = mesh.set[mesh.findset[setID]] 
end

function get_max_entity_tag(dim::Int)
    entities = gmsh.model.getEntities(dim)
    isempty(entities) && return nothing
    return maximum(t -> t[2], entities)  # entities are (dim, tag) tuples
end

function get_max_node_tag()
    tags, _, _ = gmsh.model.mesh.getNodes()
    isempty(tags) && return nothing
    return maximum(tags)
end

function get_max_element_tag()
    _, elementTags, _ = gmsh.model.mesh.getElements()
    isempty(elementTags) && return 0
    return maximum(vcat(elementTags...))  # flatten and get max
end

function get_max_physical_group_tag()
    phys = gmsh.model.getPhysicalGroups()
    isempty(phys) && return nothing
    return maximum(t -> t[2], phys)  # physical groups are (dim, tag) tuples
end

to_vector(x) = x isa AbstractVector || x isa AbstractDict ? x : [x]

function to_vector2(A, n::Int)

    A_vec = A isa AbstractVector ? A : [A]  # convert to vector if A is not a vector

    if length(A_vec) == 1                   # if A has only 1 element
        return fill(A_vec[1], n)            # repeat the elements n times
    elseif length(A_vec) == n               # if A has n elements
        return A_vec                        # return the original A vector
    else                                    # otherwise: sho an error message
        error("Length of input must be 1 or $n, but got length $(length(A_vec))")
    end
end

function to_vector3(A, n::Int)

    if A isa Array && size(A,2) > 1
        A = [A[i:i,:] for i in 1:size(A, 1)]
    elseif  A isa Symbol || (A isa Array && (A[1] isa Symbol || !(A[1] isa Array)))
        A = [A]
    end

    if length(A) == 1                   # if A has only 1 element
        return fill(A[1], n)            # repeat the elements n times
    elseif length(A) == n               # if A has n elements
        return A                        # return the original A vector
    else                                # otherwise: sho an error message
        error("Length of input must be 1 or $n, but got length $(length(A))")
    end
    
end

function to_vectorvector(A, n::Int)

    A isa Array && !isempty(A) ? A = A : A = [A]
    if !(A[1] isa Array) && size(A,1) > 1
        A = [A[i,:] for i in 1:size(A,1)]
    end
    A[1] isa Array ? A = A : A = [A]
    A = [vec(i) for i in A]

    if length(A) == 1                   # if A has only 1 element
        return fill(A[1], n)            # repeat the elements n times
    elseif length(A) == n               # if A has n elements
        return A                        # return the original A vector
    else                                # otherwise: sho an error message
        error("Length of input must be 1 or $n, but got length $(length(A_vec))")
    end

end

function mpcDepicter(Set1,Set2,Crd1,Crd2)
    numSet1 = length(Set1)
    numSet2 = length(Set2)

    if numSet1 == 1
        Set1 = repeat(Set1,numSet2,1)
        Crd1 = repeat(Crd1,numSet2,1)
    elseif numSet2 == 1
        Set2 = repeat(Set2,numSet1,1)
        Crd2 = repeat(Crd2,numSet1,1)
    end

    n = length(Set1)

    points = gmsh.model.getEntities(0)
    point_tags = [point[2] for point in points] 
    maxptag = isempty(point_tags) ? 0 : maximum(point_tags)

    curves = gmsh.model.getEntities(1)
    curve_tags = [curve[2] for curve in curves] 
    maxctag = isempty(curve_tags) ? 0 : maximum(curve_tags)

    _, elementTags, _ = gmsh.model.mesh.getElements()
    all_tags = reduce(vcat, elementTags)
    maxetag = isempty(all_tags) ? 0 : maximum(all_tags)

    for i = 1:n
        x1 = Crd1[i,1]
        y1 = Crd1[i,2]
        z1 = Crd1[i,3]
        x2 = Crd2[i,1]
        y2 = Crd2[i,2]
        z2 = Crd2[i,3]
        n1 = Set1[i]
        n2 = Set2[i]

        p1 = gmsh.model.geo.addPoint(x1, y1, z1, 1.0, maxptag+2*i-1)
        p2 = gmsh.model.geo.addPoint(x2, y2, z2, 1.0, maxptag+2*i  ) 
        l1 = gmsh.model.geo.addLine(p1, p2, maxctag+i)
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.addElements(1, l1, [1], [[maxetag+i]], [[n1,n2]])
    end
end

function manual_sparse_submatrix(K::SparseMatrixCSC{Float64,Int}, Nset::Vector{Int})
    map_idx = Dict(n => i for (i, n) in enumerate(Nset))  # original -> new index
    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for (new_col, orig_col) in enumerate(Nset)
        col_start = K.colptr[orig_col]
        col_end = K.colptr[orig_col + 1] - 1
        for idx in col_start:col_end
            orig_row = K.rowval[idx]
            if haskey(map_idx, orig_row)
                push!(rows, map_idx[orig_row])
                push!(cols, new_col)
                push!(vals, K.nzval[idx])
            end
        end
    end

    n = length(Nset)
    return sparse(rows, cols, vals, n, n)
end

function rgb_to_int(r, g, b)
    return r * 256^2 + g * 256 + b
end

# A function to generating random spheres in a box
function randomSphereGen(a, b, c, Rmin, Rmax, t, d, h, n; attempt=10000000)
    
    box = [0, 0, 0, a, b, c]
    spheres = Matrix{Float64}(undef, 0, 5)  # initialization of an empty spheres matrix

    stop_spheregen = false
    
    for i in 1:n  # for each sphere
        valid = false    # set dimensions to invalid
        j = 0
        while !valid     # while the dimensions are invalid (-> regenerate dimensions)
            j += 1
            r = rand(Rmin:Rmax) # generate a random radius between Rmin and Rmax
             
            x = rand(r+h : a-r-h) # random x coordinate of the center point
            y = rand(r+h : b-r-h) # random y coordinate of the center point
            z = rand(r+h : c-r-h) # random z coordinate of the center point
            
            center = [x, y, z]          # coordinate vector of the center point
            
            # Check if the new sphere intersects with any other spheres:
            valid = true       # dimensions are valid (temporarily)

            for (x1,y1,z1,r1) in eachrow(spheres)  # for each existing spehres
                center1 = [x1, y1, z1]  # coordinate vector of an existing center point
                dist = norm(center - center1)  # distance between new and existing cetner point
                if dist < (r + r1 + d)  # if distance is smaller than allowable distance
                    valid = false       # dimensions are invalid  
                    break               # break the calculation -> generate new random values
                end
            end

            
            if valid  # if dimensions are valid
                spheres = [spheres;x y z r t]   # add new sphere to sphere matrix
            end

            if j >= attempt
                stop_spheregen = true
                break
            end
        end
        if stop_spheregen
            break  
        end
        println("sphere $i succeded after $j attempt")
    end
    
    return box, spheres
end

function getIntprop(id::Integer)

    elems = Dict{Int, Gauss}()

    elems[1]  = Gauss(:line,  :line2,   "CompositeGauss3")
    elems[2]  = Gauss(:tria,  :tria3,   "Gauss1")
    elems[3]  = Gauss(:quad,  :quad4,   "CompositeGauss3")
    elems[4]  = Gauss(:tetra, :tetra4,  "Gauss1")
    elems[5]  = Gauss(:hexa,  :hexa8,   "CompositeGauss3")
    elems[6]  = Gauss(:prism, :prism6,  "Gauss1")
    elems[7]  = Gauss(:pyra,  :pyra5,   "Gauss1")
    elems[8]  = Gauss(:line,  :line3,   "CompositeGauss5")
    elems[9]  = Gauss(:tria,  :tria6,   "Gauss2")
    elems[10] = Gauss(:quad,  :quad9,   "CompositeGauss5")
    elems[11] = Gauss(:tetra, :tetra10, "Gauss2")
    elems[12] = Gauss(:hexa,  :hexa27,  "CompositeGauss5")
    elems[13] = Gauss(:prism, :prism18, "Gauss4")
    elems[14] = Gauss(:pyra,  :pyra14,  "Gauss4")
    elems[15] = Gauss(:point, :point1,  "Gauss1")
    elems[16] = Gauss(:quad,  :quad8,   "CompositeGauss5")
    elems[17] = Gauss(:hexa,  :hexa20,  "CompositeGauss5")
    elems[18] = Gauss(:prism, :prism15, "Gauss4")
    elems[19] = Gauss(:pyra,  :pyra13,  "Gauss4")

    return elems[id]
end

function getIntprop(shape::Symbol,points::Integer)

    elems = Dict{Tuple, String}()

    elems[(:tria,1)]  = "Gauss1"
    elems[(:tria,3)]  = "Gauss2"
    elems[(:tria,4)]  = "Gauss3"
    elems[(:tria,6)]  = "Gauss4"
    elems[(:tria,7)]  = "Gauss5"
    elems[(:tria,12)] = "Gauss6"
    elems[(:quad,1)]  = "CompositeGauss1"
    elems[(:quad,3)]  = "Gauss1"
    elems[(:quad,4)]  = "CompositeGauss3"
    elems[(:quad,7)]  = "Gauss2"
    elems[(:quad,9)]  = "CompositeGauss5"
    elems[(:quad,16)] = "CompositeGauss7"
    elems[(:quad,25)] = "CompositeGauss9"

    return(elems[(shape,points)])

end

function getIntprop(sym::Symbol)
    str = string(sym)

    if startswith(str, "GL")
        return "CompositeGauss" * str[3:end]
    elseif startswith(str, "RG")
        return "Gauss" * str[3:end]
    else
        error("Unsupported symbol prefix: $sym")
    end
end

function dsearchn(xyzfull,xyzset)

    # Convert to Float64 and transpose:
    xyzfull = Matrix{Float64}(xyzfull')
    xyzset = Matrix{Float64}(xyzset')

    tree = KDTree(xyzfull) # KDTree

    k = 1  # nearest neighbor
    idxs, _ = knn(tree, xyzset, k)

    # Convert the result:
    nset0 = vec(idxs)
    nset1 = [nset0[i][1] for i in eachindex(nset0)]
    nset2 = sort(unique(nset1))

    return nset2
end

# A function to compute signed distances beyond box faces for each point
function dbox(xyz,boxvertices)
    
    n = length(boxvertices)

    x1 = boxvertices[1]
    y1 = boxvertices[2]
    z1 = n == 4 ? 0 : boxvertices[3]
    x2 = n == 4 ? boxvertices[3] : boxvertices[4]
    y2 = n == 4 ? boxvertices[4] : boxvertices[5]
    z2 = n == 4 ? 0 : boxvertices[6]

    d = maximum(hcat(
        y1 .- xyz[:, 2],
        xyz[:, 2] .- y2,
        x1 .- xyz[:, 1],
        xyz[:, 1] .- x2,
        z1 .- xyz[:, 3],
        xyz[:, 3] .- z2,
    ), dims=2)
    return vec(d)
end

function findrow(A, B)

    # finds the rows of B in A

    # Create a dictionary mapping rows to their indices
    l = Dict{Vector{Int}, Vector{Int}}()
    for i in 1:size(A,1)
        r = collect(A[i, :])  # ensure proper vector comparison
        if haskey(l, r)
            push!(l[r], i)
        else
            l[r] = [i]
        end
    end

    # Collect the matching indices, skip if not found
    result = Int[]
    for r in eachrow(B)
        rv = collect(r)
        if haskey(l, rv) && !isempty(l[rv])
            push!(result, popfirst!(l[rv]))
        end
    end
    return result
end


@inline function normedcrossprod(ax, ay, az, bx, by, bz)
    cx = ay*bz - az*by
    cy = az*bx - ax*bz
    cz = ax*by - ay*bx
    L = sqrt(cx^2 + cy^2 + cz^2)
    return cx / L, cy / L, cz / L
end

@inline function crossprod(ax, ay, az, bx, by, bz)
    cx = ay*bz - az*by
    cy = az*bx - ax*bz
    cz = ax*by - ay*bx
    return cx, cy, cz
end

@inline function csysgen(ax, ay, az, bx, by, bz)
    e3x, e3y, e3z = normedcrossprod(ax, ay, az, bx, by, bz)
    if e3y ≈ 1 || e3y ≈ -1
        e1x = 0
        e1y = 0
        e1z = 1
    else
        a1x = e3z
        a1y = 0
        a1z = -e3x
        L1 = sqrt(a1x^2+a1z^2)
        e1x = a1x/L1
        e1y = 0
        e1z = a1z/L1
    end
    e2x, e2y, e2z = crossprod(e3x, e3y, e3z, e1x, e1y, e1z)   
    return e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z 
end

@inline function invgen3(a11, a12, a13, a21, a22, a23, a31, a32, a33)
    # matrix determinant
    det = a11 * (a22 * a33 - a32 * a23) -
          a12 * (a21 * a33 - a31 * a23) +
          a13 * (a21 * a32 - a31 * a22)

    # elements of the inverse 
    b11 = (a22 * a33 - a32 * a23) / det
    b12 = (a32 * a13 - a12 * a33) / det
    b13 = (a12 * a23 - a22 * a13) / det
    b21 = (a31 * a23 - a21 * a33) / det
    b22 = (a11 * a33 - a31 * a13) / det
    b23 = (a21 * a13 - a11 * a23) / det
    b31 = (a21 * a32 - a31 * a22) / det
    b32 = (a31 * a12 - a11 * a32) / det
    b33 = (a11 * a22 - a21 * a12) / det

    return det, b11, b12, b13, b21, b22, b23, b31, b32, b33
end

# A function to convert Scope to vector of single scopes:
function scopeProcessor(scope0,multisel)
    if scope0 isa Selector && size(scope0.selpar,1) > 1 && multisel
        scope = Vector{Selector}(undef, size(scope0.selpar,1))
        for (i, row) in enumerate(eachrow(scope0.selpar))
            row = reshape(row, 1, :)
            scope[i] = Selector(scope0.object,scope0.seldef,row,scope0.presel,scope0.id,scope0.name,scope0.tol)
        end
    elseif scope0 isa Tuple && size(scope0[end],1) > 1 && multisel
        l = length(scope0)
        scope = Vector{Tuple}(undef, size(scope0[end],1))
            for (i, row) in enumerate(eachrow(scope0[end]))
            row = reshape(row, 1, :)
            if l == 2
                scope[i] = (scope0[1],row)
            elseif l == 3
                scope[i] = (scope0[1],scope0[2],row)
            end
        end
    else
        scope = to_vector(scope0)
    end
    return scope
end
    

Base.:+(a::BoundaryConditions, b::BoundaryConditions) = mergeBoundaryConditions([a, b])





