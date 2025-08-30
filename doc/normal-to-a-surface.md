using Gmsh
using LinearAlgebra

function get_element_normal(element_tag::Int, local_coords::Vector{Float64})
    
    # 1. Háló file beolvasása (a kódban feltételezzük, hogy ez már megtörtént)
    
    # 2. Lekérdezzük az elem típusát és a csomópontjait
    element_type, node_tags, _ = Gmsh.model.mesh.getElement(element_tag)
    
    if isempty(node_tags)
        @error "Az elem ($element_tag) nem létezik vagy nem tartalmaz csomópontokat."
        return nothing
    end
    
    # 3. Lekérdezzük a csomópontok 3D-s koordinátáit
    coord, _ = Gmsh.model.mesh.getNodes(node_tags)
    
    # A koordináták tömbjét 3D-s vektorokká alakítjuk
    num_nodes = length(node_tags)
    node_coords = [coord[i:i+2] for i in 1:3:length(coord)]
    
    # 4. Lekérdezzük a formafüggvények parciális deriváltjait
    # GMSH_2D_TRIANGLE, GMSH_2D_QUADRANGLE, stb.
    # a `basisType` határozza meg, milyen típusú bázisra van szükség.
    # A GMSH dokumentációja alapján a `getBasisFunctions` adja meg a deriváltakat.
    try
        basis_derivatives_uv = Gmsh.model.mesh.getBasisFunctions(element_type, "Grad", [0.0], local_coords)[1]
    catch e
        @error "Hiba a bázisfüggvények lekérdezésekor: $e"
        return nothing
    end

    # A `getBasisFunctions` egy lapos tömböt ad vissza (u1, v1, w1, u2, v2, w2, ...)
    # Szükségünk van az u és v szerinti deriváltakra (a w itt 0, mivel 2D-s elem)
    num_basis_funs = length(basis_derivatives_uv) ÷ 3
    dNi_du = basis_derivatives_uv[1:3:end]
    dNi_dv = basis_derivatives_uv[2:3:end]
    
    # 5. Kiszámoljuk a tangensvektorokat (Jacobian oszlopvektorai)
    # T_u = sum(dNi/du * Pi)
    # T_v = sum(dNi/dv * Pi)
    T_u = [0.0, 0.0, 0.0]
    T_v = [0.0, 0.0, 0.0]

    for i in 1:num_nodes
        T_u += dNi_du[i] * node_coords[i]
        T_v += dNi_dv[i] * node_coords[i]
    end

    # 6. Kiszámoljuk a normálvektort (kereszt-szorzat) és normalizáljuk
    normal_vector = cross(T_u, T_v)
    normalized_normal = normal_vector / norm(normal_vector)
    
    return normalized_normal
end

# Példa használat
Gmsh.initialize()
# Töltse be a 3D-s hálót, ami tartalmazza a görbült felületeket.
Gmsh.open("path/to/your/curved_mesh.msh")

# Cserélje ki a 123-at egy létező görbült felületi elem (pl. GMSH_TRIANGLE_6) azonosítójára.
element_tag = 123

# A pont a lokális koordinátákban (u, v)
# Például egy sarokcsomópont: (0, 0)
local_point = [0.0, 0.0]

normal = get_element_normal(element_tag, local_point)

if normal !== nothing
    println("A normálvektor a(z) $element_tag elemen a(z) $local_point pontban:")
    println(normal)
end

Gmsh.finalize()