### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ b5cad96e-e21a-44e7-8726-942da6ee0140
using LowLevelFEM

# ╔═╡ 9cf25288-48b3-4abf-b01b-0c85f056b216
structured_box_mesh()

# ╔═╡ eefcf225-f3be-4be4-9160-e9343fa0ba7d
mat = Material("body")

# ╔═╡ b881cb08-92af-449e-ba3e-5a4f101d6efd
U = Field([mat], type=:VectorField, dim=3, fieldName=:u)

# ╔═╡ 99682979-dc69-4f54-a0a4-fa775c5199d4
u = VectorField(U, "body", [1,2,3])

# ╔═╡ 6bf9c0cc-61cb-4a32-a8b1-31ea0c5e4042
showElementResults(u)

# ╔═╡ 60ffffb4-bcf6-4946-8b88-ebb5e3d08081
openPostProcessor()

# ╔═╡ Cell order:
# ╠═b5cad96e-e21a-44e7-8726-942da6ee0140
# ╠═9cf25288-48b3-4abf-b01b-0c85f056b216
# ╠═eefcf225-f3be-4be4-9160-e9343fa0ba7d
# ╠═b881cb08-92af-449e-ba3e-5a4f101d6efd
# ╠═99682979-dc69-4f54-a0a4-fa775c5199d4
# ╠═6bf9c0cc-61cb-4a32-a8b1-31ea0c5e4042
# ╠═60ffffb4-bcf6-4946-8b88-ebb5e3d08081
