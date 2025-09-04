Szuper 🙌 Akkor itt egy **README.md**-stílusú dokumentáció a modulhoz, amit be tudsz tenni GitHubra vagy helyi jegyzetbe.

---

```markdown
# SparseSolveAuto.jl

**SparseSolveAuto.jl** egy kísérleti wrapper a Julia beépített `UmfpackLU` faktorizációjához, amely hatékonyan old meg **sparse jobb oldali** egyenletrendszereket (sparse vektor vagy mátrix esetén).  
A modul automatikusan választ **single-thread** és **multi-thread** módszer között, és a megoldás típusát is optimalizálja (dense vagy sparse visszatérés).

---

## Funkciók

- ✅ `solve_sparse_rhs(F, b::SparseVector)`  
  Sparse vektor jobb oldal megoldása.  
  A visszatérés `Vector` vagy `SparseVector`, a sparsity alapján.

- ✅ `solve_sparse_rhs(F, B::SparseMatrixCSC)`  
  Sparse mátrix jobb oldal megoldása.  
  Automatikusan eldönti, hogy single-thread vagy multi-thread ágat használjon.  
  A visszatérés `Matrix` vagy `SparseMatrixCSC`.

- ✅ **Automatikus kalibráció**  
  Az első híváskor lefut egy benchmark (`calibrate_rhs_threshold`), amely meghatározza,
  hány jobb oldal (`m`) esetén érdemes multi-thread módszerre váltani.  
  Az eredményt elmenti a fájlba:  
```

\~/.julia/umfpack\_threshold.toml

````

- ✅ **Perzisztens beállítások**  
A threshold és a sparsity-threshold értékek automatikusan betöltődnek új Julia futásnál.

- ✅ **Automatikus sparse kimenet**  
Ha a megoldás kitöltöttsége kisebb, mint `SPARSITY_THRESHOLD` (alapértelmezés: 0.1),
akkor `SparseVector` vagy `SparseMatrixCSC` formátumban tér vissza.

- ✅ **Információs logok**  
Minden futásnál látszik, hogy a solver single-thread vagy multi-thread ágat választott.

---

## Telepítés

Ez egy önálló modul.  
Másold be a `SparseSolveAuto.jl` fájlba, és töltsd be:

```julia
include("SparseSolveAuto.jl")
using .SparseSolveAuto
````

---

## Használat

```julia
using SparseArrays, LinearAlgebra
using .SparseSolveAuto

# Példa mátrix és jobb oldalak
A = sprand(2000, 2000, 0.001) + I
b = sprand(2000, 0.001)            # sparse vektor
C = sprand(2000, 50, 0.001)        # sparse mátrix (sok oszlop)

# Faktorizáció (UMFPACK)
F = lu(A)

# Sparse vektor jobb oldal
x = SparseSolveAuto.solve_sparse_rhs(F, b)
println("x típusa: ", typeof(x))

# Sparse mátrix jobb oldal
X = SparseSolveAuto.solve_sparse_rhs(F, C)
println("X típusa: ", typeof(X))
```

Első híváskor:

* Lefut a kalibráció (különböző `m` értékekre single vs multi mérések).
* A kapott küszöbérték (pl. `m ≥ 64`) elmentésre kerül.
* A további futások már a fájlból töltik be a beállítást.

---

## Konfiguráció

* **Threshold fájl helye**:
  `~/.julia/umfpack_threshold.toml`

* Például a fájl tartalma:

```toml
threshold = 64
sparsity_threshold = 0.1
```

* **Sparsity küszöb átállítása futás közben**:

```julia
SparseSolveAuto.SPARSITY_THRESHOLD[] = 0.05
```

---

## Korlátok

* Az `UmfpackLU` solver a Julia beépített SuiteSparse/UMFPACK implementációját használja.
* A jobb oldal (`b` vagy `B`) mindig teljesen megoldódik → a megoldás általában dense.
* Sparse visszatérés csak akkor várható, ha a megoldás valóban nagyon ritka.

---

## License

MIT (szabadon felhasználható és módosítható).

```

---

👉 Szeretnéd, hogy a README-hez készítsek egy kis **teljesítmény-összehasonlító ábrát** is (single vs multi), amit a `BenchmarkTools` + `Plots.jl`-lal lehet előállítani?
```




---

Nagyon jó 🙌 Akkor a README-be betehetünk egy kis teljesítmény-összehasonlító ábrát (single vs. multi), amit `BenchmarkTools` és `Plots.jl` segítségével lehet generálni. Ez jól illusztrálja, hogy a kalibráció miért fontos.

---

### Példa benchmark szkript (plots-hoz)

```julia
using SparseArrays, LinearAlgebra, Polyester, BenchmarkTools, Plots

function bench_rhs(n::Int, density::Float64, m_values::Vector{Int})
    A = sprand(n, n, density) + I
    F = lu(A)

    ts = Float64[]
    tm = Float64[]

    for m in m_values
        C = sprand(n, m, density)

        # single-thread mérés
        t_single = @belapsed begin
            n2, m2 = size(C)
            work = zeros(Float64, n2)
            for j in 1:m2
                fill!(work, 0.0)
                for k in C.colptr[j]:(C.colptr[j+1]-1)
                    work[C.rowval[k]] = C.nzval[k]
                end
                F \ work
            end
        end

        # multi-thread mérés
        t_multi = @belapsed begin
            n2, m2 = size(C)
            @batch for j in 1:m2
                work = zeros(Float64, n2)
                for k in C.colptr[j]:(C.colptr[j+1]-1)
                    work[C.rowval[k]] = C.nzval[k]
                end
                F \ work
            end
        end

        push!(ts, t_single)
        push!(tm, t_multi)

        println("m=$m → single=$(round(t_single*1000,digits=3)) ms, multi=$(round(t_multi*1000,digits=3)) ms")
    end

    return ts, tm
end

# Paraméterek
n = 2000
density = 0.001
m_values = [1, 2, 4, 8, 16, 32, 64, 128]

ts, tm = bench_rhs(n, density, m_values)

# Ábra rajzolás
plot(m_values, ts, label="single-thread", lw=2, marker=:o)
plot!(m_values, tm, label="multi-thread", lw=2, marker=:s)
xlabel!("RHS oszlopok száma (m)")
ylabel!("Idő [s]")
title!("Sparse RHS solve idő (n=$n, density=$density)")
```

---

### README-hez illeszthető példa ábra

Ha lefuttatod a fenti szkriptet, kapsz egy görbét:

* az `x` tengely: jobb oldali oszlopok száma (`m`)
* az `y` tengely: futási idő \[s]
* két görbe: single vs. multi

A képet mentsd le pl.:

```julia
savefig("rhs_benchmark.png")
```

Majd a README.md-be beillesztheted:

```markdown
## Teljesítmény példa

Az alábbi ábra a `solve_sparse_rhs` futási idejét mutatja különböző számú
jobb oldali oszlop (`m`) esetén, egy 2000×2000 méretű, 0.1% sűrűségű mátrixnál.

![RHS benchmark](rhs_benchmark.png)
```

---

👉 Szeretnéd, hogy a benchmark kódot kiegészítsem egy automatikus **küszöb-detektálással** is (pl. jelölje be az ábrán, hol vált gyorsabbá a multi a single)?
