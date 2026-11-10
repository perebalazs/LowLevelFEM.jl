module SparseSolveAuto

using SparseArrays, LinearAlgebra, Polyester, BenchmarkTools, TOML

const CONFIG_FILE = joinpath(homedir(), ".julia", "umfpack_threshold.toml")
const RHS_PARALLEL_THRESHOLD = Ref{Int}(0)
const SPARSITY_THRESHOLD = Ref{Float64}(0.1)  # ha kitöltöttség < 10%, akkor sparse visszatérés

# -------------------------------
# Threshold kezelése
# -------------------------------
function save_threshold(thr::Int)
    mkpath(dirname(CONFIG_FILE))
    open(CONFIG_FILE, "w") do io
        TOML.print(io, Dict("threshold" => thr,
                            "sparsity_threshold" => SPARSITY_THRESHOLD[]))
    end
end

function load_threshold()
    if isfile(CONFIG_FILE)
        try
            data = TOML.parsefile(CONFIG_FILE)
            thr = get(data, "threshold", 0)
            SPARSITY_THRESHOLD[] = get(data, "sparsity_threshold", SPARSITY_THRESHOLD[])
            return thr
        catch
            return 0
        end
    end
    return 0
end

# -------------------------------
# Kalibráció
# -------------------------------
function calibrate_rhs_threshold(n::Int, density::Float64; m_values=[1,2,4,8,16,32,64,128])
    A = sprand(n, n, density) + I
    F = lu(A)

    best_threshold = typemax(Int)

    @info "Kalibráció indul (n=$n, density=$density)"
    for m in m_values
        C = sprand(n, m, density)

        t_single = @elapsed begin
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

        t_multi = @elapsed begin
            n2, m2 = size(C)
            @batch for j in 1:m2
                work = zeros(Float64, n2)
                for k in C.colptr[j]:(C.colptr[j+1]-1)
                    work[C.rowval[k]] = C.nzval[k]
                end
                F \ work
            end
        end

        @info "m=$m → single=$(round(t_single*1000,digits=3)) ms, multi=$(round(t_multi*1000,digits=3)) ms"

        if t_multi < t_single
            best_threshold = min(best_threshold, m)
        end
    end

    thr = best_threshold == typemax(Int) ? typemax(Int) : best_threshold
    RHS_PARALLEL_THRESHOLD[] = thr
    save_threshold(thr)

    if thr == typemax(Int)
        @info "→ Nincs olyan m, ahol multi gyorsabb. Marad single-thread."
    else
        @info "→ Ajánlott küszöb: m ≥ $thr → multi-thread"
    end

    return thr
end

# -------------------------------
# Threshold betöltése
# -------------------------------
function ensure_threshold!(n::Int=2000, density::Float64=0.001)
    if RHS_PARALLEL_THRESHOLD[] == 0
        thr = load_threshold()
        if thr != 0
            @info "Betöltött threshold: $thr (fájlból: $CONFIG_FILE)"
            RHS_PARALLEL_THRESHOLD[] = thr
        else
            @info ">> Első hívás: kalibráció indul..."
            calibrate_rhs_threshold(n, density)
        end
    end
    return RHS_PARALLEL_THRESHOLD[]
end

# -------------------------------
# Megoldók
# -------------------------------
function solve_sparse_rhs(F::SparseArrays.UMFPACK.UmfpackLU{Float64,Int}, b::SparseVector{Float64,Int})
    n = F.n
    work = zeros(Float64, n)
    for (i,v) in zip(b.nzind, b.nzval)
        work[i] = v
    end
    @info "solve_sparse_rhs: sparse vektor RHS → single-thread"
    x = F \ work
    return count(!iszero, x) / n < SPARSITY_THRESHOLD[] ? sparsevec(x) : x
end

function solve_sparse_rhs(F::SparseArrays.UMFPACK.UmfpackLU{Float64,Int}, B::SparseMatrixCSC{Float64,Int})
    n, m = size(B)
    thr = ensure_threshold!(n, nnz(B) / (n*m))

    X = Matrix{Float64}(undef, n, m)

    if m < thr
        @info "solve_sparse_rhs: sparse mátrix RHS → single-thread (m=$m, threshold=$thr)"
        work = zeros(Float64, n)
        for j in 1:m
            fill!(work, 0.0)
            for k in B.colptr[j]:(B.colptr[j+1]-1)
                work[B.rowval[k]] = B.nzval[k]
            end
            X[:, j] = F \ work
        end
    else
        @info "solve_sparse_rhs: sparse mátrix RHS → multi-thread (m=$m, threshold=$thr)"
        @batch for j in 1:m
            work = zeros(Float64, n)
            for k in B.colptr[j]:(B.colptr[j+1]-1)
                work[B.rowval[k]] = B.nzval[k]
            end
            X[:, j] = F \ work
        end
    end

    return count(!iszero, X) / length(X) < SPARSITY_THRESHOLD[] ? sparse(X) : X
end

end # module

