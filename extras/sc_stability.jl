# Stability of soft joints on a hanging particle chain, Gauss-Seidel vs a direct joint solve.
# Particle 1 hangs from the ground and each other particle from the one above, all in 1D.
# Without gravity and rest lengths one step is a linear map of positions, velocities, and impulses,
# so its spectral radius decides stability. There is no rotation, so no geometric stiffness.

using LinearAlgebra
using Printf

struct Chain
    invMasses::Vector{Float64}
    J::Matrix{Float64}
    D::Matrix{Float64}
    km::Vector{Float64}
end

function Chain(n, linkMass, endMass)
    invMasses = fill(1.0 / linkMass, n)
    invMasses[n] = 1.0 / endMass
    J = zeros(n, n)
    J[1, 1] = 1.0
    for j = 2:n
        J[j, j-1] = -1.0
        J[j, j] = 1.0
    end
    D = J * Diagonal(invMasses) * J'
    Chain(invMasses, J, D, diag(D))
end

struct Soft
    bias::Float64
    massScale::Float64
    impulseScale::Float64
end

# Infinite hertz is rigid with full Baumgarte
function soft_coefficients(hertz, zeta, h)
    if isinf(hertz)
        return Soft(1.0 / h, 1.0, 0.0)
    end
    omega = 2.0 * pi * hertz
    c = h * omega * (2.0 * zeta + h * omega)
    Soft(omega / (2.0 * zeta + h * omega), c / (1.0 + c), 1.0 / (1.0 + c))
end

# solver: :gs or :direct
# relax: :none, :gs (one rigid sweep), :exact (D dλ = -J v), :gamma (Γ kept on the left only)
struct Variant
    name::String
    solver::Symbol
    relax::Symbol
end

struct Config
    dt::Float64
    substeps::Int
    iterations::Int
    zeta::Float64
end

function step!(x, chain::Chain, variant::Variant, soft::Soft, config::Config, biasedFactor, relaxFactor)
    n = length(chain.km)
    y = view(x, 1:n)
    v = view(x, n+1:2n)
    lambda = view(x, 2n+1:3n)
    J = chain.J
    invM = chain.invMasses
    km = chain.km
    h = config.dt / config.substeps

    # Soft Step's leak as a compliance, Γ = K / c
    gamma = km .* (soft.impulseScale / soft.massScale)

    for s = 1:config.substeps
        v .+= invM .* (J' * lambda)

        if variant.solver == :gs
            for iteration = 1:config.iterations
                for j = 1:n
                    Jj = view(J, j, :)
                    cdot = dot(Jj, v)
                    C = dot(Jj, y)
                    dl = -soft.massScale / km[j] * (cdot + soft.bias * C) - soft.impulseScale * lambda[j]
                    lambda[j] += dl
                    v .+= invM .* Jj .* dl
                end
            end
        else
            dl = biasedFactor \ (-(J * v .+ soft.bias .* (J * y)) .- gamma .* lambda)
            lambda .+= dl
            v .+= invM .* (J' * dl)
        end

        y .+= h .* v

        if variant.relax == :gs
            for j = 1:n
                Jj = view(J, j, :)
                dl = -dot(Jj, v) / km[j]
                lambda[j] += dl
                v .+= invM .* Jj .* dl
            end
        elseif variant.relax == :exact || variant.relax == :gamma
            dl = relaxFactor \ (-(J * v))
            lambda .+= dl
            v .+= invM .* (J' * dl)
        end
    end
end

function spectral_radius(chain::Chain, variant::Variant, hertz, config::Config)
    h = config.dt / config.substeps
    soft = soft_coefficients(hertz, config.zeta, h)
    gamma = chain.km .* (soft.impulseScale / soft.massScale)
    biasedFactor = cholesky(Symmetric(chain.D + Diagonal(gamma)))
    relaxFactor = variant.relax == :gamma ? biasedFactor : cholesky(Symmetric(chain.D))

    n3 = 3 * length(chain.km)
    A = zeros(n3, n3)
    for k = 1:n3
        x = zeros(n3)
        x[k] = 1.0
        step!(x, chain, variant, soft, config, biasedFactor, relaxFactor)
        A[:, k] = x
    end

    if all(isfinite, A) == false
        return Inf
    end
    maximum(abs.(eigvals(A)))
end

# Lowest unstable hertz, Inf if stable everywhere up to rigid
function threshold(chain::Chain, variant::Variant, config::Config)
    stable(hz) = spectral_radius(chain, variant, hz, config) <= 1.0 + 1e-9
    grid = 10.0 .^ range(0, 5, length=101)
    lo = 0.0
    for hz in grid
        if stable(hz) == false
            hi = hz
            for k = 1:30
                mid = 0.5 * (lo + hi)
                stable(mid) ? (lo = mid) : (hi = mid)
            end
            return hi
        end
        lo = hz
    end
    stable(Inf) ? Inf : -1.0
end

fmt_hz(hz) = isinf(hz) ? "stable" : (hz < 0 ? "rigid unstable" : @sprintf("%.1f", hz))

gsRelax = Variant("gs relax", :gs, :gs)
gsNoRelax = Variant("gs", :gs, :none)
directExact = Variant("direct exact relax", :direct, :exact)
directGamma = Variant("direct gamma relax", :direct, :gamma)
directNoRelax = Variant("direct", :direct, :none)

# sc_gs.jl: n = 2, mass ratio 100, h = 1/60, 8 iterations, zeta 1
# Its simulation diverges between 45 and 48 Hz without relax, 50 and 52 Hz with relax
println("Validation against sc_gs.jl")
let chain = Chain(2, 1.0, 100.0), config = Config(1 / 60, 1, 8, 1.0)
    for variant in (gsNoRelax, gsRelax)
        @printf("  %-20s %s\n", variant.name, fmt_hz(threshold(chain, variant, config)))
    end
end

# Soft Step: 60 Hz step, 4 substeps, one biased sweep, joint zeta 10
softStep = Config(1 / 60, 4, 1, 10.0)
variants = (gsNoRelax, gsRelax, directNoRelax, directExact, directGamma)

println()
println("Lowest unstable joint hertz, Soft Step schedule (dt 1/60, 4 substeps, zeta 10)")
@printf("  %-6s %-6s", "links", "ratio")
for variant in variants
    @printf(" %20s", variant.name)
end
println()
for n in (2, 5, 10, 41, 160), ratio in (1.0, 100.0, 670.0)
    chain = Chain(n, 1.0, ratio)
    @printf("  %-6d %-6.0f", n, ratio)
    for variant in variants
        @printf(" %20s", fmt_hz(threshold(chain, variant, softStep)))
    end
    println()
end

println()
println("Spectral radius, Ball & Chain like (41 particles, end mass 670)")
let chain = Chain(41, 1.0, 670.0)
    @printf("  %-8s", "hertz")
    for variant in variants
        @printf(" %20s", variant.name)
    end
    println()
    for hz in (15.0, 30.0, 60.0, 120.0, 240.0, 480.0, 960.0, Inf)
        @printf("  %-8s", isinf(hz) ? "rigid" : @sprintf("%.0f", hz))
        for variant in variants
            @printf(" %20.6f", spectral_radius(chain, variant, hz, softStep))
        end
        println()
    end
end

println()
println("Lowest unstable joint hertz by substep count, 41 particles")
@printf("  %-9s %-6s %-6s %12s %12s %16s\n", "substeps", "ratio", "zeta", "gs", "gs relax", "gs relax / inv_h")
for substeps in (1, 2, 4, 8, 16), ratio in (1.0, 670.0), zeta in (1.0, 10.0)
    config = Config(1 / 60, substeps, 1, zeta)
    chain = Chain(41, 1.0, ratio)
    a = threshold(chain, gsNoRelax, config)
    b = threshold(chain, gsRelax, config)
    @printf("  %-9d %-6.0f %-6.0f %12s %12s %16.3f\n", substeps, ratio, zeta, fmt_hz(a), fmt_hz(b), b * config.dt / substeps)
end
