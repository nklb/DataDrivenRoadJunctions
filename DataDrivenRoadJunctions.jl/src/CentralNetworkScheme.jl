using LinearAlgebra

@enum BoundaryCondition ZeroFlux HomNeumann Outgoing Periodic FixedFluxOutgoing
@enum Scheme FirstOrder SecondOrder


mutable struct NumSolution11
    u₁::Vector{<:Real}
    u₂::Vector{<:Real}
    N::Integer
    t::Real
    xc::LinRange
end


mutable struct NumSolution21
    u₁::Vector{<:Real}
    u₂::Vector{<:Real}
    u₃::Vector{<:Real}
    N::Integer
    t::Real
    xc₁::LinRange
    xc₂::LinRange
    xc₃::LinRange
    oF::Real
end


mutable struct Problem21
    f₁
    f₂
    f₃
    λ::Real
    cpl
    bc::BoundaryCondition
    nodeslope::Bool
    Fe
    Fm
    fcut::Bool
end


function u(s::NumSolution11)
    [s.u₁; s.u₂]
end


function dx(s::NumSolution11)
    dx = (s.xc.stop - s.xc.start) / s.xc.lendiv
end


function constantinitialdata21(u_1::Real, u_2::Real, u_3::Real, dstretch::Real=1, N::Integer=201)
    u₁ = zeros(N-1)
    u₂ = zeros(N-1)
    u₃ = zeros(N-1)
    u₁ .= u_1
    u₂ .= u_2
    u₃ .= u_3
    t = 0
    xc₁ = LinRange(-dstretch, 0, N-1)
    xc₂ = LinRange(-dstretch, 0, N-1)
    xc₃ = LinRange(0, dstretch, N-1)
    NumSolution21(u₁, u₂, u₃, N, t, xc₁, xc₂, xc₃, 0)
end


function Problem21(f₁, f₂,  f₃, λ::Real, cpl, bc::BoundaryCondition)
    Problem21(f₁, f₂, f₃, λ, cpl, bc, 0.5, false, 0, 0)
end


function numflux(u_l::Real, u_r::Real, v_l::Real, v_r::Real, s_l⁺::Real, s_r⁻::Real,
                 λ::Real)
    0.5 * (v_r + v_l) - 0.5 * λ * (u_r - u_l) - .5 * (s_r⁻ - s_l⁺)
end

# modified for external coupling models
function coupling(p::Problem21, u_tr1::Real, u_tr2::Real, u_tr3::Real)
    out = vec(flux(p.cpl, reshape([u_tr1, u_tr2, u_tr3], 3, 1)))
    v_R1 = out[1]
    v_R2 = out[2]
    v_L3 = out[3]
    u_R1 = 0
    u_R2 = 0
    u_L3 = 0
    u_R1, u_R2, u_L3, v_R1, v_R2, v_L3
end


function mcslopes(v::Vector{<:Real})
    dv = 2 * diff(v)
    central = (v[3:end] - v[1:end-2]) / 2
    sdv = sign.(dv)

    (sdv[2:end] .== sdv[1:end-1]) .* sdv[2:end] .*
        min.(abs.(dv[1:end-1]), abs.(central), abs.(dv[2:end]))
end


function charslopes(u::Vector{<:Real}, v::Vector{<:Real}, λ::Real)
    mcslopes(.5 * (v - λ .* u)), mcslopes(.5 * (v + λ .* u))
end


function centralcoupling(p::Problem21, s::NumSolution21, T::Real, sc::Scheme=FirstOrder, CFL::Real = 0.9)
    dx = (s.xc₁.stop - s.xc₁.start) / s.xc₁.lendiv
    dtref = CFL * dx / p.λ
    N = s.N
    warned = false
    nfx₁ = zeros(N)
    nfx₂ = zeros(N)
    nfx₃ = zeros(N)

    F₁ = zeros(N-1)
    F₂ = zeros(N-1)
    F₃ = zeros(N-1)
    
    s₁⁻ = zeros(N-1)
    s₁⁺ = zeros(N-1)
    s₂⁻ = zeros(N-1)
    s₂⁺ = zeros(N-1)
    s₃⁻ = zeros(N-1)
    s₃⁺ = zeros(N-1)

    while s.t<T
        
        u_R1, u_R2, u_L3, v_R1, v_R2, v_L3 = coupling(p, s.u₁[end], s.u₂[end], s.u₃[1])

        if p.fcut
            v_R1 = positivepart(v_R1)
            v_R2 = positivepart(v_R2)
            v_L3 = positivepart(v_L3)
        end
        
        λ = p.λ

        for k=2:N-1
            λ = max(λ, abs(p.f₁(s.u₁[k]) - p.f₁(s.u₁[k-1])),
                         abs(p.f₂(s.u₂[k]) - p.f₂(s.u₂[k-1])),
                         abs(p.f₃(s.u₃[k]) - p.f₃(s.u₃[k-1])))
        end
        λ = max(λ, abs(p.f₃(s.u₃[1]) - p.f₁(s.u₁[N-1])),
                abs(p.f₃(s.u₃[1]) - p.f₂(s.u₂[N-1])))

        λ = max(λ, 2 * abs(v_R1 - p.f₁(s.u₁[N-1])),
                2 * abs(v_R2 - p.f₂(s.u₂[N-1])),
                2 * abs(p.f₃(s.u₃[1]) - v_L3))
        
        if p.bc == FixedFluxOutgoing
            λ = max(λ, abs(p.f₁(s.u₁[1]) - p.Fe(s.t)),
                    abs(p.f₂(s.u₂[1])) - p.Fm(s.t))
        end
        
        if sc == SecondOrder
            s₁⁻, s₁⁺ = charslopes(vcat(s.u₁[1], s.u₁, u_R1), vcat(p.f₁(s.u₁[1]), p.f₁.(s.u₁), v_R1), λ)
            s₂⁻, s₂⁺ = charslopes(vcat(s.u₂[1], s.u₂, u_R2), vcat(p.f₂(s.u₂[1]), p.f₂.(s.u₂), v_R2), λ)
            s₃⁻, s₃⁺ = charslopes(vcat(u_L3, s.u₃, s.u₃[end]), vcat(v_L3, p.f₃.(s.u₃), p.f₃(s.u₃[end])), λ)
        end
        
        #@show s₂⁺[end], s₁⁺[end], s₃⁻[1]
        if p.bc == HomNeumann
            nfx₁[1] = numflux(s.u₁[1], s.u₁[1], p.f₁(s.u₁[1]), p.f₁(s.u₁[1]), 0, s₁⁻[1], λ)
            nfx₂[1] = numflux(s.u₂[1], s.u₂[1], p.f₂(s.u₂[1]), p.f₂(s.u₂[1]), 0, s₂⁻[1], λ)
        elseif p.bc == FixedFluxOutgoing
            # @show nfx₁[1] = p.Fe(s.t)
            # nfx₂[1] = p.Fm(s.t)
            nfx₁[1] = numflux(s.u₁[1], s.u₁[1], p.Fe(s.t), p.f₁(s.u₁[1]), 0, s₁⁻[1], λ)
            nfx₂[1] = numflux(s.u₂[1], s.u₂[1], p.Fm(s.t), p.f₂(s.u₂[1]), 0, s₂⁻[1], λ)
        end
            
        if p.bc in [HomNeumann, Outgoing, FixedFluxOutgoing]
            nfx₃[N] = numflux(s.u₃[end], s.u₃[end], p.f₃(s.u₃[end]), p.f₃(s.u₃[end]), s₃⁺[end], 0, λ)
        end

        if !p.nodeslope
            s₁⁻[end]= 0
            s₂⁻[end]= 0
            s₃⁺[1] = 0
        end

        # fluxes away from the junction
        for k=2:N-1
            nfx₁[k] = numflux(s.u₁[k-1], s.u₁[k], p.f₁(s.u₁[k-1]), p.f₁(s.u₁[k]), s₁⁺[k-1], s₁⁻[k], λ)
            nfx₂[k] = numflux(s.u₂[k-1], s.u₂[k], p.f₂(s.u₂[k-1]), p.f₂(s.u₂[k]), s₂⁺[k-1], s₂⁻[k], λ)
            nfx₃[k] = numflux(s.u₃[k-1], s.u₃[k], p.f₃(s.u₃[k-1]), p.f₃(s.u₃[k]), s₃⁺[k-1], s₃⁻[k], λ)
        end

        # fluxes junction

        nfx₁[N] = v_R1
        nfx₂[N] = v_R2
        nfx₃[1] = v_L3

        #@show λ

        dt = CFL * dx / (2 * λ)
        
        # println("interface fluxes: street 1: ", nfx₁[[N-1, N]], ", street 2: ",
        #         nfx₂[[N-1, N]], ", street 3: ", nfx₃[[1, 2]])
        if abs(nfx₁[N] + nfx₂[N] - nfx₃[1]) > 1e-12 && ~warned
            println("Scheme not conservative")
            warned = true
        end
        
        # state update
        s.u₁ .-= dt/dx * diff(nfx₁)
        s.u₂ .-= dt/dx * diff(nfx₂)
        s.u₃ .-= dt/dx * diff(nfx₃)
        s.t += dt
        println("t=", s.t)
        
        if p.bc == FixedFluxOutgoing
            s.oF = p.f₃(s.u₃[end])
        end
    end
    s
end

