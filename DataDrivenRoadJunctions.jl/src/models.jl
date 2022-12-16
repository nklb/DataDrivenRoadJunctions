using BlackBoxOptim
using Printf
using Interpolations
using MATLAB
using StatsBase
using EasyGit
using NNlib: softplus

# function softplus(x::Real)
#     log(exp(x) + 1)
# end


function rms(x)
    sqrt(mean(x.*x))
end


function rmssp(x::Vector{<:Real})
    y = softplus.(x)
    rms(y)
end


function tomatrixrows(input::Vector{Vector{Float64}})
    permutedims(reduce(hcat, input))
end


struct GreenshieldModel
    rhomax::Real
    vmax::Real
end


function greenshieldestimate(rho::Vector{<:Real}, vel::Vector{<:Real})
    filter = (vel .!= 0)
    # A = [ones(sum(filter)) rho[filter]]
    # leastsq = (transpose(A) * A)\(transpose(A) * vel[filter])
    
    # vmax = leastsq[1]
    # rhomax = -1/(leastsq[2]/vmax)



    function vdiff(x::Vector)
        v = x[1] .* (1 .- 1/x[2] .* rho[filter])
        return norm((v - vel[filter])) #rho[filter] .*
    end
    gloptimum = bboptimize(vdiff; SearchRange = [(0, 45), (0, 0.4)],
                           MaxFuncEvals = 2_000)
    vmax = best_candidate(gloptimum)[1]
    rhomax = best_candidate(gloptimum)[2]
    GreenshieldModel(rhomax, vmax)
end


function greenshieldestimate(D::MacroData, shifted::Bool=true)
    M = GreenshieldModel[]
  
    if shifted
        ρ = shiftedrho(D)
        v = shiftedv(D)
    else
        ρ = D.ρ
        v = D.v
    end

    for k=1:3
        push!(M, greenshieldestimate(ρ[k, :], v[k, :]))
    end
    M
end


function greenshieldestimate(D::Vector{MacroData}, shifted::Bool=true)
    M = GreenshieldModel[]
  
    if shifted
        ρ = shiftedrho(D)
        v = shiftedv(D)
    else
        error("to be implemented")
    end

    for k=1:3
        push!(M, greenshieldestimate(ρ[k, :], v[k, :]))
    end
    M
end


function velocity(M::GreenshieldModel, rho)
    M.vmax .* (1 .- rho ./ M.rhomax)
end


function flux(M::GreenshieldModel, rho)
    rho .* velocity(M, rho)
end


function cutoffflux(M::GreenshieldModel, rho)
    if rho > M.rhomax
        return 0
    else
        return flux(M, rho)
    end
end


function sigma(M::GreenshieldModel)
    0.5 * M.rhomax
end


function demand(M::GreenshieldModel, ρ)
    σ = sigma(M)
    maxdemand = flux(M, σ)

    [ρ[j] <= σ ? flux(M, min(ρ[j], M.rhomax)) : maxdemand for j=1:length(ρ)]
end


function supply(M::GreenshieldModel, ρ)
    σ = sigma(M)
    maxsupply = flux(M, σ)

    [ρ[j] <= σ ? maxsupply : flux(M, min(ρ[j], M.rhomax)) for j=1:length(ρ)] 
end


function inversefluxupper(M::GreenshieldModel, f)
    0.5 .* sqrt(M.rhomax / M.vmax) .* sqrt.(M.rhomax * M.vmax .- 4 .* f) .+ sigma(M)
end


function inversefluxlower(M::GreenshieldModel, f)
    -0.5 .* sqrt(M.rhomax / M.vmax) .* sqrt.(M.rhomax * M.vmax .- 4 .* f) .+ sigma(M)
end

    
function firstordercouplingETO(ρ::Matrix{<:Real}, M::Vector{GreenshieldModel},
                               β::Real, verbose::Bool=false)
    @assert size(ρ, 1) == 3
    @assert length(M) == 3

    N = size(ρ, 2)
    ρ_out = similar(ρ)
    fluxE = zeros(N)
    fluxT = zeros(N)
    fluxO = zeros(N)
    
    demandE = demand(M[1], ρ[1, :])
    demandT = demand(M[2], ρ[2, :])
    supplyO = supply(M[3], ρ[3, :])

    lowDemand = (demandE + demandT) .<= supplyO 
    #fluxO = min.(supplyO, demandE + demandT)
    fluxO[lowDemand] = demandE[lowDemand] + demandT[lowDemand]
    fluxE[lowDemand] = demandE[lowDemand]
    fluxT[lowDemand] = demandT[lowDemand]

    fluxO[.!lowDemand] = supplyO[.!lowDemand]
    
    Elow = [lowDemand[j] ? false : demandE[j] < β .* fluxO[j] for j=1:N]
    Tlow = [lowDemand[j] ? false : demandT[j] < (1-β) .* fluxO[j] for j=1:N]
    @assert !any(Elow .& Tlow)
    fluxE[Elow] = demandE[Elow]
    fluxT[Tlow] = demandT[Tlow]
    fluxE[Tlow] = (fluxO[Tlow] .- demandT[Tlow])
    fluxT[Elow] = (fluxO[Elow] .- demandE[Elow])
    row = .!lowDemand .& .!Elow .& .!Tlow   
    fluxE[row] = β .* fluxO[row]
    fluxT[row] = (1-β) .* fluxO[row]
    
    ρ_out[1, lowDemand .| Elow] = ρ[1, lowDemand .| Elow]
    ρ_out[2, lowDemand .| Tlow] = ρ[2, lowDemand .| Tlow]
    ρ_out[1, .!(lowDemand .| Elow)] =
        inversefluxupper(M[1], fluxE[.!(lowDemand .| Elow)]) 
    ρ_out[2, .!(lowDemand .| Tlow)] =
        inversefluxupper(M[2], fluxT[.!(lowDemand .| Tlow)]) 

    
    lowDemand = (demandE + demandT) .< supplyO 
    ρ_out[3, lowDemand] = inversefluxlower(M[3], fluxO[lowDemand])
    ρ_out[3, .!lowDemand] = fluxO[.!lowDemand]

    if verbose
        @printf "From %d inputs %d have low demand " N sum(lowDemand)
        @printf "and %d have low supply " sum(.!lowDemand)
        @printf "(%d limited by demand in E, \n" sum(Elow)
        @printf " %d limited by demand in T, " sum(Tlow)
        @printf " %d flux adjustments by right of way rule)\n" sum(row)
    end
    ρ_out, [transpose(fluxE); transpose(fluxT); transpose(fluxO)], [Elow, Tlow, row]
end


function ETOerror(D, ρ_out::Matrix{<:Real}, f_out::Matrix{<:Real},
                  shifted::Bool=true)

    if shifted
        rho_in = shiftedrho(D)
        f_in = shiftedflux(D)
        #t = referencetimerange(D)
    else
        rho_in = D.ρ
        f_in = flux(D)
        t = D.t
    end
    
    ρ_error = [rmsd(ρ_out[j, :], rho_in[j, :]) for j=1:3]
    f_error = [rmsd(f_out[j, :], f_in[j, :]) for j=1:3]


    push!(ρ_error, rmsd(vec(ρ_out), vec(rho_in))),
    push!(f_error, rmsd(vec(f_out), vec(f_in)))
end


function firstorderETOerror(D, M::Vector{GreenshieldModel}, β::Real,
                            shifted::Bool=true)

    if shifted
        ρ = shiftedrho(D)
    else
        ρ = D.ρ
    end
    
    ρ_out, f_out, _ = firstordercouplingETO(ρ, M, β)
    ETOerror(D, ρ_out, f_out, shifted)
end


function pressureAR(M::GreenshieldModel, rho, w)
    w ./ M.rhomax .* rho
end


function velocityAR(M::GreenshieldModel, rho, w)
    w .- pressureAR(M, rho, w)
end


function fluxAR(M::GreenshieldModel, rho, w)
    rho .* velocityAR(M, rho, w)
end

struct FlowMaxRS
    M::Vector{GreenshieldModel}
    β::Real
end


struct LinearRS
    M::Vector{GreenshieldModel}
    A::Matrix{<:Real}
end


function flux(RS::FlowMaxRS, ρ_in::Matrix{<:Real})
    _, f_out, _ = firstordercouplingETO(ρ_in, RS.M, RS.β)
    f_out
end


function rho(RS::FlowMaxRS, ρ_in::Matrix{<:Real})
    ρ_out, _, _ = firstordercouplingETO(ρ_in, RS.M, RS.β)
    ρ_out
end


function flux(RS, D)
    flux(RS, shiftedrho(D))
end


function rho(RS, D)
    flux(RS, shiftedrho(D))
end


function LinearRS(M::Vector{GreenshieldModel}, x::Vector{<:Real})
    B = reshape(x, 2, 6)
    A = [B[1,:]'; (B[2,:]-B[1,:])'; B[2, :]']
    LinearRS(M, A)
end


function flux(FM::LinearRS, ρ_in::Matrix{<:Real})
    f_in = [flux(FM.M[k], ρ_in[k, :]) for k=1:3]
    FM.A *  [ρ_in; permutedims(reduce(hcat, f_in))]
end


function fluxerror(FM, D)
    f_in = shiftedflux(D)
    f_out = flux(FM, D)

    f_error = [rmsd(f_out[j, :], f_in[j, :]) for j=1:3]

    push!(f_error, rmsd(vec(f_out), vec(f_in)))
end


function rhoerror(FM, ρ_in::Matrix{<:Real})
    ρ_out = rho(FM, ρ_in)
    
    ρ_error = [rmsd(ρ_out[j, :], ρ_in[j, :]) for j=1:3]

    push!(ρ_error, rmsd(vec(ρ_out), vec(ρ_in)))
end


function rhoerror(RS, D)
    rhoerror(RS, shiftedrho(D))
end


function shiftedflux(D::Vector{MacroData})
    data = Array{Float64}(undef, 3, 0)
    for d in D
        data = hcat(data, shiftedflux(d))
    end
    data
end


function shiftedrho(D::Vector{MacroData})
    data = Array{Float64}(undef, 3, 0)
    for d in D
        data = hcat(data, shiftedrho(d))
    end
    data
end


function shiftedv(D::Vector{MacroData})
    data = Array{Float64}(undef, 3, 0)
    for d in D
        data = hcat(data, shiftedv(d))
    end
    data
end


function rawshifteddata(D::Vector{MacroData})
    [shiftedrho(D); shiftedflux(D)]
end


function varyFD(M::Vector{GreenshieldModel}, x::Vector{<:Real})
    if length(x) == 3
        return [GreenshieldModel(M[k].rhomax, x[k] * M[k].vmax) for k = 1:3]
    elseif length(x) == 6
        return [GreenshieldModel(x[3+k] * M[k].rhomax, x[k] * M[k].vmax) for k = 1:3]
    else
        error("Input vector must have length 3 or 6.")
    end
end


function fdpenalty(RS, D)
    ρ = shiftedrho(D)
    f = shiftedflux(D)
    sum([rmsd(flux(RS.M[k], ρ[k, :]), f[k, :]) for k=1:3])
end


function negpenalty(RS, ρ::Matrix{<:Real}, reg::Bool)
    outflux = flux(RS, ρ)

    if !reg
        return rms(positivepart(-vec(outflux)))
    else
        return rmssp(-vec(outflux))
    end
end


function dspenalty(RS, ρ::Matrix{<:Real}, reg::Bool)
    outflux = flux(RS, ρ)
    dsdiff = vcat(outflux[1, :] - demand(RS.M[1], ρ[1, :]),
                  outflux[2, :] - demand(RS.M[2], ρ[2, :]),
                  outflux[3, :] - supply(RS.M[3], ρ[3, :]))

    if !reg
        return rms(positivepart(dsdiff))
    else
        return rmssp(dsdiff)
    end
end
