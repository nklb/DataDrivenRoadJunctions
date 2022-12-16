using MAT
using BSON
using LinearAlgebra
using PolygonOps
using StatsBase
using NLsolve
using KernelDensity
using Interpolations


struct MacroData
    ρ::Matrix{Float64}
    v::Matrix{Float64}
    t::LinRange{Float64, Int64}
    shift::Vector{Int64}
end


function getvolume(series::String, varname::String)::Vector{Vector{Real}}
    Lmat = matopen(string(gitroot(), "/FCD/", series, "/Volumes.mat"))
    Vol = read(Lmat, varname)
    [Vol[i, :] for i = [collect(1:size(Vol, 1)); 1]]
end


function getvolumes(series::String)
    EntryVol = getvolume(series, "APoly")
    MainInVol = getvolume(series, "BPoly")
    OutVol = getvolume(series, "CPoly")

    EntryVol, MainInVol, OutVol
end


function diameter(Vol::Vector{Vector{Real}})
    n = length(Vol)
    dist = zeros(n, n)
    for j=1:n
        for k=1:n
            dist[j, k] = norm(Vol[j] - Vol[k], 2)
        end
    end
    maximum(vec(dist))
end

    
function xinout(series::String)
    EntryVol, MainInVol, OutVol = getvolumes(series)
    
    xio = (maximum([MainInVol[k][1] for k=1:length(MainInVol)]),
           minimum([OutVol[k][1] for k=1:length(OutVol)]))
    yrange = (maximum([MainInVol[k][2] for k=1:length(MainInVol)]),
              minimum([OutVol[k][2] for k=1:length(OutVol)]))

    xio, sqrt( (xio[1]-xio[2])^2 + (yrange[1] - yrange[2])^2 )
end


function volumetoshape(Vol::Vector{Vector{Real}})
    N = size(Vol, 1)
    Shape([Vol[i][1] for i = 1:N], [Vol[i][2] for i=1:N])
end


function passingvehicles(Vol::Vector{Vector{Real}},
                         vehicles::Vector{Vehicle},
                         t_range::LinRange{Float64, Int64})
    N = length(vehicles)
    M = length(t_range)
    n = zeros(N)
    vx, vy = vehiclepositions(vehicles, t_range, false)

    for j = 1:N
        indicator = [inpolygon([vx[j, k], vy[j, k]], Vol) for k=1:M]
        n[j] = sum(indicator) / M
        indicator .= false
    end
    n
end


function passingvehicles(Vol::Vector{Vector{Real}},
                         vehicles::Vector{Vehicle})
    passingvehicles(Vol, vehicles, defaulttrange(vehicles))
end


function classifyvehicles(VolA::Vector{Vector{Real}}, VolB::Vector{Vector{Real}},
                          vehicles::Vector{Vehicle})

    na = passingvehicles(VolA, vehicles)
    nb = passingvehicles(VolB, vehicles)

    N = length(vehicles)
    pa = Array{Bool}(undef, N)
    pb = Array{Bool}(undef, N)
    pa .= false
    pb .= false
    
    pa[na .> nb] .= true
    pb[na .< nb] .= true

    pa, pb
end


function inouttime(V::Vehicle, x_in::Real, x_out::Real)
    function inoutdist!(F, t)
        F[1] = V.traj_EtGlob(t[1]) - x_in
        F[2] = V.traj_EtGlob(t[2]) - x_out
    end

    function j!(J, t)
        J[1, 1] = derivative(V.traj_EtGlob)(t[1])
        J[2, 2] = derivative(V.traj_EtGlob)(t[2])
        J[1, 2] = 0
        J[2, 1] = 0
    end
        
    y = nlsolve(inoutdist!, j!, [V.fitted_trajectoryGlob[1, 3],
                                 V.fitted_trajectoryGlob[end, 3]])
    y.zero
end


function disappearing(series::String, vehicles::Vector{Vehicle})
    EntryVol, MainInVol, OutVol = getvolumes(series)
    
    E, M = classifyvehicles(EntryVol, MainInVol, vehicles)
    O = passingvehicles(OutVol, vehicles)
    O .== 0
end


function inoutfluxes(series::String, vehicles::Vector{Vehicle})
    filter = .! disappearing(series, vehicles)
    xio, dist = xinout(series)
    iostretch = [inouttime(vehicles[k], xio[1], xio[2]) for k=1:length(vehicles)]

    tstart = minimum([iot[1] for iot in iostretch[filter]])
    tend = maximum([iot[2] for iot in iostretch[filter]]) + 5

    EntryVol, MainInVol, OutVol = getvolumes(series)
    E, M = classifyvehicles(EntryVol, MainInVol, vehicles)

    ne = sum(filter .&& E)
    nm = sum(filter .&& M)
    no = sum(filter)

    DE = kde([iot[1] for iot in iostretch[filter .&& E]], bandwidth=.75)
    DM = kde([iot[1] for iot in iostretch[filter .&& M]], bandwidth=.75)

    DO = kde([iot[2] for iot in iostretch[filter]], bandwidth=.75)

    trange = LinRange(tstart, tend, 5000)

    kinterp(K, nk, x) = nk * pdf(K::UnivariateKDE, x)
    trange, kinterp(DE, ne, trange), kinterp(DM, nm, trange), kinterp(DO, no, trange)
end
    

function stats(series::String, vehicles::Vector{Vehicle})
    EntryVol, MainInVol, OutVol = getvolumes(series)

    E, M = classifyvehicles(EntryVol, MainInVol, vehicles)
    O = passingvehicles(OutVol, vehicles)
    fltr = O .!= 0

    println("Sorted out ", sum(.!fltr), " trajectories.")
    evehicles = vehicles[E] 
    mvehicles = vehicles[M]
    xio, dist = xinout(series)

    iostretch = [inouttime(vehicles[k], xio[1], xio[2]) for k=1:length(vehicles)]
    meanspeed = [ 1e-3 * dist / ( 1/3600 * (iostretch[k][2] - iostretch[k][1]))
                  for k=1:length(vehicles)]

    println(sum(fltr .&& E), " entering cars with mean speed ", mean(meanspeed[fltr .&& E]))
    println(sum(fltr .&& M), " passing cars with mean speed ", mean(meanspeed[fltr .&& M]))

    T0, T, starttime = observationtime(vehicles)

    [sum(.!fltr), sum(fltr .&& E), mean(meanspeed[fltr .&& E]), sum(fltr .&& M),
     mean(meanspeed[fltr .&& M]), starttime, round(T/60), mod(T, 60)]
end


function numberdensities(Vols::Vector{Vector{Vector{Real}}},
                         vehicles::Vector{Vehicle},
                         t_range::LinRange{Float64, Int64})
    M = length(t_range)
    vn = length(Vols)

    ρ = zeros(vn, M)
    v = zeros(vn, M)
    vx, vy = vehiclepositions(vehicles, t_range)
    vvx, vvy = vehiclespeed(vehicles, t_range)

    diam = zeros(vn)

    for i=1:vn
        diam[i] = diameter(Vols[i])
    end    

    for i=1:M
        for j=1:vn
            indicator = [inpolygon([vx[k, i], vy[k, i]],
                                   Vols[j])
                         for k=1:size(vx,1)]
            ρ[j, i] = sum(indicator)
            v[j, i] = dot(indicator,
                          sqrt.(vvx[:, i] .* vvx[:, i] .+
                              vvy[:, i] .* vvy[:,i])) ./
                      max(ρ[j, i], 1)
        end
    end

    for j=1:vn
        ρ[j, :] = ρ[j, :] / diam[j] #abs(PolygonOps.area(Vols[j]))
    end

    ρ, v
end

      
function positivepart(x::Real)
    max(x, 0)
end


function negativepart(x::Real)
    max(-x, 0)
end


function positivepart(x::Vector{<:Real})
    positivepart.(x)
end


function negativepart(x::Vector{<:Real})
    negativepart.(x)
end


function estimateshift(ρ::Matrix{Float64}, v::Matrix{Float64}, t::LinRange{Float64, Int64})

    fluxes = [ρ[j, :] .* v[j, :] for j=1:3]
    
    dt = t[2] - t[1]
    boundT = round(Integer, 5/dt)
    boundO = round(Integer, 25/dt)
    shiftrangeT=-boundT:boundT
    shiftrangeO=-boundO:boundO
    N = (length(shiftrangeT), length(shiftrangeO))
    deviation = zeros(N)

    for i=1:N[1]
        for j=1:N[2]
            s2 = shiftrangeT[i]
            s3 = shiftrangeO[j]
            trange = (1 + max(negativepart(s2), negativepart(s3))):(t.len-max(positivepart(s2),
                                                                      positivepart(s3)))
            deviation[i, j] = rmsd(fluxes[1][trange] + fluxes[2][trange .+ s2],
                                   fluxes[3][trange .+ s3])
        end
    end

    best = argmin(deviation)

    [0, shiftrangeT[best[1]], shiftrangeO[best[2]]], minimum(deviation)
end


function defaulttrange(vehicles::Vector{Vehicle})
    t0, T, _ = observationtime(vehicles)
    M = ceil(Integer, (T - t0) * 4)
    LinRange(t0, T, M)
end


function loadmacrodata(dshash::String, dsclass)
    @assert dsclass in [:FD, :training, :test, :application]
    if dsclass == :FD
        return convert(Vector{GreenshieldModel}, BSON.load(string(gitroot(), "/out/data_", dshash, ".bson"))[:FD])
    else
        return convert(Vector{MacroData}, BSON.load(string(gitroot(), "/out/data_", dshash, ".bson"))[dsclass])
    end
end


function MacroData(series::String, dset::String)
    vehicles = getvehicledset(series, dset)
    EntryVol, MainInVol, OutVol = getvolumes(series)

    t = defaulttrange(vehicles)
    ρ, v = numberdensities([EntryVol, MainInVol, OutVol], vehicles, t)
    shift, _ = estimateshift(ρ, v, t)
    
    MacroData(ρ, v, t, shift)
end


function flux(D::MacroData)
    D.ρ .* D.v
end


function shiftedrange(D::MacroData)
    (1 + max(negativepart(D.shift[2]),
             negativepart(D.shift[3]))):(D.t.len-max(positivepart(D.shift[2]),
                                               positivepart(D.shift[3])))
end


function referencetimerange(D::MacroData)
    trange = shiftedrange(D)
    D.t[trange]
end


function shiftedrho(D::MacroData)
    trange = shiftedrange(D)
    [D.ρ[1, trange]'; D.ρ[2, trange .+ D.shift[2]]'; D.ρ[3, trange .+ D.shift[3]]']
end


function shiftedv(D::MacroData)
    trange = shiftedrange(D)
    [D.v[1, trange]'; D.v[2, trange .+ D.shift[2]]'; D.v[3, trange .+ D.shift[3]]']
end


function shiftedflux(D::MacroData)
    ρ = shiftedrho(D)
    v = shiftedv(D)
    ρ .* v
end


function delay(D::MacroData)
    D.shift .* (D.t[2] - D.t[1])
end


function idstodata(myids, series)
    [MacroData(series, id) for id in myids]
end

