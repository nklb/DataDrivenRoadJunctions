using DataFrames
using Polynomials
using Serialization
using TranscodingStreams
using CodecXz
using Dates
using EasyGit

struct VehicleTemplate
    width::Real
    length::Real
    contour::DataFrame
end


struct Vehicle
    vehicle_id::Integer
    discrete_trajectory::DataFrame

    fitTypeE::String
    traj_EtGlob::Polynomial
    traj_tEGlob::Polynomial

    fitTypeN::String
    traj_NtGlob::Polynomial
    traj_tNGlob::Polynomial

    fitted_trajectoryGlob::Array{Real, 2}
    fitted_speedGlob::Array{Real, 2}
    fitted_accGlob::Array{Real, 2}

    template_id::Integer
end


function loadcompressed(filepath::String)
    io = open(filepath, "r")
    io = TranscodingStream(XzDecompressor(), io)
    content = deserialize(io)
    close(io)
    content
end


function getvehicledset(series::String, id::String)::Vector{Vehicle}
    loadcompressed(string(gitroot(), "/data/", series, "/", id, ".jls.xz"))
end


function getvehicletemplates(series::String)
    loadcompressed(string(gitroot(), "/data/", series, "/VehicleTemplates.jls.xz"))
end


function getdsetids(series::String)
    dsetfld = string(gitroot(), "/data/", series, "/")
    ids = String[]
    for f in readdir(dsetfld)
        length(f) < 7 && continue
        if f[end-6:end] == ".jls.xz"
            push!(ids, f[1:end-7])
        end
    end
    ids
end


function vehiclepositions(vehicles::Vector{Vehicle}, t_range::LinRange{Float64, Int64}, cut::Bool=true)
    N = length(vehicles)
    M = length(t_range)

    vx = zeros(Real, N, M)
    vy = zeros(Real, N, M)

    t_cut = ones(length(t_range))

    for i=1:N
        if cut
            t_cut = (t_range .>= vehicles[i].fitted_trajectoryGlob[1,3]) .& 
                (t_range .<= vehicles[i].fitted_trajectoryGlob[end,3])
        end
        vx[i, :] =  t_cut .* [vehicles[i].traj_EtGlob(t) for t in t_range]
        vy[i, :] =  t_cut .* [vehicles[i].traj_NtGlob(t) for t in t_range]
    end
    vx, vy
end


function vehiclespeed(vehicles::Vector{Vehicle}, t_range, cut::Bool=true)
    N = length(vehicles)
    M = length(t_range)

    vvx = zeros(Real, N, M)
    vvy = zeros(Real, N, M)

    t_cut = ones(length(t_range))

    for i=1:N
        if cut
            t_cut = (t_range .>= vehicles[i].fitted_trajectoryGlob[1,3]) .& 
                (t_range .<= vehicles[i].fitted_trajectoryGlob[end,3])
        end
        vvx[i, :] =  t_cut .* [derivative(vehicles[i].traj_EtGlob)(t)
                               for t in t_range]
        vvy[i, :] =  t_cut .* [derivative(vehicles[i].traj_NtGlob)(t)
                               for t in t_range]
    end
    vvx, vvy
end


function observationtime(vehicles::Vector{Vehicle})
    N = length(vehicles)
    
    t0 = minimum([vehicles[i].fitted_trajectoryGlob[1, 3] for i = 1:N])
    T = maximum([vehicles[i].fitted_trajectoryGlob[end, 3] for i = 1:N])
    stamp = vehicles[1].discrete_trajectory[1, 5]
    zeroDateTime = DateTime(string(stamp[1:10], "T", stamp[12:end-2])) -
        Microsecond(round(vehicles[1].discrete_trajectory[1,7] * 1e6))
    t0, T, zeroDateTime
end
