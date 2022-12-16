using Plots
using LinearAlgebra
using Statistics
using StatsBase
using DelimitedFiles
using DataDrivenRoadJunctions


function positivepart(x::Real)
    max(x, 0)
end


function negativepart(x::Real)
    max(-x, 0)
end


series = "BB-global-incl-LaneIDs"

vehicles = getvehicledset(series, "10324")
t0, T, zs = observationtime(vehicles)

EntryVol, MainInVol, OutVol = getvolumes(series)
t = defaulttrange(vehicles)
M = length(t)
ρ, v = numberdensities([EntryVol, MainInVol, OutVol], vehicles, t)
fluxes = [ρ[j, :] .* v[j, :] for j=1:3]

dt = t[2] - t[1]
shiftrangeT=-20:20
shiftrangeO=-100:100
N = (length(shiftrangeT), length(shiftrangeO))
deviation = zeros(N)

for i=1:N[1]
    for j=1:N[2]
        s2 = shiftrangeT[i]
        s3 = shiftrangeO[j]
        xrange = (1 + max(negativepart(s2), negativepart(s3))):(M-max(positivepart(s2),
                                                                      positivepart(s3)))
        deviation[i, j] = rmsd(fluxes[1][xrange] + fluxes[2][xrange .+ s2], fluxes[3][xrange .+ s3])
    end
end

outfile = string(gitroot(), "/out/delay")
HM = heatmap(dt * shiftrangeO, dt * shiftrangeT, deviation, xlabel="τ³", ylabel="τ²")
savefig(HM, string(outfile, "_hmap.pdf"))

xx = [i for i in (dt * shiftrangeO), j in 1:length(shiftrangeT)]
yy = [j for i in 1:length(shiftrangeO), j in (dt*shiftrangeT)]

open(string(outfile, ".dat"), "w") do io
    writedlm(io, [vec(xx) vec(yy) vec(transpose(deviation))])
end
P = []

push!(P, heatmap(dt * shiftrangeO, dt * shiftrangeT, deviation, xlabel="τ_3", ylabel="τ_2"))

best = argmin(deviation)
s2 = shiftrangeT[best[1]]
s3 = shiftrangeO[best[2]]
xrange = (1 + max(negativepart(s2), negativepart(s3))):(M-max(positivepart(s2),
                                                              positivepart(s3)))
influx = fluxes[1][xrange] + fluxes[2][xrange .+ s2]
outflux = fluxes[3][xrange .+ s3]
err = influx-outflux

open(string(outfile, "_bestfit.dat"), "w") do io
    writedlm(io, [t[xrange] influx outflux err])
end

PF = plot(t[xrange], fluxes[1][xrange] + fluxes[2][xrange .+ s2], label="f_E(.)+f_T(.+τ2)")
plot!(PF, t[xrange], fluxes[3][xrange .+ s3], label="f_O(.+τ3)")
push!(P, PF)

push!(P, plot(t[xrange], err, label="conservation error"))
println("Lowest deviation (", minimum(deviation), ") for τ_2=", shiftrangeT[best[1]] * dt, " and τ_3=", shiftrangeO[best[2]] * dt, ".")

savefig(plot(P ..., layout=(3,1), size=(3200,1800)), string(gitroot(),"/flux_conservation_vs_delay_", gitshorthead(), ".pdf"))
