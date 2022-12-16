using DataDrivenRoadJunctions
using Plots
using DelimitedFiles

dshash = "53e0af7"

D = loadmacrodata(dshash, :training)
M = loadmacrodata(dshash, :FD)

ρ = 1000 .* shiftedrho(D)
v = shiftedv(D)
st = shiftedflux(D)

filter = (v .!= 0)

ofile = string(gitroot(), "/out/FD_", gitshorthead())

plots_p1 = []
for j=1:3
    GS = scatter(ρ[j, filter[j, :]], st[j, filter[j, :]], label = "data",
                 title = string("f=ρv on road ",j), legend=:topright,
                 xlabel="ρ", ylabel="f")
    x = LinRange(0, M[j].rhomax, 1000)
    GS = plot!(GS, 1000 * x, flux(M[j], x), label = "f")

    open(string(ofile, "_road",  j, ".dat"), "w") do io
        writedlm(io, [1000*x flux(M[j], x)])
    end

    open(string(ofile, "_data_road",  j, ".dat"), "w") do io
        writedlm(io, [ρ[j, filter[j, :]] st[j, filter[j, :]]])
    end
    push!(plots_p1, GS)
end

savefig(plot(plots_p1 ..., layout=(1,3), size=(3200, 600)),
        string(gitroot(), "/out/FD_", gitshorthead(), ".pdf"))	

