using Plots
using DataDrivenRoadJunctions    

series = "BB-global-incl-LaneIDs"

vehicles = getvehicledset(series, "10314")

trange, ein, min, out = inoutfluxes(series, vehicles)

plt = plot(trange, ein, label="influx on-ramp")
plot!(plt, trange, min, label="influx main lanes")
savefig(plt, string("../out/influxes_", gitshorthead(), ".pdf"))

plt2 = plot(trange, out, label="outflux")
savefig(plt2, string("../out/outfluxes_", gitshorthead(), ".pdf"))
