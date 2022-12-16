using Printf
using DataDrivenRoadJunctions

series = "BB-global-incl-LaneIDs"
ids = getdsetids(series)
N = length(ids)

vehiclestats = [stats(series, getvehicledset(series, ids[k])) for k=1:N]

for k=1:N
    vstats = vehiclestats[k]
    @printf("%s & %d min, %d sec & %d & %f & %d & %f\\\\ \n", vstats[6],
            vstats[7], vstats[8],
            vstats[4], vstats[5], vstats[2], vstats[3])
end

AD = [MacroData(series, id) for id in ids]

for k=1:N
    dt = AD[k].t[2] - AD[k].t[1]
    @printf("dset %d: %.2f & %.2f & %.2f\n", k, dt * AD[k].shift[1], dt * AD[k].shift[2], dt * AD[k].shift[3])
end
