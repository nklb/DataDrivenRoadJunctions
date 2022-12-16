using Random
using DataDrivenRoadJunctions

rng = Xoshiro(98761);

series = "BB-global-incl-LaneIDs"

ids = getdsetids(series)

N = length(ids)
q = round(Integer, N/4)

shuffle!(rng, ids)

train_ids = sort(ids[1:q])
test_ids = sort(ids[q+1:2*q])
applic_ids = sort(ids[2*q+1:end])

ofile = string(gitroot(), "/out/dset_selection_", gitshorthead(), ".jls")
io = open(ofile, "w")
serialize(io, Dict([("train_ids", train_ids), ("test_ids", test_ids),
                    ("applic_ids", applic_ids), ("series", series)]))
close(io)
