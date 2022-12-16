using DataDrivenRoadJunctions

dshash = "53e0af7"

dset_selection = string(gitroot(), "/out/dset_selection_", dshash, ".jls")
series = deserialize(dset_selection)["series"]
train_ids = deserialize(dset_selection)["train_ids"]
test_ids = deserialize(dset_selection)["test_ids"]
applic_ids = deserialize(dset_selection)["applic_ids"]

D_train = idstodata(train_ids, series)
D_test = idstodata(test_ids, series)
D_applic = idstodata(applic_ids, series)
M = greenshieldestimate(D_train)

ofile = string(gitroot(), "/out/data_", dshash, ".bson") 
bson(ofile, FD = M, training = D_train, test = D_test, application = D_applic)


