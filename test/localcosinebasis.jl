using WaveletsExt:LocalCosineBasis

data = rand(Float64, 1024)
fit = LocalCosineBasis.analysis(1024)
results = fit(data);
cost = LocalCosineBasis.cost(results);
@info "tree size:$(size(results)), cost size $(size(cost))"
@info cost[1:10]
