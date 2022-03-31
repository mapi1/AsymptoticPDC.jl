using PDC, Test

# Many test build on results from the original Matlab toolbox which serves as the gold standard here

# Prepare some data
y = get_sunspot_melanoma_data()
melanoma, _ = detrend(y[:, 3], 1)
sunspots, _ = detrend(y[:, 4], 1)
u = [sunspots melanoma]

for file in
    sort([file for file in readdir(@__DIR__) if match(r"^test_.*\.jl$", file) !== nothing])
    include(file)
end