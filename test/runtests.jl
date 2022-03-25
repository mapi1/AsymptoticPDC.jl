using PDC, Test

# Many test build on results from the original Matlab toolbox which serves as the gold standard here

for file in
    sort([file for file in readdir(@__DIR__) if match(r"^test_.*\.jl$", file) !== nothing])
    include(file)
end