using PDC, Test

@testset "Minimal" begin
    @test true
end

for file in
    sort([file for file in readdir(@__DIR__) if match(r"^test_.*\.jl$", file) !== nothing])
    include(file)
end