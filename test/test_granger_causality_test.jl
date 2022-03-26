# Prepare some data
y = get_sunspot_melanoma_data()
melanoma, _ = detrend(y[:, 3], 1)
sunspots, _ = detrend(y[:, 4], 1)
u = [sunspots melanoma]

mattol = 0.001 # tolerance from matlab

@testset "Granger Causality Test" begin
    (_, pf, A, _, _, _) = mvar(u, maxorder=2, method="NS", criterion=nothing)
    connectivity, pvalues = granger_causality_test(u, A, pf)
      
    @test connectivity[1,2] == 0.0
    @test connectivity[2,1] == 1.0

    @test pvalues[1,2] ≈ 0.2666 atol = mattol
    @test pvalues[2,1] ≈ 0.0 atol = mattol
end

@testset "Instantaneous Granger Causality Test" begin
    (_, pf, A, _, _, _) = mvar(u, maxorder=2, method="NS", criterion=nothing)
    connectivity, pvalues = instantaneous_granger_causality_test(u, A, pf)
      
    @test connectivity[1,2] == 0.0
    @test connectivity[2,1] == 0.0

    @test pvalues[1,2] ≈ 0.1160 atol = mattol
    @test pvalues[2,1] ≈ 0.1160 atol = mattol
end