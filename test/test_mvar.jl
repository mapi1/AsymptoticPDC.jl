# Prepare some data
y = get_sunspot_melanoma_data()
melanoma, _ = detrend(y[:, 3], 1)
sunspots, _ = detrend(y[:, 4], 1)
u = [sunspots melanoma]

mattol = 0.001 # tolerance from matlab

@testset "MVAR - Nuttal Strand" begin
    (order, pf, A, ef, vic, Vicv) = mvar(u, maxorder=2, method="NS", criterion=nothing)
    
    @test order == 2

    A_NS_M = [
        [0.9910  -23.4414
            0.0008   -0.0754];;;
        [-0.3250  -17.7061
            0.0044   -0.0807]
    ]   
    @test A ≈ A_NS_M atol = mattol

    pf_NS_M = [792.6836   -1.6557
        -1.6557    0.0483]
    @test pf ≈ pf_NS_M atol = mattol
end

@testset "MVAR - LS" begin
    (order, pf, A, ef, vic, Vicv) = mvar(u, maxorder=2, method="LS", criterion=nothing)

    @test order == 2
    A_LS_M = [
        [0.9515 -27.2894
            0.0009 -0.0722];;;
        [-0.2654 -18.8789
            0.0043 -0.0780]
    ]
    @test A ≈ A_LS_M atol = mattol

    pf_LS_M = [952.0619 -2.1079
        -2.1079 0.0562]
    @test pf ≈ pf_LS_M atol = mattol

    
end