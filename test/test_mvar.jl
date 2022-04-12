mattol = 0.001 # tolerance from matlab

@testset "MVAR - Nuttal Strand" begin
    (model, vic, Vicv) = mcar(u, maxorder=2, method="NS", criterion=nothing)
    
    @test model.order == 2

    A_NS_M = [
        [0.9910  -23.4414
            0.0008   -0.0754];;;
        [-0.3250  -17.7061
            0.0044   -0.0807]
    ]   
    @test model.A ≈ A_NS_M atol = mattol

    pf_NS_M = [792.6836   -1.6557
        -1.6557    0.0483]
    @test model.pf ≈ pf_NS_M atol = mattol
end

@testset "MVAR - Vieira Morf" begin
    (model, vic, Vicv) = mcar(u, maxorder=2, method="VM", criterion=nothing)
    
    @test model.order == 2

    A_VM_M = [
        [0.9821  -25.3824
            0.0008   -0.0738];;;
        [-0.3072  -15.8521
            0.0043   -0.0823]
    ]   
    @test model.A ≈ A_VM_M atol = mattol

    pf_VM_M = [791.3738   -1.6861
        -1.6116    0.0483]
    @test model.pf ≈ pf_VM_M atol = mattol
end

@testset "MVAR - LS" begin
    (model, vic, Vicv) = mcar(u, maxorder=2, method="LS", criterion=nothing)

    @test model.order == 2
    A_LS_M = [
        [0.9515 -27.2894
            0.0009 -0.0722];;;
        [-0.2654 -18.8789
            0.0043 -0.0780]
    ]
    @test model.A ≈ A_LS_M atol = mattol

    pf_LS_M = [952.0619 -2.1079
        -2.1079 0.0562]
    @test model.pf ≈ pf_LS_M atol = mattol

    
end