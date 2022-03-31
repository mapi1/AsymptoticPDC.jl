(res, vic, Vicv) = mvar((u), maxorder = 2, method = "NS", criterion = nothing)

@testset "Spectra" begin
    
    spec = spectra(res, nFreqs=2)
    spec_m = [[3612.8692960852+0im 10.5775742584943+0im
        10.5775742584943+0im 0.0583557154306129+0im];;;
    [439.214351139546+0im -1.42438768963195+1.65703566679711im
        -1.42438768963195-1.65703566679711im 0.0641046201218521-3.46944695195361e-18im]]
    
    @test spec ≈ spec_m

end

@testset "Coherence" begin
    coh = coherence(res, nFreqs=2)

    coh_m = [
        [1+0im 0.728481049938064+0im
            0.728481049938064+0im 1+0im];;;
        [1+0im -0.268438797758722+0.312283422186324im
            -0.268438797758722-0.312283422186324im 1+0im]
    ]
    
    @test coh ≈ coh_m
end

