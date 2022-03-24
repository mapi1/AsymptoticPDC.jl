@testset "Sunspot Melanoma" begin
    y = get_sunspot_melanoma_data()
    @test size(y) == (37, 4)
end