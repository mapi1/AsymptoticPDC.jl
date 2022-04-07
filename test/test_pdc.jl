(model, vic, Vicv) = mvar(u, maxorder=2, method="NS", criterion=nothing)

@testset "Original PDC" begin
    org_pdc = pdc(model, u; nFreqs=2, metric="euc")

    org_pdc_m = [
        [0.999758994313119   0.999211301898553
        0.000241005686881   0.000788698101447];;;
        [0.999986332442013   0.999015012487769
        0.000013667557987   0.000984987512231]
    ]
    @test org_pdc.coherence ≈ org_pdc_m
end

@testset "Asymptotic Original PDC" begin
    org_pdc = pdc(model, u; nFreqs=2, α = 0.01, metric="euc")

    org_pvalue_m = [
        [0.003869198351588   0.104451981460104
        0.000000009250769   0.000000005155287];;;
        [0.000000391761377   0.290133015071853
        0.019921952043625   0.000000003985120]
    ]
    @test org_pdc.pvalues ≈ org_pvalue_m

    org_threshold_m = [
        [0.794958785072552   2.514839682571546
        0.000048466781525   0.000153323905283];;;
        [0.271372595736710   3.800455885612148
        0.000016544953722   0.000231704924285]
    ]
    @test org_pdc.thresholds ≈ org_threshold_m

    org_lower_conf_m = [
        [0.999332765565927   0.996442718554940
        -0.000185223060311  -0.001979885242167];;;
        [0.999964681002367   0.995827616758760
        -0.000007983881660  -0.002202408216779]
    ]
    @test org_pdc.lower_conf ≈ org_lower_conf_m

    org_upper_conf_m = [
        [1.000185223060311   1.001979885242167
        0.000667234434073   0.003557281445060];;;
        [1.000007983881660   1.002202408216779
        0.000035318997633   0.004172383241240]
    ]
    @test org_pdc.upper_conf ≈ org_upper_conf_m
end

@testset "Generalized PDC" begin
    gen_pdc = pdc(model, u; nFreqs=2, metric="diag")

    gen_pdc_m = [
        [0.201858666396383   0.071702345640446
        0.798141333603617   0.928297654359554];;;
        [0.816873188142837   0.058234914675719
        0.183126811857164   0.941765085324281]
    ]
    @test gen_pdc.coherence ≈ gen_pdc_m 
end

@testset "Asymptotic Generalied PDC" begin
    gen_pdc = pdc(model, u; nFreqs=2, α = 0.01, metric="diag")

    gen_pvalue_m = [
        [0.003869198351588   0.104451981460104
        0.000000009250769   0.000000005155287];;;
        [0.000000391761377   0.290133015071853
        0.019921952043625   0.000000003985120]
    ]
    @test gen_pdc.pvalues ≈ gen_pvalue_m

    gen_threshold_m = [
        [0.160508003536476   0.180462234371686
        0.160508003536351   0.180462234371686];;;
        [0.221680027278671   0.221537435835247
        0.221680027278671   0.221537435835247]
    ]
    @test gen_pdc.thresholds ≈ gen_threshold_m

    gen_lower_conf_m = [
        [-0.112007874258400  -0.168358899773962
        0.484274792948835   0.688236408945146];;;
        [0.550300119970548  -0.124963732233905
        -0.083446256315125   0.758566438414657]
    ]
    @test gen_pdc.lower_conf ≈ gen_lower_conf_m

    gen_upper_conf_m = [
        [0.515725207051165   0.311763591054854
        1.112007874258400   1.168358899773962];;;
        [1.083446256315125   0.241433561585343
        0.449699880029452   1.124963732233905]
    ]
    @test gen_pdc.upper_conf ≈ gen_upper_conf_m
end

@testset "Information PDC" begin
    info_pdc = pdc(model, u; nFreqs=2, metric="info")

    info_pdc_m = [
        [0.238666456001199   0.058497174242366
        0.943677905337930   0.757336306794671];;;
        [0.661488474892879   0.057956212201207
        0.148292632498802   0.937257956548512]
    ]
    @test info_pdc.coherence ≈ info_pdc_m 
end

@testset "Asymptotic Information PDC" begin
    info_pdc = pdc(model, u; nFreqs=2, α = 0.01, metric="info")

    info_pvalue_m = [
        [0.003869198351588   0.104451981460104
        0.000000009250769   0.000000005155287];;;
        [0.000000391761377   0.290133015071853
        0.019921952043625   0.000000003985120]
    ]
    @test info_pdc.pvalues ≈ info_pvalue_m

    info_threshold_m = [
        [0.189775732931153   0.147227132863174
        0.189775732931005   0.147227132863174];;;
        [0.179512297976340   0.220477195051720
        0.179512297976340   0.220477195051720]
    ]
    @test info_pdc.thresholds ≈ info_threshold_m

    info_lower_conf_m = [
        [-0.172477143579408  -0.127903474081998
        0.667412185399146   0.301198609152465];;;
        [0.310307056247436  -0.129050103432401
        -0.078816775945015   0.702085022276534]
    ]
    @test info_pdc.lower_conf ≈ info_lower_conf_m

    info_upper_conf_m = [
        [0.649810055581806   0.244897822566730
        1.219943625276715   1.213474004436876];;;
        [1.012669893538322   0.244962527834814
        0.375402040942618   1.172430890820489]
    ]
    @test info_pdc.upper_conf ≈ info_upper_conf_m
end