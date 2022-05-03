var documenterSearchIndex = {"docs":
[{"location":"Utilities/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"Utilities/","page":"Utilities","title":"Utilities","text":"detrend","category":"page"},{"location":"Utilities/#AsymptoticPDC.detrend","page":"Utilities","title":"AsymptoticPDC.detrend","text":"detrend(signal::AbstractVector, order::Int=1; verbose::Bool=true)\n\nDetrend a signal by removing polynomial trend of order p using build in least squares .\n\nArgs:\n\nsignal Data Vector containing te signal\norder=1: order of polynomial\n\nReturn:\n\n(detrended_signal, coefs): The detrended signal and the estimated coefficients from order 0 to  p\n\nExamples\n\njulia> signal = sin.([1:100;]) + 0.03 .* [1:100;] .+ 1\njulia> _, coefficients = detrend(signal)\njulia> coefficients\n2-element Vector{Float64}:\n 1.0583465387418984\n 0.028819440616267268\n\n\n\n\n\ndetrend(signals::AbstractArray, order::Int = 1; verbose = false, dims::Int = 1)\n\nDeterend an array signals along dimension dims.\n\n\n\n\n\n","category":"function"},{"location":"MCAR/#Multichannel-Autoregression-Model","page":"MCAR Model","title":"Multichannel Autoregression Model","text":"","category":"section"},{"location":"MCAR/","page":"MCAR Model","title":"MCAR Model","text":"Definition: A multichannel or multivariate autoregressive model (MCAR) of order p is defined as","category":"page"},{"location":"MCAR/","page":"MCAR Model","title":"MCAR Model","text":"X(t) = sum_i=1^p A(i)X(t-i) + E(t)","category":"page"},{"location":"MCAR/","page":"MCAR Model","title":"MCAR Model","text":"where X(t) = X_1(t) X_2(t) X_k(t) is a k-channel set of signals and E(t) a vector of k withe noises at time t. A is a k  times k times p matrix holding the model parameters.","category":"page"},{"location":"MCAR/","page":"MCAR Model","title":"MCAR Model","text":"mcar","category":"page"},{"location":"MCAR/#AsymptoticPDC.mcar","page":"MCAR Model","title":"AsymptoticPDC.mcar","text":"mcar(u; maxorder::Union{Nothing,Int}=nothing, criterion::Union{Nothing,String}=\"AIC\", method::String=\"NS\", verbose::Bool=true)\n\nCompute a multichannel AR or vector AR model of the input matrix u containing the signals/channels xi, u = [x1 x2 ... xn] \n\nArgs\n\nu: input Matrix containing signals 1 to n u = [x1 x2 ... xn]\n\nKeywords\n\nmaxorder::Union{Nothing, Int} = nothing: The maximal order of the AR model, defaults to nothing where the order is chosen based on a simple heuristic (maxorder = 3√samples/nChannels; Nuttall 1976)\ncriterion = \"AIC\": The information criterion used to choose the model order. Use one of the following:\n\"AIC\": Akaike's Informaion Criterion\n\"HQ\": Hannan Quinn\n\"BIC\": Bayesian Information Criterion, Schwarz 1978\n\"FPE\": Final prediction error, Akaike, 1970\nnothing: maxorder becomes the fixed order\nmethod = \"LS\": Method used for etsimation. Use one of:\n\"LS\" least squares based on \\ \n\"NS\" Nuttall-Strand Method (multi-channel generalization of the single-channel Burg lattice algorithm)\n\"VM\" Vieira-Morf Method (multi-channel generalization of the single-channel geometric lattice algorithm)\n\nReturn\n\nReturns a tuple (model, besticvalue, ic_values)\n\nResult model is of type MCAR_Model, with following fields:\n\norder: is the (chosen) model order\nnChannels: number of channels\nsamples: number of samples per channel\nA: contains the AR coefficients [n x n x order]\npf: is the covariance matrix [order x order]\nef: the residuals\n\n\n\n\n\n","category":"function"},{"location":"MCAR/","page":"MCAR Model","title":"MCAR Model","text":"MCAR_Model","category":"page"},{"location":"MCAR/#AsymptoticPDC.MCAR_Model","page":"MCAR Model","title":"AsymptoticPDC.MCAR_Model","text":"MCAR_Model with following fields:\n\norder: is the (chosen) model order\nnChannels: number of channels\nsamples: number of samples per channel\nA: contains the AR coefficients [n x n x order]\npf: is the covariance matrix [order x order]\nef: the residuals    \n\n\n\n\n\n","category":"type"},{"location":"PDC/#Partial-Directed-Coherence","page":"Partial Directed Coherence","title":"Partial Directed Coherence","text":"","category":"section"},{"location":"PDC/","page":"Partial Directed Coherence","title":"Partial Directed Coherence","text":"pdc","category":"page"},{"location":"PDC/#AsymptoticPDC.pdc","page":"Partial Directed Coherence","title":"AsymptoticPDC.pdc","text":"pdc([model], u; nFreqs::Int = 128, α = 0.0, fs = 1, metric::String=\"euc\", maxorder::Int=30, criterion::Union{Nothing,String}=\"AIC\", method::String=\"NS\", verbose::Bool=true)\n\nComputes the partial directed coherence pdc based on a multivariate AR model of the input matrix u containing the signals/channels xi, u = [x1 x2 ... xn]. If you already fitted a model just put it as first parameter to avoid recalculation. \n\nArgs\n\nu: input Matrix containing signals 1 to n u = [x1 x2 ... xn]\n\nKeywords\n\nnFreqs = 128: number of frequencies for which the pdc shall be calculated \nmetric = \"euc\": If \"euc\" the basic pdc is returned, else if \"diag\" a generalized (normalized for different covariances) pdc is returned\nfs = 1: sampling frequency in Hz\nα = 0.0: significance level for asymptotic statistics. If 0 no statistics are computed which is faster\n\nThe following Keywords are inherited from mcar used for model estimation:\n\n* maxorder::Union{Nothing, Int} = `nothing`: The maximal order of the AR model, defaults to `nothing` where the order is chosen based on a simple heuristic (maxorder = 3√samples/nChannels; Nuttall 1976)\n* criterion = `\"AIC\"`: The information criterion used to choose the model order. Use one of the following:\n    - `\"AIC\"`: Akaike's Informaion Criterion\n    - `\"HQ\"`: Hannan Quinn\n    - `\"BIC\"`: Bayesian Information Criterion, Schwarz 1978\n    - `\"FPE\"`: Final prediction error, Akaike, 1970\n    - nothing: maxorder becomes the fixed order\n* method = `\"LS\"`: Method used for etsimation. Use one of:\n    - `\"LS\"` least squares based on \\ \n    - `\"NS\"` Nuttall-Strand Method (multi-channel generalization of the single-channel Burg lattice algorithm)\n    - `\"VM\"` Vieira-Morf Method (multi-channel generalization of the single-channel geometric lattice algorithm)\n\nReturn\n\nReturns an object that is a subtype of AbstractPartialDirectedCoherence depending on if asymptotic statistics have been calculated.\n\n\n\n\n\n","category":"function"},{"location":"Examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Some examples for demonstration purposes.","category":"page"},{"location":"Examples/#Sunspots-and-Melanoma","page":"Examples","title":"Sunspots and Melanoma","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"The etiology of melanoma is complex and may include the influences of trauma, heredity and hormonal activity. In particular, exposure to solar radiation may be involved in the pathogenesis of melanoma. Melanoma is more common in fair-skinned individuals and most frequent in skin sites exposed to the sun. In white populations melanoma is more common in areas closer to the equator where the intensity of solar radiation is higher. Data from various parts of the world suggest that the incidence of melanoma is increasing. The data below, giving age-adjusted melanoma incidence, are from the Connecticut Tumor Registry from 1936-1972. Connecticut has the longest record of state population-based cancer statistics in the United States of America. The data also includes the sunspot relative number. Houghton, Munster and Viola (1978) have shown that the age-adjusted incidence rate for malignant melanoma in the state of Connecticut has risen since 1935 and that superimposed on the rise are 3-5 year periods in which the rise in the rate of incidence is excessive. These periods have a cycle of 8-11 years and follow times of maximum sunspot activity. The relationship between solar cycles and melanoma supports the hypothesis that melanoma is related to sun exposure and provides evidence that solar radiation may trigger the development of clinically apparent melanoma.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"The data can be obtained through the following helper function","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"get_sunspot_melanoma_data","category":"page"},{"location":"Examples/#AsymptoticPDC.get_sunspot_melanoma_data","page":"Examples","title":"AsymptoticPDC.get_sunspot_melanoma_data","text":"get_sunspot_melanoma_data()\n\nThis data is from Andrews and Herzberg. D. F. Andrews, A. M. Herzberg. (1985) Data: A Collection of Problems from Many Fields for the Student and Research Worker. Springer, New York.\n\nReturns u:\n\nu[:, 1]: Year 1936-1972 \nu[:, 2]: Annual male melanoma incidence (age-adjusted per 10e5) in Connecticut\nu[:, 3]: Annual total melanoma incidence (age-adjusted per 10e5) in Connecticut\nu[:, 4]: Annual sunspot relative number \n\n\n\n\n\n","category":"function"},{"location":"Examples/#Setup","page":"Examples","title":"Setup","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"A first look at the detrended data","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using Plots, AsymptoticPDC\ndefault(dpi = 200)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"u = get_sunspot_melanoma_data()\nt = u[:, 1]\nmelanoma, _ = detrend(u[:, 3])\nsunspots, _ = detrend(u[:, 4])\np1 = plot(t, melanoma, title = \"Detrended annual melanoma incidence Connecticut\")\np2 = plot(t, sunspots, title = \"Detrended annual sunspot number\")\nplot(p1, p2, layout = (2, 1), legend = :none)","category":"page"},{"location":"Examples/#Estimation","page":"Examples","title":"Estimation","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Estimate the model parameters using mcar","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"y = [melanoma sunspots]\nmodel, _, _ = mcar(y)","category":"page"},{"location":"Examples/#Connectivity","page":"Examples","title":"Connectivity","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Perform a Granger causality test to get a connectivity matrix","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"connectivity, pvalues = granger_causality_test(model, y);\nprintln(\"Connectivity matrix:\") # hide\nshow(stdout, \"text/plain\", connectivity) # hide\nprintln(\"Granger causality test p-values:\") # hide\nshow(stdout, \"text/plain\", pvalues) # hide","category":"page"},{"location":"Examples/#PDC-Analysis","page":"Examples","title":"PDC Analysis","text":"","category":"section"},{"location":"Examples/#Original-PDC-analysis","page":"Examples","title":"Original PDC analysis","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"The time series are not standardized and the sunspots series has a larger variance which explains the counterintuitive results below","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"original_pdc = pdc(model, y; metric = \"euc\")\npdcplot(original_pdc, cnames = [\"Melanoma\", \"Sunspots\"])","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"But when calculating the asymptotic statistics for the original PDC analysis, only the coherence from sunspots -> melanoma is considered significant (grey shaded area). ","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"original_pdc_asymp = pdc(model, y; metric = \"euc\", α = 0.01)\npdcplot(original_pdc_asymp, cnames = [\"Melanoma\", \"Sunspots\"])","category":"page"},{"location":"Examples/#Generalized-PDC-analysis","page":"Examples","title":"Generalized PDC analysis","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"The generalized PDC analysis gives a better view inside the \"true\" interaction. The ~11 year sun activity cycle is clearly visible","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"generalized_pdc = pdc(model, y; metric = \"diag\", α = 0.01)\npdcplot(generalized_pdc, cnames = [\"Melanoma\", \"Sunspots\"])","category":"page"},{"location":"Examples/#Information-PDC-analysis","page":"Examples","title":"Information PDC analysis","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"The information PDC analysis gives similar results","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"info_pdc = pdc(model, y; metric = \"info\", α = 0.01)\npdcplot(info_pdc, cnames = [\"Melanoma\", \"Sunspots\"])","category":"page"},{"location":"Stat_Tests/#Statistical-Tests","page":"Statistical Tests","title":"Statistical Tests","text":"","category":"section"},{"location":"Stat_Tests/","page":"Statistical Tests","title":"Statistical Tests","text":"granger_causality_test","category":"page"},{"location":"Stat_Tests/#AsymptoticPDC.granger_causality_test","page":"Statistical Tests","title":"AsymptoticPDC.granger_causality_test","text":"grangercausalitytest(model, u; α = 0.01, verbose::Bool = true)\n\n\n\n\n\n","category":"function"},{"location":"Stat_Tests/","page":"Statistical Tests","title":"Statistical Tests","text":"instantaneous_granger_causality_test","category":"page"},{"location":"Stat_Tests/#AsymptoticPDC.instantaneous_granger_causality_test","page":"Statistical Tests","title":"AsymptoticPDC.instantaneous_granger_causality_test","text":"instantaneousgrangercausality_test(model, u; α = 0.01, verbose::Bool = true)\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = AsymptoticPDC","category":"page"},{"location":"#AsymptoticPDC.jl","page":"Home","title":"AsymptoticPDC.jl","text":"","category":"section"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"MCAR.md\",\n    \"PDC.md\",\n    \"Stat_Tests.md\",\n    \"Utilities.md\",\n    \"Examples.md\"\n]\nDepth = 2","category":"page"}]
}