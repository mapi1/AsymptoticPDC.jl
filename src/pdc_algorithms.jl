"""
    pdc([model], u; nFreqs::Int = 128, α = 0.0, fs = 1, metric::String="euc", maxorder::Int=30, criterion::Union{Nothing,String}="AIC", method::String="NS", verbose::Bool=true)

Computes the partial directed coherence pdc based on a multivariate AR model of the input matrix u containing the signals/channels xi, u = [x1 x2 ... xn].
If you already fitted a model just put it as first parameter to avoid recalculation. 

# Args

* `u`: input Matrix containing signals 1 to n u = [x1 x2 ... xn]

# Keywords

* `nFreqs` = 128: number of frequencies for which the pdc shall be calculated 
* `metric` = "euc": If "euc" the basic pdc is returned, else if "diag" a generalized (normalized for different covariances) pdc is returned
* `fs` = 1: sampling frequency in Hz
* `α` = 0.0: significance level for asymptotic statistics. If 0 no statistics are computed which is faster

The following Keywords are inherited from `mcar` used for model estimation:
    
    * maxorder::Union{Nothing, Int} = `nothing`: The maximal order of the AR model, defaults to `nothing` where the order is chosen based on a simple heuristic (maxorder = 3√samples/nChannels; Nuttall 1976)
    * criterion = `"AIC"`: The information criterion used to choose the model order. Use one of the following:
        - `"AIC"`: Akaike's Informaion Criterion
        - `"HQ"`: Hannan Quinn
        - `"BIC"`: Bayesian Information Criterion, Schwarz 1978
        - `"FPE"`: Final prediction error, Akaike, 1970
        - nothing: maxorder becomes the fixed order
    * method = `"LS"`: Method used for etsimation. Use one of:
        - `"LS"` least squares based on \\ 
        - `"NS"` Nuttall-Strand Method (multi-channel generalization of the single-channel Burg lattice algorithm)
        - `"VM"` Vieira-Morf Method (multi-channel generalization of the single-channel geometric lattice algorithm)

# Return 

Returns an object that is a subtype of `AbstractPartialDirectedCoherence` depending on if asymptotic statistics have been calculated.
"""
function pdc(u; nFreqs::Int = 128, α = 0.0, fs = 1, metric::String="euc", maxorder::Int=30, criterion::Union{Nothing,String}="AIC", method::String="NS", verbose::Bool=true)
    model, _, _ = mcar(u, maxorder=maxorder, criterion=criterion, method=method, verbose=verbose)
    return pdc(model, u; nFreqs=nFreqs, α=α, fs=fs, metric=metric, verbose=verbose)
end

function pdc(model::MCAR_Model, u; nFreqs::Int = 128, α = 0.0, fs = 1, metric::String="euc", verbose::Bool=true)
    if lowercase(metric) == "euc"
        return original_pdc(model, u; nFreqs=nFreqs, α = α, fs = fs)
    elseif lowercase(metric) == "diag"
        return generalized_pdc(model, u; nFreqs=nFreqs, α = α, fs = fs)
    elseif lowercase(metric) == "info"
        return information_pdc(model, u; nFreqs=nFreqs, α = α, fs = fs)
    else
        throw(DomainError(metric, "Invalid metric"))
    end
end

abstract type AbstractPartialDirectedCoherence end

struct PartialDirectedCoherence{T,F,A<:AbstractArray{T,3},B<:AbstractArray{Complex{T},3}} <: AbstractPartialDirectedCoherence
    freq::F
    coherence::A
    spectra::B
end

struct AsymptoticPartialDirectedCoherence{T,F,A<:AbstractArray{T,3},B<:AbstractArray{Complex{T},3}} <: AbstractPartialDirectedCoherence
    freq::F
    coherence::A
    spectra::B
    lower_conf::A
    upper_conf::A
    pvalues::A
    thresholds::A
end

###
function original_pdc(model, u; nFreqs=128, α=0.0, fs=1)
    Af = A2f(model.A, nFreqs)
    nChannels = model.nChannels
    pdc_res = zeros(nChannels, nChannels, nFreqs)

    # get spectra as well, as AR is already fitted
    spec = _spectra(Af, model.pf)

    # asymptotic statistics
    calculate_asymptotic = !iszero(α)
    if calculate_asymptotic
        intermediary = init_asymptotic_properties(model, nFreqs, u, α)
    end

    for f in 1:nFreqs
        a = vec(Af[f, :, :])
        a = [real(a); imag(a)]

        # asymptotic statistics
        calculate_asymptotic && update_on_f!(intermediary, f, model, nFreqs)
        for j in 1:nChannels
            Ije = fIj(j, nChannels)
            for i in 1:nChannels
                Iije = fIij(i, j, nChannels)
                num = a' * Iije * a
                den = a' * Ije * a
                pdc_res[i, j, f] = num / den

                # asymptotic statistics
                if calculate_asymptotic
                    dpdc_dev = zeros(1, div(nChannels * (nChannels + 1), 2))
                    calc_asymptotic_properties!(intermediary, model, α, pdc_res, Iije, Ije, i, j, f, num, den, a, dpdc_dev)
                end
            end
        end
    end
    freq = range(0, 0.5fs, length = nFreqs)
    if calculate_asymptotic
        return AsymptoticPartialDirectedCoherence(freq, pdc_res, spec, intermediary.lower_conf, intermediary.upper_conf, intermediary.pvalues, intermediary.threshold) 
    else
        return PartialDirectedCoherence(freq, pdc_res, spec) 
    end
end

# struct that holds intermediate results while calculating the asymptotic statistics
mutable struct IntermediateAsymptotic
    pvalues
    threshold
    lower_conf
    upper_conf
    varass1
    varass2
    patdfr
    patdenr
    Ca
    omega2
    L
    gamma
    omega
    omega_evar
    icdf_norm_alpha
end

function init_asymptotic_properties(model, nFreqs, u, α)
    # depend on f,i,j
    pvalues = Array{Float64}(undef, model.nChannels, model.nChannels, nFreqs)
    threshold = similar(pvalues)
    lower_conf = similar(pvalues)
    upper_conf = similar(pvalues)
    varass1 = similar(pvalues)
    varass2 = similar(pvalues)
    patdfr = similar(pvalues)
    patdenr = similar(pvalues)

    # independent
    gamma = bigautocorr(u, model.order)
    omega = kron(inv(gamma), model.pf)
    inv_D = pinv(dmatrix(model.nChannels))
    omega_evar = 2 * inv_D * kron(model.pf, model.pf) * inv_D'
    icdf_norm_alpha = invlogcdf(Normal(0, 1), log(1 - 0.5α))

    # depend on f 
    f_ind = 0
    Ca = fCa(f_ind, model.order, model.nChannels)
    omega2 = Ca * omega * Ca'
    L = fChol(omega2)
    return IntermediateAsymptotic(pvalues, threshold, lower_conf, upper_conf, varass1, varass2, patdfr, patdenr, Ca, omega2, L, gamma, omega, omega_evar, icdf_norm_alpha)
end

function update_on_f!(intermediary, f, model, nFreqs)
    f_ind = (f - 1) / (2 * nFreqs) #%Corrected 7/25/2011, f starting at 0 rad/s.
    intermediary.Ca = fCa(f_ind, model.order, model.nChannels)
    intermediary.omega2 = intermediary.Ca * intermediary.omega * intermediary.Ca'
    intermediary.L = fChol(intermediary.omega2)
end

function calc_asymptotic_properties!(intermediary, model, α, pdc_res, Iije, Ije, i, j, f, num, den, a, dpdc_dev)
    G1a = 2a' * Iije / den - 2num * a' * Ije / (den^2)
    G1 = -G1a * intermediary.Ca
    varalpha = G1 * intermediary.omega * G1'
    varevar = dpdc_dev * intermediary.omega_evar * dpdc_dev'
    intermediary.varass1[i, j, f] = (varalpha + varevar[1]) / model.samples

    intermediary.lower_conf[i, j, f] = pdc_res[i, j, f] - sqrt(intermediary.varass1[i, j, f]) * intermediary.icdf_norm_alpha
    intermediary.upper_conf[i, j, f] = pdc_res[i, j, f] + sqrt(intermediary.varass1[i, j, f]) * intermediary.icdf_norm_alpha

    G2a = Iije / den
    d = fEig(real(intermediary.L), real(G2a))

    patdf = (sum(d) .^ 2) ./ sum(d .^ 2)
    patden = sum(d) ./ sum(d .^ 2)

    test_dist = Chisq(patdf)
    intermediary.threshold[i, j, f] = invlogcdf(test_dist, log(1 - α)) / (patden * model.samples)
    intermediary.pvalues[i, j, f] = 1 - cdf(test_dist, pdc_res[i, j, f] * patden * model.samples)

    intermediary.varass2[i, j, f] = patdf / (patden * model.samples) .^ 2
    intermediary.patdfr[i, j, f] = patdf
    intermediary.patdenr[i, j, f] = patden
end
###
function generalized_pdc(model, u; nFreqs=128, α=0.0, fs=1)
    nChannels = model.nChannels
    Af = A2f(model.A, nFreqs)
    pdc_res = zeros(nChannels, nChannels, nFreqs)

    evar_d = diagm(diag(model.pf))
    evar_d_big = kron(I(2nChannels), evar_d)
    pinv_evar_d_big = pinv(evar_d_big)

    # get spectra as well, as AR is already fitted
    spec = _spectra(Af, model.pf)

    # asymptotic statistics
    calculate_asymptotic = !iszero(α)
    if calculate_asymptotic
        de_deh = dmatrix(nChannels)
        debig_de = fdebig_de(nChannels)
        dedinv_dev = spdiagm(sparse(vec(-pinv_evar_d_big*pinv_evar_d_big)))
        dedinv_deh = dedinv_dev*debig_de*de_deh
        intermediary = init_asymptotic_properties(model, nFreqs, u, α)
    end

    for f in 1:nFreqs
        a = vec(Af[f, :, :])
        a = [real(a); imag(a)]

        # asymptotic statistics
        calculate_asymptotic && update_on_f!(intermediary, f, model, nFreqs)
        for j in 1:nChannels
            Ij = fIj(j, nChannels)
            Ije = Ij * pinv_evar_d_big

            for i in 1:nChannels
                Iij = fIij(i, j, nChannels)
                Iije = Iij * pinv_evar_d_big
                num = a' * Iije * a
                den = a' * Ije * a
                pdc_res[i, j, f] = num / den

                # asymptotic statistics
                if calculate_asymptotic
                    dnum_dev = kron((Iij*a)', a') * dedinv_deh
                    dden_dev = kron((Ij*a)', a') * dedinv_deh
                    dpdc_dev = (den*dnum_dev - num*dden_dev)/(den^2)
                    calc_asymptotic_properties!(intermediary, model, α, pdc_res, Iije, Ije, i, j, f, num, den, a, dpdc_dev)
                end
            end
        end
    end
    freq = range(0, 0.5fs, length = nFreqs)
    if calculate_asymptotic
        return AsymptoticPartialDirectedCoherence(freq, pdc_res, spec, intermediary.lower_conf, intermediary.upper_conf, intermediary.pvalues, intermediary.threshold) 
    else
        return PartialDirectedCoherence(freq, pdc_res, spec) 
    end 
end

function information_pdc(model, u; nFreqs=128, α=0.0, fs=1)
    nChannels = size(model.pf, 1)

    Af = A2f(model.A, nFreqs)
    pdc_res = zeros(nChannels, nChannels, nFreqs)

    evar_d = diagm(diag(model.pf))
    evar_d_big = kron(I(2nChannels), evar_d)
    pinv_evar_d_big = pinv(evar_d_big)

    evar_big = kron(I(2nChannels), model.pf)
    pinv_evar_big = pinv(evar_big)

    # get spectra as well, as AR is already fitted
    spec = _spectra(Af, model.pf)

    # asymptotic statistics
    calculate_asymptotic = !iszero(α)
    if calculate_asymptotic
        de_deh = dmatrix(nChannels)
        debig_de = fdebig_de(nChannels)
        dedinv_devd = spdiagm(sparse(vec(-pinv_evar_d_big*pinv_evar_d_big)))
        dedinv_dehd = dedinv_devd*debig_de*de_deh
        inv_e = sparse(pinv_evar_big) 
        dedinv_dev = -kron(inv_e', inv_e)
        dedinv_deh = sparse(dedinv_dev*debig_de*de_deh)
        intermediary = init_asymptotic_properties(model, nFreqs, u, α)
    end
    
    for f in 1:nFreqs
        a = vec(Af[f, :, :])
        a = [real(a); imag(a)]

        # asymptotic statistics
        calculate_asymptotic && update_on_f!(intermediary, f, model, nFreqs)
        for j in 1:nChannels
            Ij = fIj(j, nChannels)
            Ije = Ij * pinv_evar_big * Ij

            for i in 1:nChannels
                Iij = fIij(i, j, nChannels)
                Iije = Iij * pinv_evar_d_big
                num = a' * Iije * a
                den = a' * Ije * a
                pdc_res[i, j, f] = num / den

                 # asymptotic statistics
                 if calculate_asymptotic
                    dnum_dev = kron((Iij*a)', a') * dedinv_dehd
                    dden_dev = kron((Ij*a)', a') * dedinv_deh
                    dpdc_dev = (den*dnum_dev - num*dden_dev)/(den^2)
                    calc_asymptotic_properties!(intermediary, model, α, pdc_res, Iije, Ije, i, j, f, num, den, a, dpdc_dev)
                end
            end
        end
    end
    freq = range(0, 0.5fs, length = nFreqs)
    if calculate_asymptotic
        return AsymptoticPartialDirectedCoherence(freq, pdc_res, spec, intermediary.lower_conf, intermediary.upper_conf, intermediary.pvalues, intermediary.threshold) 
    else
        return PartialDirectedCoherence(freq, pdc_res, spec) 
    end 
end