# method structure outline

struct PartialDirectedCoherence{T,F,A<:AbstractArray{T,3}}
    coherence::A
    freq::F
end

struct AsymptoticPartialDirectedCoherence{T,F,A<:AbstractArray{T,3}}
    freq::F
    coherence::A
    lower_conf::A
    upper_conf::A
    pvalues::A
    thresholds::A
end

struct CrossPowerSpectra{T,F,A<:AbstractArray{T,3}}
    power::A
    freq::F
end

function original_pdc(model; nFreqs=128, α=0.0)
    nChannels, _, order = size(model.A)

    Af = A2f(model.A, nFreqs)
    pdc_res = zeros(nChannels, nChannels, nFreqs)

    # get spectra as well, as AR is already fitted
    spec = _spectra(Af, model.pf)
    coh = _coherence(spec)

    # asymptotic statistics
    calculate_asymptotic = !iszero(α)
    calculate_asymptotic && intermediary = init_asymptotic_properties(model, nFreqs, u, α)

    for f in 1:nFreqs
        a = vec(Af[f, :, :])
        a = [real(a); imag(a)]

        # pdc asymp
        f_ind = (f - 1) / (2 * nFreqs) #%Corrected 7/25/2011, f starting at 0 rad/s.
        Ca = fCa(f_ind, order, nChannels)
        omega2 = Ca * omega * Ca'
        L = fChol(omega2)
        for j in 1:nChannels
            Ije = fIj(j, nChannels)
            for i in 1:nChannels
                Iije = fIij(i, j, nChannels)
                num = a' * Iije * a
                den = a' * Ije * a
                pdc_res[i, j, f] = num / den

                # asymptotic statistics
                if calculate_asymptotic
                    ## euc
                    dpdc_dev = zeros(div(nChannels * (nChannels + 1), 2))
                    ## diag
                    
                    ## info
                    
                   calc_asymptotic_properties!(intermediary, α, pdc_res, Iije, Ije, num, den, a, Ca, L, dpdc_dev)
                end
            end
        end
    end
    return (pdc_res)
end

function init_asymptotic_properties(model, nFreqs, u, α)
    pvalues = Array{Float64}(undef, nChannels, nChannels, nFreqs)
    threshold = Array{Float64}(undef, nChannels, nChannels, nFreqs)
    lower_conf = Array{Float64}(undef, nChannels, nChannels, nFreqs)
    upper_conf = Array{Float64}(undef, nChannels, nChannels, nFreqs)
    varass1 = Array{Float64}(undef, nChannels, nChannels, nFreqs)
    varass2 = Array{Float64}(undef, nChannels, nChannels, nFreqs)
    patdfr = Array{Float64}(undef, nChannels, nChannels, nFreqs)
    patdenr = Array{Float64}(undef, nChannels, nChannels, nFreqs)

    gamma = bigautocorr(u, model.order)
    omega = kron(inv(gamma), model.pf)
    inv_D = pinv(dmatrix(nChannels))
    omega_evar = 2 * inv_D * kron(model.pf, model.pf) * inv_D'
    icdf_norm_alpha = invlogcdf(Normal(0, 1), log(1 - 0.5α))
    return IntermediateAsymptotic(pvalues, threshold, lower_conf, upper_conf, varass1, varass2, patdfr, patdenr, gamma, omega, omega_evar, icdf_norm_alpha)
end

struct IntermediateAsymptotic
    pvalues
    threshold
    lower_conf
    upper_conf
    varass1
    varass2
    patdfr
    patdenr
    gamma
    omega
    omega_evar
    icdf_norm_alpha
end

function calc_asymptotic_properties!(intermediary, α, pdc_res, Iije, Ije, num, den, a, Ca, L, dpdc_dev)
    G1a = 2a' * Iije / den - 2num * a' * Ije / (den^2)
    G1 = -G1a * Ca
    varalpha = G1 * intermediary.omega * G1'
    varevar = dpdc_dev' * omega_evar * dpdc_dev
    intermediary.varass1[i, j, f] = (varalpha + varevar) / samples

    intermediary.lower_conf[i, j, f] = pdc_res[i, j, f] - sqrt(intermediary.varass1[i, j, f]) * intermediary.icdf_norm_alpha
    intermediary.upper_conf[i, j, f] = pdc_res[i, j, f] + sqrt(intermediary.varass1[i, j, f]) * intermediary.icdf_norm_alpha

    G2a = Iije / den
    d = fEig(real(L), real(G2a))

    patdf = (sum(d) .^ 2) ./ sum(d .^ 2)
    patden = sum(d) ./ sum(d .^ 2)

    test_dist = Chisq(patdf)
    intermediary.threshold[i, j, f] = invlogcdf(test_dist, log(1 - α)) / (patden * samples)
    intermediary.pvalues[i, j, f] = 1 - cdf(test_dist, pdc_res[i, j, f] * patden * samples)

    intermediary.varass2[i, j, f] = patdf / (patden * samples) .^ 2
    intermediary.patdfr[i, j, f] = patdf
    intermediary.patdenr[i, j, f] = patden
end
