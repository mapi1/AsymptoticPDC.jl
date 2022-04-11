"""
    spectra(model::MCAR_Model; nFreqs::Int=128, fs=1)

Calculate the cross power spectra from an existing MCAR model.
"""
function spectra(model::MCAR_Model; nFreqs::Int=128, fs=1)
    Af = A2f(model.A, nFreqs)
    spec = _spectra(Af, model.pf)
    freq = range(0, 0.5fs, length=nFreqs)
    return DSP.Periodograms.CrossPowerSpectra(spec, freq)
end

function _spectra(Af, pf)
    nChannels = size(Af, 2)
    nFreqs = size(Af, 1)
    spec = zeros(eltype(Af), nChannels, nChannels, nFreqs)
    for f in 1:nFreqs # could be Symmetric
        H = inv(Af[f, :, :])
        spec[:, :, f] = H * pf * H'
    end
    return spec
end

"""
    coherence(model::MCAR_Model; nFreqs::Int=128, fs=1)
    
Calculate the coherence from an existing MCAR model.
"""
function coherence(model::MCAR_Model; nFreqs::Int=128, fs=1)
    spec = spectra(model; nFreqs=nFreqs, fs=fs)
    coh = _coherence(spec.power)
    return DSP.Periodograms.Coherence(coh, spec.freq)

end

function _coherence(spec)
    nChannels = size(spec, 1)
    nFreqs = size(spec, 3)

    coh = similar(spec)
    for f in 1:nFreqs
        for ii in 1:nChannels, jj in 1:nChannels # could be Symmetric
            coh[jj, ii, f] = spec[jj, ii, f] / sqrt(spec[jj, jj, f] * spec[ii, ii, f])
        end
    end
    return coh
end