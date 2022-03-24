
"""
pdc(u::Matrix{<:Real}; nFreqs::Int = 124, metric::String = "euc", maxorder::Int = 30, criterion::Union{Nothing, String} = "AIC", method::String = "NS")

Computes the partial directed coherence pdc based on a multivariate AR model of the input matrix u containing the signals/channels xi, u = [x1 x2 ... xn] 

# Args

* u::Matrix{<:Real}: input Matrix containing signals 1 to n u = [x1 x2 ... xn]

# Keywords

* nFreqs::Int = 124: number of frequencies for which the pdc shall be calculated 
* metric::String = "euc": If "euc" the basic pdc is returned, else if "diag" a generalized (normalized for different covariances) pdc is returned
* maxorder::Union{Nothing, Int} = nothing: The maximal order of the AR model, defaults to nothing where the order is choosen based on length and the number of channels (Nuttall 1976)
* criterion::String = "AIC": The information criterion used to choose the model order. Use one of the following:
    - "AIC": Akaike's Informaion Criterion
    - "HQ": Hannan Quinn
    - "BIC": Bayesian Information Criterion, Schwarz 1978
    - "FPE": Final prediction error (Akaike, 1970)
    - nothing: maxorder becomes the fixed order
* method::String = "LS": Method used for etsimation. Use one of:
    - "LS" least squares based on \\ 
    - "NS" Nuttall Strand

# Return 

Return the tuple (pdc , spectra, coh), where:

* pdc: partial directed coherence [n x n x nFreqs]
* spectra: (Cross-) spectra of the signals [n x n x nFreqs]
* coh: normal coherence [n x n x nFreqs]
"""
function pdc(u::Matrix{<:Real}; nFreqs::Int = 124, metric::String = "euc", maxorder::Int = 30, criterion::Union{Nothing, String} = "AIC", method::String = "NS")
    len, nChannels = size(u) 
    order,pf,A,ef,vic,Vicv = mvar(u, maxorder = maxorder, criterion = criterion, method = method)
    
    Af = A2f(A, nFreqs)
    pdc_res = zeros(nChannels,nChannels,nFreqs)
    # if metric == "euc"
    #     dpdc_dev = zeros(1, Int((nChannels*(nChannels+1))/2));
    # else
    #     throw(ArgumentError("Invalid metric, use 'euc' for original pdc"))
    # end
    if metric == "diag"
        evar_d = diagm(diag(pf))
        evar_d_big = kron(Matrix(I,2nChannels, 2nChannels), evar_d);
        pinv_evar_d_big = pinv(evar_d_big);   
    end
    # gamma = bigautocorr(u, order)
    # omega = kron(inv(gamma), pf)
    # omega_evar = 2*pinv(Dup(nChannels)) * kron(pf, pf) * pinv(Dup(nChannels))'
    
    # get spectra as well, as AR is already fitted
    spectra = zeros(Complex, nChannels,nChannels,nFreqs)
    coh = copy(spectra)
    for f in 1:nFreqs
        freq = (f-1) / 2nFreqs
        # Ca = fCa(freq, p, nChannels)
        a = vec(Af[f,:,:]) 
        a = [real(a); imag(a)]
        # println(a)
        # omega2 = Ca * omega * Ca'
        # L = fChol(omega2)
        
        # pdc
        for jj in 1:nChannels
            Ij = fIj(jj,nChannels);
            if metric == "euc" # basic pdc
                Ije = Ij;
            elseif metric == "diag" # gerneralized pdc
                Ije = Ij * pinv_evar_d_big
            else
                throw(ArgumentError("Invalid metric, use 'euc' for original pdc or 'diag' for a generalized pdc."))
            end
            
            for ii in 1:nChannels
                Iij = fIij(ii,jj,nChannels);
                if metric == "euc"
                    Iije = Iij;
                elseif metric == "diag"
                    Iije = Iij * pinv_evar_d_big
                else
                end
                num = a'*Iije*a
                den = a'*Ije*a
                pdc_res[ii,jj,f] = num/den
            end
        end
        
        # spectra
        H = inv(reshape(Af[f,:,:],nChannels, nChannels))
        spectra[:,:, f] = H*pf*H'
        
        # coherence (Sij/sqrt(Sii * Sij))
        
        for jj in 1:nChannels, ii in 1:nChannels
            coh[jj, ii, f] = spectra[jj, ii, f] / sqrt(spectra[jj, jj, f] * spectra[ii, ii, f])
        end
    end
    
    return (pdc_res, spectra, coh)
end