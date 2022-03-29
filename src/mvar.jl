"""
mvar(u::Matrix{<:Real}; maxorder::Union{Nothing, Int} = nothing, criterion::Union{Nothing, String} = "AIC", method::String = "LS"

Compute a multivariate AR or vector AR model of the input matrix u containing the signals/channels xi, u = [x1 x2 ... xn] 

# Args

* u::Matrix{<:Real}: input Matrix containing signals 1 to n u = [x1 x2 ... xn]

# Keywords

* maxorder::Union{Nothing, Int} = `nothing`: The maximal order of the AR model, defaults to `nothing` where the order is chosen based on length and the number of channels (Nuttall 1976)
* criterion::String = `AIC`: The information criterion used to choose the model order. Use one of the following:
    - `AIC`: Akaike's Informaion Criterion
    - `HQ`: Hannan Quinn
    - `BIC`: Bayesian Information Criterion, Schwarz 1978
    - `FPE`: Final prediction error, Akaike, 1970
    - nothing: maxorder becomes the fixed order
* method::String = `LS`: Method used for etsimation. Use one of:
    - `LS` least squares based on \\ 
    - `NS` Nuttall-Strand Method (multi-channel generalization of the single-channel Burg lattice algorithm)
    - `VM` Vieira-Morf Method (multi-channel generalization of the single-channel geometric lattice algorithm)

# Return 

Return the tuple (order, pf, A, ef, vic, Vicv), where:

* order: is the (chosen) model order
* pf: is the covariance matrix [order x order]
* A: contains the AR coefficients [n x n x order]
* ef: the residuals
* vic: the 'optimal' IC value
* Vicv: Vector with all IC values
"""
function mvar(u; maxorder::Union{Nothing,Int}=nothing, criterion::Union{Nothing,String}="AIC", method::String="NS", verbose::Bool=true)
    len, nChannels = size(u)
    verbose && @info "$nChannels channels with $len samples."
    # len > nChannels || throw(ArgumentError("The input u is probably transposed, the number of channels is larger than the channel lengths."))

    if maxorder === nothing
        maxorder = Int(round(3 * sqrt(len) / nChannels)) # Suggested by Nuttall, 1976.
        verbose && @info "The maxorder was choosen as $maxorder."
    else
        maxorder >= 1 || throw(DomainError("The maximal order needs to be >= 1, is: $maxorder."))
    end

    # assign method
    if method == "LS"
        verbose && @info "Using least squares method"
        method_ar = (u, order) -> mc_ar_ls(u, order)
    elseif method == "NS"
        verbose && @info "Using Nuttall-Strand method"
        method_ar = (u, order) -> mc_ar_lc(u, order; method = method)
    elseif method == "VM"
        verbose && @info "Using Vieira-Morf method"
        method_ar = (u, order) -> mc_ar_lc(u, order; method = method)
    else
        throw(DomainError(method, "Invalid method, see documentation for valid methods."))
    end

    # assign information criterion ic
    if criterion == "AIC" # Akaike's Informaion Criterion (AIC)
        ic = (len, pf, nChannels, order) -> len * log(abs(det(pf))) + 2order * nChannels^2
    elseif criterion == "HQ" # Hannan-Quinn
        ic = (len, pf, nChannels, order) -> len * log(abs(det(pf))) + 2log(log(len)) * order * nChannels^2
    elseif criterion == "BIC" # Bayesian Information Criterion (Schwarz, 1978)
        ic = (len, pf, nChannels, order) -> len * log(abs(det(pf))) + log(len) * order * nChannels^2
    elseif criterion == "FPE" # Final prediction error (Akaike, 1970)
        ic = (len, pf, nChannels, order) -> log(abs(det(pf)) * ((len + nChannels * order + 1) / (len - nChannels * order - 1))^nChannels)
    elseif criterion === nothing
        verbose && @info "Order: $maxorder is fix."
        pf, A, ef = method_ar(u, maxorder)
        vic = nothing
        Vicv = nothing
        return (maxorder, pf, A, ef, vic, Vicv)
    else
        throw(DomainError(criterion, "Not a valid information criterion for model selection."))
    end

    # estimate correct order
    vicv = Inf
    order = 1
    pf, A, ef = method_ar(u, order)
    Vicv = zeros(maxorder)
    while order <= maxorder
        _pf, _A, _ef = method_ar(u, order)
        vic = ic(len, _pf, nChannels, order)
        Vicv[order] = vic
        if vic > vicv
            order -= 1
            verbose && @info "Order: $order was chosen by $criterion."
            break
        end
        pf, A, ef = _pf, _A, _ef
        vicv = vic
        order += 1
    end
    # B  = [ ] 
    # eb = [ ] 
    # pb = [ ]
    # return (order,pf,A,pb,B,ef,eb,vic,Vicv)
    return (order, pf, A, ef, vicv, Vicv)
end

# least squares based estimation of multichannel AR model
function mc_ar_ls(u, order)
    len, nChannels = size(u)
    # Build up regressor matrix
    Z = get_Z(u, order)

    # regression
    Gamma = Z * Z'
    U1 = Gamma \ Z

    SU = u' * u - u' * Z' * U1 * u
    pf = SU / (len - nChannels * order - 1)
    b = kron(U1, I(nChannels)) * reshape(u', 1, nChannels * len)'
    ef = reshape(reshape(u', len * nChannels, 1) - kron(Z', I(nChannels)) * b, len, nChannels)
    A = reshape(b, nChannels, nChannels, order)

    return pf, A, ef
end

# multichannel AR based on Nuttal-Strand or Vieira-Morf method
function mc_ar_lc(u, order::Int; method = "NS")
    len, nChannels = size(u)

    pf = u' * u         # Eq. (15.90)
    pb = copy(pf)       # Eq. (15.90)
    ef = copy(u)        # Eq. (15.91)
    eb = copy(u)        # Eq. (15.91)

    # Allocate forward and backward autoregressive coefficient matrices
    A = zeros(nChannels, nChannels, order)
    B = zeros(nChannels, nChannels, order)

    # Main Loop - Levinson Recursion 
    for M in 1:order
        #  Update estimated covariance errors             - Eq. (15.89)
        pfhat = ef[M+1:len, :]' * ef[M+1:len, :]
        pbhat = eb[M:len-1, :]' * eb[M:len-1, :]
        pfbhat = ef[M+1:len, :]' * eb[M:len-1, :]

        if method == "NS"
            #  Calculate estimated partial correlation matrix - Eq. (15.98)
            RHO = sylvester(pfhat * inv(pf), inv(pb) * pbhat, -2 * pfbhat)
            
            #  Update forward and backward reflection coefficients
            AM = -RHO * inv(pb)     # Eq. (15.73)
            BM = -RHO' * inv(pf)    # Eq. (15.74) & (15.78)
        end

        if method == "VM"
            # Eq. (15.88)
            Spfhat = pfhat^0.5
            Spbhat = pbhat^0.5
            ISpfhat = inv(Spfhat)
            ISpbhat = inv(Spbhat)
            RHO = ISpfhat * pfbhat * ISpbhat'
         
            AM = -Spfhat * RHO * ISpbhat    # Eq. (15.82)
            BM = -Spbhat * RHO' * ISpfhat   # Eq. (15.83)
        end
        A[:, :, M] = AM
        B[:, :, M] = BM

        #  Update forward and backward covariance error
        pf = pf - AM * BM * pf  # Eq. (15.75)
        pb = pb - BM * AM * pb  # Eq. (15.76)

        #  Update forward and backward predictor coefficients 
        for K in 1:M-1
            temp1 = A[:, :, K]
            A[:, :, K] = A[:, :, K] + AM * B[:, :, M-K]     # Eq. (15.71)
            B[:, :, M-K] = B[:, :, M-K] + BM * temp1        # Eq. (15.71)
        end

        Tef = copy(ef)
        ef[M+1:len, :] = ef[M+1:len, :] + (AM * eb[M:len-1, :]')'   # Eq.(15.84)
        eb[M+1:len, :] = eb[M:len-1, :] + (BM * Tef[M+1:len, :]')'  # Eq.(15.85
    end
    pf /= len # normalization 
    # return pf, -A, pb, -B, ef, eb #,ISTAT
    return pf, -A, ef
end