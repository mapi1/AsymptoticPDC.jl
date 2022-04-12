"""
    mcar(u; maxorder::Union{Nothing,Int}=nothing, criterion::Union{Nothing,String}="AIC", method::String="NS", verbose::Bool=true)

Compute a multichannel AR or vector AR model of the input matrix u containing the signals/channels xi, u = [x1 x2 ... xn] 

# Args

* u: input Matrix containing signals 1 to n `u = [x1 x2 ... xn]`

# Keywords

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
Returns a tuple (model, best_ic_value, ic_values)

Result model is of type `MCAR_Model`, with following fields:
* order: is the (chosen) model order
* nChannels: number of channels
* samples: number of samples per channel
* A: contains the AR coefficients [n x n x order]
* pf: is the covariance matrix [order x order]
* ef: the residuals
"""
function mcar(u; maxorder::Union{Nothing,Int}=nothing, criterion::Union{Nothing,String}="AIC", method::String="NS", verbose::Bool=true)
    samples, nChannels = size(u)
    verbose && @info "$nChannels channels with $samples samples."

    if maxorder === nothing
        maxorder = Int(round(3 * sqrt(samples) / nChannels)) # Suggested by Nuttall, 1976.
        verbose && @info "The maxorder was set to $maxorder."
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
        return (MCAR_Model(maxorder, nChannels, samples, A, pf, ef), nothing, nothing)
    else
        throw(DomainError(criterion, "Not a valid information criterion for model selection."))
    end

    # estimate correct order
    best_ic_value = Inf
    order = 0
    ic_values = zeros(maxorder)
    while order < maxorder
        order += 1
        pf, A, ef = method_ar(u, order)
        current_ic = ic(samples, pf, nChannels, order)
        ic_values[order] = current_ic
        if current_ic > best_ic_value
            order -= 1
            break
        end
        best_ic_value = current_ic
    end
    verbose && @info "Order: $order was selected by $criterion."

    return (MCAR_Model(order, nChannels, samples, A, pf, ef), best_ic_value, ic_values)
end

"""
`MCAR_Model` with following fields:
* order: is the (chosen) model order
* nChannels: number of channels
* samples: number of samples per channel
* A: contains the AR coefficients [n x n x order]
* pf: is the covariance matrix [order x order]
* ef: the residuals    
"""
struct MCAR_Model{T, U<:AbstractArray{T,3}, V<:AbstractArray{T,2}}
    order::Int
    nChannels::Int
    samples::Int
    A::U
    pf::V
    ef::V
end

# least squares based estimation of multichannel AR model
function mc_ar_ls(u, order)
    samples, nChannels = size(u)
    # Build up regressor matrix
    Z = get_Z(u, order)

    # regression
    Gamma = Z * Z'
    U1 = Gamma \ Z

    SU = u' * u - u' * Z' * U1 * u
    pf = SU / (samples - nChannels * order - 1)
    b = kron(U1, I(nChannels)) * reshape(u', 1, nChannels * samples)'
    ef = reshape(reshape(u', samples * nChannels, 1) - kron(Z', I(nChannels)) * b, samples, nChannels)
    A = reshape(b, nChannels, nChannels, order)

    return pf, A, ef
end

# multichannel AR based on Nuttal-Strand or Vieira-Morf method
function mc_ar_lc(u, order::Int; method = "NS")
    samples, nChannels = size(u)

    input_type = !(eltype(u) <: LinearAlgebra.BlasFloat) ? Float64 : eltype(u) # sylvester() only accepts BlasFloat
    pf = u' * u         # Eq. (15.90)
    pb = copy(pf)       # Eq. (15.90)
    ef = input_type.(u)        # Eq. (15.91)
    eb = input_type.(u)        # Eq. (15.91)

    # Allocate forward and backward autoregressive coefficient matrices
    A = zeros(input_type, nChannels, nChannels, order)
    B = zeros(input_type, nChannels, nChannels, order)

    # Main Loop - Levinson Recursion 
    for M in 1:order
        #  Update estimated covariance errors             - Eq. (15.89)
        pfhat = ef[M+1:samples, :]' * ef[M+1:samples, :]
        pbhat = eb[M:samples-1, :]' * eb[M:samples-1, :]
        pfbhat = ef[M+1:samples, :]' * eb[M:samples-1, :]

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
        ef[M+1:samples, :] = ef[M+1:samples, :] + (AM * eb[M:samples-1, :]')'   # Eq.(15.84)
        eb[M+1:samples, :] = eb[M:samples-1, :] + (BM * Tef[M+1:samples, :]')'  # Eq.(15.85
    end
    pf /= samples # normalization 
    # return pf, -A, pb, -B, ef, eb #,ISTAT
    return pf, -A, ef
end