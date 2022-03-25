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
    - `NS` Nuttall Strand

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
        method_ar = (u, order) -> mvarls(u, order)
    elseif method == "NS"
        verbose && @info "Using method by Nuttall Strand"
        method_ar = (u, order) -> mvarns(u, order)
    else
        throw(DomainError(method, "Invalid method, see docummentation for vaild methods."))
    end

    # assign information criterion ic
    if criterion == "AIC" # Akaike's Informaion Criterion (AIC)
        ic = (len, pf, nChannels, order) -> len * log(abs(det(pf))) + 2order * nChannels^2
    elseif criterion == "HQ" # Hannan-Quinn
        ic = (len, pf, nChannels, order) -> len * log(abs(det(pf))) + 2log(log(len)) * order * nChannels^2
    elseif criterion == "BIC" # Bayesian Information Criterion, Schwarz 1978
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

# least squares based estimation of mvar model
function mvarls(u, order)
    len, nChannels = size(u)
    # Build up regressor matrix
    Z = zeros(nChannels * order, len)
    for i in 1:order
        inds = (1:nChannels) .+ nChannels * (i - 1)
        Z[inds, i + 1 : end] = u'[:, 1 : end - i]
    end

    # regression
    Gamma = Z * Z'
    U1 = Gamma \ Z

    SU = u' * u - u' * Z' * U1 * u
    pf = SU / (len - nChannels * order - 1)
    b = kron(U1, I(nChannels)) * reshape(u', 1, nChannels * len)'
    e = reshape(reshape(u', len * nChannels, 1) - kron(Z', I(nChannels)) * b, len, nChannels)
    # e = reshape(reshape(u, len * nChannels, 1) - kron(Z', I(nChannels)) * b, len, nChannels)
    # e = reshape(u, len * nChannels, 1) - kron(Z', Matrix(I, nChannels, nChannels)) * b

    A = reshape(b, nChannels, nChannels, order)
    ef = e

    return pf, A, ef
end

# mv ar, based on Nuttal Strand Algorithm
function mvarns(u::Matrix{<:Real}, order::Int)
    len, nChannels = size(u)

    ef = copy(u)        # Eq. (15.91)
    eb = copy(u)        # Eq. (15.91)
    pf = u' * u         # Eq. (15.90)
    pb = copy(pf)       # Eq. (15.90)
    # return pf,pb,ef,eb#,ISTAT

    M = 0
    #   ** Main Loop **
    A = zeros(nChannels, nChannels, order)
    B = zeros(nChannels, nChannels, order)
    while true
        #  Update estimated covariance errors               Eq. (15.89)
        pfhat = ef[M+2:len, :]' * ef[M+2:len, :]
        pbhat = eb[M+1:len-1, :]' * eb[M+1:len-1, :]
        pfbhat = ef[M+2:len, :]' * eb[M+1:len-1, :]

        M = M + 1
        #  Calculate estimated partial correlation matrix - Eq. (15.98)
        #             (Nuttall-Strand algorithm only)
        RHO = sylvester(pfhat * inv(pf), inv(pb) * pbhat, -2 * pfbhat)

        #  Update forward and backward reflection coeficients
        #  Eqs. (15.73),(15.74),(15.78) (algoritmo de Nuttall-Strand)
        AM = -RHO * inv(pb)
        BM = -RHO' * inv(pf)
        #   if exist('A')~=2, A=[]; end;
        A[:, :, M] = AM
        #   if exist('B')~=2, B=[]; end;
        B[:, :, M] = BM

        #  Update forward and backward covariance error  - Eqs. (15.75),(15.76)
        pf = pf - AM * BM * pf
        pb = pb - BM * AM * pb

        #  Update forward and backward predictor coefficients - Eqs.(15.84),(15.85)
        if M != 1
            for K in 1:M-1
                temp1 = A[:, :, K]
                A[:, :, K] = A[:, :, K] + AM * B[:, :, M-K]
                B[:, :, M-K] = B[:, :, M-K] + BM * temp1
            end
        end
        #  Update residues
        Tef = copy(ef)
        # ef[:,len:-1:M+1] = ef[:,len:-1:M+1] + AM * eb[:,len-1:-1:M];
        # ef[M+1:len,:] = ef[M+1:len,:] + eb[M:len-1,:] * AM;
        ef[M+1:len, :] = ef[M+1:len, :] + (AM * eb[M:len-1, :]')'
        # eb[:,len:-1:M+1] = eb[:,len-1:-1:M] + BM * Tef[:,len:-1:M+1];
        eb[M+1:len, :] = eb[M:len-1, :] + (BM * Tef[M+1:len, :]')'
        # eb[M+1:len,:] = eb[M:len-1,:] + Tef[M+1:len,:] * BM;

        #  Verify if model order is adequate
        if M == order
            A = -A
            B = -B
            break
        end
    end
    pf /= len # normalization 
    # return pf, A, pb, B, ef, eb #,ISTAT
    return pf, A, ef
end

function fIij(i::Int, j::Int, nChannels::Int)
    Iij = zeros(nChannels^2)
    Iij[nChannels*(j-1)+i] = 1
    Iij = diagm(Iij)
    c = kron(I(2), Iij)
    return c
end

function fIj(j::Int, nChannels)
    Ij = zeros(nChannels)
    Ij[j] = 1
    Ij = diagm(Ij)
    Ij = kron(Ij, Matrix(I, nChannels, nChannels))
    c = kron(I(2), Ij)
    return c
end


# Helper converts A into the frequency domain
function A2f(A::AbstractArray{<:Real}, nFreqs::Int)
    nChannels = size(A, 1)
    order = size(A, 3)

    exponents = -reshape((-im * π * kron(0:(nFreqs-1), (1:order)) / nFreqs), order, nFreqs)'
    Areshaped = reshape(A, nChannels, nChannels, 1, order)
    Af = zeros(Complex, nChannels, nChannels, nFreqs, order)

    for kk = 1:nFreqs
        Af[:, :, kk, :] = Areshaped
    end
    for i = 1:nChannels
        for k = 1:nChannels
            Af[i, k, :, :] = reshape(Af[i, k, :, :], nFreqs, order) .* exp.(exponents)
        end
    end
    Af = permutedims(Af, [3, 1, 2, 4])
    AL = zeros(Complex, nFreqs, nChannels, nChannels)

    for kk = 1:nFreqs
        temp = zeros(nChannels, nChannels)
        for k = 1:order
            temp = temp + reshape(Af[kk, :, :, k], nChannels, nChannels)
        end
        temp = I - temp
        AL[kk, :, :] = reshape(temp, 1, nChannels, nChannels)
    end
    return AL
end



# Helper that tries Cholesky factorization an diagonalizes otherwise
function fChol(omega)
    L = 0
    try
        L = cholesky(omega).L
        # If there's a small negative eigenvalue, diagonalize
    catch e
        println(e)
        #   disp('linalgerror, probably IP = 1.')
        ei = eigen(omega)
        vals, vecs = ei
        vecs = reverse(vecs, dims=2)
        reverse!(vals)
        L = zeros(size(vecs))
        for i = 1:length(vals)
            if ei.values[i] < 0
                ei.values[i] = eps()
            end
            L[:, i] = vecs[:, i] .* sqrt(vals[i])
        end
    end
    return L
end

function fCa(f::Real, order::Int, nChannels::Int)
    C1 = cos.(-2π * f .* (1:order))
    S1 = sin.(-2π * f .* (1:order))
    C2 = [C1 S1]'
    d = kron(C2, Matrix(I, nChannels^2, nChannels^2))
    return d
end

function Dup(n)
    d = zeros(n * n, Int((n * (n + 1)) / 2))
    count = 1
    for j = 1:n, i = 1:n
        if i >= j
            d[(j-1)*n+i, count] = 1
            count = count + 1
        else
            d[(j-1)*n+i, :] = d[(i-1)*n+j, :]
        end
    end
    return d
end

# Helper 
# Autocorrelation. Data in row-wise orientation. From order 0 to p-1.
# Output: n x n blocks of autocorr of lags i. (Nuttall Strand matrix)
function bigautocorr(u, order)
    len, nChannels = size(u)
    gamma = zeros(nChannels * order, nChannels * order)
    for i = 1:order, j = 1:order
        gamma[((i-1)*nChannels+1):i*nChannels, ((j-1)*nChannels+1):j*nChannels] = lag(u, i - 1, default=zero(u[1]))' * lag(u, j - 1, default=zero(u[1]))
    end
    return gamma ./ len
end

