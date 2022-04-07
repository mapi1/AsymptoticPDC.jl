
# """
# function detrend(signal::AbstractVector; p::Int = 1, coefs::Union{AbstractVector,Nothing,Real} = nothing)
# Detrend a signal by removing polynomial trend of order p using build in least squares.
# Choose p = 0 to remove only mean or input coefficients from previous detrending to detrend by those.
# # Args:
# *   signal::Vector: Data Vector containing te signal
# *   p::Int: order of polynomial
# *   coefs::Union{Vector{<:Real}, Nothing, Real}: Coefficients to do the same detrending on different signal
# *   return_coefs::Bool: if true returns coefficients detrended by
# # Return:
# * (newSignal, coefs): The detrended signal and the estimated coefficients from order 0 to  p
# # Examples
# ```julia
# julia> signal = sin.([1:100;]) + 0.03 .* [1:100;]
# julia> detrend(signal)
# Vector{Float}
# ```
# """
# function detrend(signal::AbstractVector; p::Int=1, coefs::Union{AbstractVector,Nothing,Real}=nothing)
#     p >= 0 || throw(DomainError("Order p hast to be positive or zero, is: $p"))
#     # Build up regression matrix x
#     x = ones(length(signal))
#     xt = cumsum(x)
#     for i in 1:p
#         x = hcat(x, xt .^ i)
#     end
#     # If no coefs given estimate them
#     if coefs === nothing
#         coefs = x \ signal
#     else
#         length(coefs) == p + 1 || throw(DomainError("p is $p, but length of coefs is $(length(coefs))"))
#     end
#     trend = x * coefs

#     return (signal .- trend), coefs
# end
# detrend(x, 1) # detrend order one
# detrend(x, [1]) # detrend coeff one

function detrend(signal::AbstractVector, order::Int=1; verbose::Bool=true)
    order >= 0 || throw(DomainError("Order p hast to be positive or zero, is: $order"))
    verbose && @info "polynomial detrending of order $order"
    
    # Build up regression matrix x
    x = ones(length(signal))
    xt = cumsum(x)
    for i in 1:order
        x = hcat(x, xt .^ i) # Risk for Int overflow for large orders and 
    end

    coefs = x \ signal # Estimate coefficients using LS
    trend = x * coefs

    return (signal .- trend), coefs
end

function detrend(signal::AbstractVector, coeffs::AbstractVector; verbose::Bool=true)
    verbose && @info "polynomial detrending with given coefficients"
    order = length(coeffs) - 1
    # Build up regression matrix x
    x = ones(length(signal))
    xt = cumsum(x)
    for i in 1:order
        x = hcat(x, xt .^ i)
    end

    trend = x * coeffs

    return (signal .- trend), coeffs
end

function detrend!(signal::AbstractVector; p::Int=1)

end
function retrend()
    #TODO
end

function get_Z(u, order)
    len, nChannels = size(u)
    Z = zeros(eltype(u), nChannels * order, len)
    for i in 1:order
        inds = (1:nChannels) .+ nChannels * (i - 1)
        Z[inds, i + 1 : end] = u'[:, 1 : end - i]
    end
    return Z
end

# function dmatrix(m)
#     D = zeros(m^2, Int(m * (m + 1) / 2))
#     u = hcat([[i, j] for i in 1:m, j in 1:m]...)'
#     w = hcat([[j, i] for i in 1:m, j in 1:m if i >= j]...)'
#     v = reverse(w, dims = 2)
#     for i=1:m*m
#         for j=1:Int(m*(m+1)/2)
#             if sum(u[i,:] .== v[j,:]) == 2
#                 D[i,j]=1
#             end
#         end
#         for j=1:Int(m*(m+1)/2)
#             if sum(u[i,:] .== w[j,:]) == 2
#                 D[i,j]=1
#             end
#         end
#     end
#     return D
# end

function dmatrix(n)
    D = zeros(n * n, Int((n * (n + 1)) / 2))
    count = 1
    for j = 1:n, i = 1:n
        if i >= j
            D[(j-1)*n+i, count] = 1
            count = count + 1
        else
            D[(j-1)*n+i, :] = D[(i-1)*n+j, :]
        end
    end
    return D
end

function vech(Y)
    m, n = size(Y);
    y = eltype(Y)[]
    for i=1:m
        push!(y, Y[i:n,i]...)
    end
    return y
end

function fIij(i::Int, j::Int, nChannels::Int)
    Iij = spzeros(nChannels^2)
    Iij[nChannels*(j-1)+i] = 1
    Iij = spdiagm(Iij)
    c = kron(I(2), Iij)
    return c
end

function fIj(j::Int, nChannels)
    Ij = spzeros(nChannels)
    Ij[j] = 1
    Ij = spdiagm(Ij)
    Ij = kron(Ij, I(nChannels))
    c = kron(I(2), Ij)
    return c
end


# Helper converts A into the frequency domain
function A2f(A, nFreqs::Int)
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



# Helper that tries Cholesky factorization and diagonalizes otherwise
function fChol(omega)
    L = 0
    try
        L = cholesky(omega).L
        # If there's a small negative eigenvalue, diagonalize
    catch e
        println(e)
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

function fCa(f, order::Int, nChannels::Int)
    C1 = cos.(-2π * f .* (1:order))
    S1 = sin.(-2π * f .* (1:order))
    C2 = [C1 S1]'
    d = kron(sparse(C2), I(nChannels^2))
    return d
end

# Helper 
# Autocorrelation. Data in row-wise orientation. From order 0 to p-1.
# Output: n x n blocks of autocorr of lags i. (Nuttall Strand matrix)
function bigautocorr(u, order)
    len, nChannels = size(u)
    gamma = zeros(nChannels * order, nChannels * order)
    @inbounds for i = 1:order, j = 1:order
        gamma[((i-1)*nChannels+1):i*nChannels, ((j-1)*nChannels+1):j*nChannels] = lag(u, i - 1, default=zero(u[1]))' * lag(u, j - 1, default=zero(u[1]))
    end
    return gamma ./ len
end

function fEig(L, G2)
    D = L' * G2 * L
    d = svd(D).S
    # d = eigvals!(D, 7:8)
    return sort(d)[end-1:end]
end

function TT(a,b)
    t = spzeros(a*b, a*b)
    for i in 1:a, j in 1:b
        t[(i-1)*b+j,(j-1)*a+i] = 1
    end
    return t
end

function fdebig_de(n)
    A = kron(TT(2n, n), sparse(I(2n*n)))
    B = kron(vec(sparse(I(2*n))), sparse(I(n)))
    c = A * kron(sparse(I(n)), B)
    return c
end