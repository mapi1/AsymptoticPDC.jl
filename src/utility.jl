
"""
function detrend(signal::AbstractVector; p::Int = 1, coefs::Union{AbstractVector,Nothing,Real} = nothing)
Detrend a signal by removing polynomial trend of order p using build in least squares.
Choose p = 0 to remove only mean or input coefficients from previous detrending to detrend by those.
# Args:
*   signal::Vector: Data Vector containing te signal
*   p::Int: order of polynomial
*   coefs::Union{Vector{<:Real}, Nothing, Real}: Coefficients to do the same detrending on different signal
*   return_coefs::Bool: if true returns coefficients detrended by
# Return:
* (newSignal, coefs): The detrended signal and the estimated coefficients from order 0 to  p
# Examples
```julia
julia> signal = sin.([1:100;]) + 0.03 .* [1:100;]
julia> detrend(signal)
Vector{Float}
```
"""
function detrend(signal::AbstractVector; p::Int=1, coefs::Union{AbstractVector,Nothing,Real}=nothing)
    p >= 0 || throw(DomainError("Order p hast to be positive or zero, is: $p"))
    # Build up regression matrix x
    x = ones(length(signal))
    xt = cumsum(x)
    for i in 1:p
        x = hcat(x, xt .^ i)
    end
    # If no coefs given estimate them
    if coefs === nothing
        coefs = x \ signal
    else
        length(coefs) == p + 1 || throw(DomainError("p is $p, but length of coefs is $(length(coefs))"))
    end
    trend = x * coefs

    return (signal .- trend), coefs
end
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