

"""
pdcplot(apdc<:AbstractPartialDirectedCoherence)
"""
@userplot PDCplot

@recipe function f(h::PDCplot; f_range = (0, 0.5), cnames = nothing)
    typeof(h.args[1]) <: AbstractPartialDirectedCoherence || throw(DomainError("Input must be <: AbstractPartialDirectedCoherence"))
    
    data = h.args[1]
    coherence = data.coherence 
    spectra = data.spectra
    f = data.freq
    with_asymptotic = typeof(data) <: AsymptoticPartialDirectedCoherence
    
    nChannels = size(coherence, 1)
    
    if cnames === nothing
        cnames = "x" .* string.(collect(1:nChannels))
    else
        length(cnames) == nChannels || throw(ArgumentError("There needs to be a name for each channel."))
    end

    # plot design
    layout := (nChannels, nChannels)
    
    for j in 1:nChannels, i in 1:nChannels
        if i == 1
            title := "$(cnames[j])"
        else
            title := ""
        end
        if i != j 
            @series begin
                if with_asymptotic
                    @series begin
                        subplot := ((i-1)*nChannels + j)
                        seriestype := :vspan
                        seriescolor := :grey
                        linewidth := 0
                        alpha := 0.5
                        label := ""
                        significant = [false; data.pvalues[i,j,:] .< data.α]
                        inds = findall(x -> !iszero(x), diff(significant))
                        f[inds]
    
                    end
                end

                seriescolor := :orangered2
                subplot := ((i-1)*nChannels + j)
                label := "$(cnames[i]) <- $(cnames[j])" 
                ylims := (0,1)
                if with_asymptotic
                    ribbon := (abs.(data.lower_conf[i,j,:] .- coherence[i,j,:]), data.upper_conf[i,j,:] .- coherence[i,j,:])
                end
                if i == nChannels
                    xlabel := "f"
                end
                
                f, vec(coherence[i,j,:])
            end
        else
            @series begin
                subplot := ((i-1)*nChannels + j)
                label := "S $(cnames[i])"
                seriescolor := :black
                linewidth := 2
                # if j == 1 
                #     ylabel := cnames[i]
                # end
                if i == nChannels
                    xlabel := "f"
                end
                
                f, vec(2abs.(spectra[i,j,:]))
            end
        end
    end
end

@userplot parPSDplot

@recipe function f(h::parPSDplot; sf = nothing, f_range = (0, 0.5), cnames = nothing)
    # input processing
    length(h.args) == 2 || throw(ArgumentError("pdcplot needs pdc and spectra as input"))
    length(f_range) == 2 || throw(ArgumentError("if the frequency range shall be specified a start and stop frequency are needed: (start, stop)"))
    pdc, spectra = h.args
    size(pdc) == size(spectra) || throw(DomainError("pdc and spectra have not the same size, are thex from the same data?"))
    
    nChannels = size(pdc, 1)
    nFreqs = size(pdc, 3)
    
    if cnames === nothing
        cnames = "x" .* string.(collect(1:nChannels))
    else
        length(cnames) == nChannels || throw(ArgumentError("There needs to be a name for each channel."))
    end
    
    f = range(f_range[1], stop = f_range[2], length = nFreqs)
    if sf === nothing
        f_flag = "c/b"
    else
        f = f .* sf
        f_flag = "Hz"
    end
    # plot design
    # title := "PDC"
    layout := (nChannels, nChannels)
    
    for j in 1:nChannels, i in 1:nChannels
        if i != j
            @series begin
                subplot := ((i-1)*nChannels + j)
                label := "$(cnames[i]) <- $(cnames[j])" 
                # seriescolor := :red
                if j == 1 
                    ylabel := cnames[i]
                end
                if i == nChannels
                    xlabel := "f [$f_flag]"
                end
                
                f, vec(pdc[i,j,:]) .* vec(2abs.(spectra[i,i,:]))
            end
        else
            @series begin
                subplot := ((i-1)*nChannels + j)
                label := "S$(cnames[i])(f)"
                seriescolor := :black
                linewidth := 2
                if j == 1 
                    ylabel := cnames[i]
                end
                if i == nChannels
                    xlabel := "f [$f_flag]"
                end
                
                f, vec(2abs.(spectra[i,j,:]))
            end
        end
    end
end