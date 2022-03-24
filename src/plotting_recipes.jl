

@userplot PDCplot

@recipe function f(h::PDCplot; sf = nothing, f_range = (0, 0.5), cnames = nothing)
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
                ylims := (0,1)
                if j == 1 
                    ylabel := cnames[i]
                end
                if i == nChannels
                    xlabel := "f [$f_flag]"
                end
                
                f, vec(pdc[i,j,:])
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