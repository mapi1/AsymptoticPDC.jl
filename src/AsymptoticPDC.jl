module AsymptoticPDC

using DSP
using LinearAlgebra
using ShiftedArrays
using RecipesBase
using Distributions
using SparseArrays

include("mvar.jl")
include("pdc_algorithms.jl")
include("spectral_algorithms.jl")
include("plotting_recipes.jl")
include("granger_causality_test.jl")
include("utility.jl")
include("../examples/sunspots_melanoma.jl")

export pdc
export mvar, MCAR_Model
export pdcplot, pdcplot!
export detrend, detrend!
export get_sunspot_melanoma_data
export granger_causality_test, instantaneous_granger_causality_test
export spectra
export coherence

end # module
