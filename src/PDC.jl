module PDC

# using DSP
using LinearAlgebra
using ShiftedArrays
using RecipesBase
using Distributions

include("pdc_algorithm.jl")
include("mvar.jl")
include("plotting_recipes.jl")
include("granger_cuasality_test.jl")
include("utility.jl")
include("../examples/sunspots_melanoma.jl")

export pdc
export mvar
export pdcplot, pdcplot!
export detrend, detrend!
export get_sunspot_melanoma_data
export granger_causality_test, instantaneous_granger_causality_test
end # module
