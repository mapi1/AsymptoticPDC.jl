module PDC

using DSP
using LinearAlgebra
using ShiftedArrays
using RecipesBase

include("pdc_algorithm.jl")
include("mvar.jl")
include("plotting_recipes.jl")
include("utility.jl")
include("../examples/sunspots_melanoma.jl")

export pdc
export mvar
export pdcplot, pdcplot!
export detrend, detrend!
export get_sunspot_melanoma_data

end # module
