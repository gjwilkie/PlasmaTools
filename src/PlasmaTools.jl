# TODO: add all exports
module PlasmaTools

include("constants.jl")
include("species.jl")
include("fusion.jl")
include("grids.jl")
include("collisions.jl")
include("gs2ls.jl")
#include("gs2reader.jl")

using .constants
using .species
using .fusion
using .grids
using .collisions
using .gs2ls
#using .gs2plot

export el,me,mp,kB,ep0
export SpeciesData
export Grid,createFVgrid,createAdvecDiffOperator
export slowingDownTime, nu_s, nu_par, collisionOperator, lnLambda, nuhat, gfunc
export createInputFile, GS2species, GS2resolution, GS2params, GS2geometry

end
