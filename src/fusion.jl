# Various routines relating specifically to fusion reactions
#TODO: Fix Ealpha, malpha, and valpha to exact values
#TODO: Get correct values for fusion cross section and Maxwellian sigma_v
module fusion
using ..constants
using ..species

const Ealpha = 3.5e6*el	# alpha birth energy
const malpha = 4.0*mp
const valpha = sqrt(2.0*Ealpha/(malpha))

"""
fusionSourceDT()

Calculates the fusion cross section for Maxwellian D and T.
"""
function fusionSourceDT(D::SpeciesData,T::SpeciesData)

   # Get these from standard sources
   sigma_v_avg = 1.0

   return D.density*T.density*sigma_v_avg
end
function fusionSourceDT(f_D::Function,f_T::Function)

   # Calculate directly from distribution function f_D and f_T
   return stuff
end


"""
DTcrossSection(v::Float64)

Calculates the fusion cross section sigma(v)
"""
function crossSectionDT(v::Float64)

end

export Ealpha,malpha,valpha, fusionSourceDT, crossSectionDT

end
