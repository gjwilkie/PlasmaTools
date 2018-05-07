# TODO: replace lnLambda with exact version. Use CODEs
# TODO: make tests
# TODO: Verify nus, nupar, etc. against H&S
# TODO: GEt rid of nuhatvta3 so that some temperature is defined
# TODO: populate solveCollisionalEquation function

# Contains various routines related to collision operators
module collisions
using ..constants
using ..species
using ..grids

export slowingdowntime, nu_s, nu_par, collisionOperator, lnLambda, nuhat, gfunc, nu_D

"""
nu_s(v::Real, a::SpeciesData, b::SpeciesData)

Returns νₛᵃᵇ(v) from Helander and Sigmar (2002), equation 3.46 (without the mass ratio factor that is cancelled in 3.40), where a and b are the test particle and target particle respectively
"""
function nu_s(v::Float64,a::SpeciesData,b::SpeciesData)
  vtb = sqrt(2.0 * b.temperature / b.mass)

  return nuhatvta3(a,b) * (a.mass/b.temperature) * Goverx(v/vtb)/vtb
end

function nu_s(v::Vector,a::SpeciesData,b::SpeciesData)
  vtb = sqrt(2.0 * b.temperature / b.mass)
  return nuhatvta3(a,b)* (a.mass/b.temperature) * Goverx(v/vtb)/vtb
end

"""
nu_par(v::Real, a::SpeciesData, b::SpeciesData)

Returns ν∥ᵃᵇ(v) from Helander and Sigmar (2002), equation 3.47, where a and b are the test particle and target particle respectively
"""
function nu_par(v::Float64,a::SpeciesData,b::SpeciesData)
  vtb = sqrt(2.0 * b.temperature / b.mass)
  if abs(v) .> eps(Float64)
     return 0.0
  else
     return nuhatvta3(a,b)* 2.0 * gfunc(v/vtb) / v^3
  end 
end

function nu_par(v::Vector,a::SpeciesData,b::SpeciesData)
  vtb = sqrt(2.0 * b.temperature / b.mass)
  y=zeros(length(v))
  idx=abs(v) .> eps(Float64)
  y[idx] = nuhatvta3(a,b)* 2.0 * gfunc(v[idx]/vtb) ./ v[idx].^3
  return y
end

function nu_D(v::Vector,a::SpeciesData,b::SpeciesData)
  vtb = sqrt(2.0 * b.temperature / b.mass)
  vta = sqrt(2.0 * a.temperature / a.mass)
  y=zeros(length(v))
  idx=abs(v) .> eps(Float64)
  y[idx] = nuhatvta3(a,b)* (erf(v[idx]/vtb) - gfunc(v[idx]/vtb)) ./ (v[idx]/vta).^3
  return y
end


"""
Chandresekhar G function
"""
function gfunc(x::Float64)
   if abs(x) > eps(Float64)
      return (erf(x)-(2.0/sqrt(pi))*x.*exp(-x.^2))./(2.0*x.^2)
   else
      return 0.0
   end
end
function gfunc(x::Vector)
   y=zeros(length(x))
   idx=abs(x) .> eps(Float64)
   y[idx] = (erf(x[idx])-(2.0/sqrt(pi))*x[idx].*exp(-x[idx].^2))./(2.0*x[idx].^2)
   return y
end

function Goverx(x::Float64)
   if abs(x) > eps(Float64)
      return (erf(x)-(2.0/sqrt(pi))*x.*exp(-x.^2))./(2.0*x.^3)
   else
      return 2.0/(3.0*sqrt(pi))    
   end
end
function Goverx(x::Vector)
   y=zeros(length(x))
   idx=abs(x) .> eps(Float64)
   y[idx] = (erf(x[idx])-(2.0/sqrt(pi))*x[idx].*exp(-x[idx].^2))./(2.0*x[idx].^3)
   y[abs(x) .<= eps(Float64)] = 2.0/(3.0*sqrt(pi)) 
   return y
end

"""
collisionOperator(vgrid::Grid, a::SpeciesData, b::SpeciesData)

Constructs a finite volume form of the collision operator
"""
function collisionOperator(vgrid::Grid, a::SpeciesData, b::SpeciesData)
   Nv = length(vgrid.ctr)

   if false
   # Replace a zero element with epsilon*next largest element
   vgrid[abs(vgrid).<eps(Float64)] = eps(Float64)*sort(abs(vgrid))[2]

   if ddv_in == [0.0]
      size(ddv_in) == (Nv,Nv) || error("ddv_in must be a square matrix the same size as v")
      ddv = ddv_in
   else
      ddv = createDifferentiationMatrix(vgrid)
   end

   return vgrid.^(-2) * ddv * ( nu_s(vgrid,a,b).*vgrid.^3  + 0.5*nu_par(vgrid,a,b).*vgrid.^4*ddv)
   end

   # Define functions that are required by operator constructor
   jacobian(v) = v.^2
   nu_s_term(v) = nu_s(v,a,b).*v
   nu_par_term(v) = 0.5*nu_par(v,a,b).*v.^2
   zerofunc(v) = 0.0
   # Since boundary conditions for collision operator are zero, we can worry about the source term later
   source = zeros(Nv)                        

   collop = createAdvecDiffOperator(vgrid,jacobian,nu_s_term,nu_par_term,[0.0,0.0],source)

   return collop
end
function collisionOperator(vgrid::Grid, a::SpeciesData, b::Array{SpeciesData,1})
   Nv = length(vgrid.ctr)
   collop = spzeros(Nv,Nv)
   for is in 1:length(b)
      collop = collop + collisionOperator(vgrid,a,b[is])
   end
   return collop
end

"""
solveCollisionalEquation(a::SpeciesData,b::Array{SpeciesData,1};rhs::Array{Float64,1}=[0.0];extraDrag=[0.0];extraDiff=[0.0])

Solves equations of the form C[F] = S with additional terms modifying collision operator accordingly. Subject to constraint that the solution integrates to a.dens.
"""
function solveCollisionalEquation(a::SpeciesData,b::Array{SpeciesData,1};rhs::Array{Float64,1}=[0.0],extraDrag::Array{Float64,1}=[0.0],extraDiff::Array{Float64,1}=[0.0])

end

"""
nuhat*v_ta^3. Useful for non-Maxwellian species, where v_ta is not defined.
"""
function nuhatvta3(a::SpeciesData,b::SpeciesData)
  return b.density * b.charge^2 * a.charge^2 * lnLambda(a,b) / (4.0*pi*ep0^2 * a.mass^2)
end

"""
nuhat from Helander and Sigmar.
"""
function nuhat(a::SpeciesData,b::SpeciesData)
   return nuhatvta3(a,b) / (2.0*a.temperature/a.mass)^1.5
end

"""
lnLambda(a::SpeciesData, b::SpeciesData)

Calculates the Coulomb logarithm ln(Λᵃᵇ ) using approximations in NRL plasma formulary"
"""
function lnLambda(s1::SpeciesData, s2::SpeciesData)
   local me_threshold = 0.05*mp	# Maximum mass that can be considered an "electron", in case one needs reduced mass ratio

   # Check for electron-electron collisions
   if (s1.mass < me_threshold) && (s2.mass < me_threshold)
      # Convert to eV and cm^3
      Te = s1.temperature/el 	
      ne = s1.density*1.0e-6
      return 23.5 - log(sqrt(ne)/Te^1.25) - sqrt(1.0e-5 + (log(Te)-2.0)^2/16.0)
   # Check for electron-ion collisions
   elseif (s1.mass < me_threshold) || (s2.mass < me_threshold)
      if (s1.mass < me_threshold) 
         me = s1.mass; Te = s1.temperature/el; ne = 1.e-6*s1.density
         mi = s2.mass; Ti = s2.temperature/el; ni = 1.e-6*s2.density; Z=s2.charge/el
      else 
         mi = s1.mass; Ti = s1.temperature/el; ni = 1.e-6*s1.density; Z=s1.charge/el
         me = s2.mass; Te = s2.temperature/el; ne = 1.e-6*s2.density
      end 
      if Te < Ti*Z*me/mi
         return 30.0 - log(sqrt(ni)*Z^2*(mp/mi)/Ti^1.5)
      elseif Te > 10.0*Z^2
         return 24.0 - log(sqrt(ne)/Te) 
      else
         return 23.0 - log(sqrt(ne)*Z/Te^1.5)
      end
  # Ion-ion collision
  else
    Z1 = s1.charge/el; Z2 = s2.charge/el
    m1 = s1.mass/mp; m2 = s2.mass/mp
    n1 = s1.density*1.0e-6; n2 = s2.density*1.0e-6
    T1 = s1.temperature/el; T2 = s2.temperature/el
    return 23.0 - log( abs(Z1*Z2)* (m1+m2)*sqrt( (n1*Z1^2/T1) + (n2*Z2^2/T2)) /(m1*T2 + m2*T1 ))
  end
end

"""
slowingDownTime(a::SpeciesData, el::SpeciesData)

Calculates the slowing-down time of species a against electrons where v_ta << v_te.
"""
function slowingdowntime(a::SpeciesData, electron::SpeciesData)
   num = (3.0/(16.0*sqrt(pi)))*electron.mass*a.mass*(2.0*electron.temperature/electron.mass)^1.5
   return num/(a.charge*el^2*electron.density*lnLambda(a,electron))
end

end

