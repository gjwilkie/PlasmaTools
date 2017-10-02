#!/bin/julia
using NetCDF
using PyPlot
using Dierckx
using ..constants:mp,me
using ..gs2ls: GS2specs, GS2geometry


"""
plotSummary(file::String)

Displays a set of time-dependent plots of relevant plasma quantities. Useful for determining the overall health of the run and ideal time samples.
"""
function plotTranspSummary(file::String)

   t = ncread(file,"TIME")
   Ti0 = ncread(file,"TI0")
   Te0 = ncread(file,"TE0")
   q95 = ncread(file,"Q95")
   ne0 = ncread(file,"NE")[1,:]'

   #TODO: Add some kind of global diagnostic of fast ions


end

"""
plotSummary(file::String,time::Float64)

Displays a set of radius-dependent plots of relevant plasma quantities at the given time.
"""
function plotTranspSummary(file::String,time::Float64)

   # TODO: Make radial plots
   kappa = ncread(file,"ELONG")[1,:]
   tri = ncread(file,"TRIANG")[1,:]
   
end

function getValFrom2Darray(rho::Float64,time::Float64,rgrid::Vector,tgrid::Vector,y::Matrix)
   yfunc = Spline2D(rgrid,tgrid,y,bc="nearest")
   return evaluate(yfunc,rho,time)
end

function getGradient(rho::Float64,time::Float64,rgrid::Vector,tgrid::Vector,y::Matrix)
   dr = 0.01
   yfunc = Spline2D(rgrid,tgrid,y,bc="nearest")
   y2 = evaluate(yfunc,rho+dr,time)
   y1 = evaluate(yfunc,rho-dr,time)
   return (y2-y1)/(2.0*dr)
end

function getArrayAtTime(t::Float64,tgrid::Vector,ygrid::Vector)
   nx = length(ygrid[:,1])
   x = zeros(Float64,nx)
   for i in 1:nx
      yfunc = Spline1D(tgrid,y[i,:]',bc="nearest")
      x[i]  = evaluate(yfunc,t)
   end
   return x
end


function getLocalData(file::String,rho::Float64,t0::Float64,Zi::Float64,Ai::Float64;Zz::Float64=-999.9,Az::Float64=-999.9)

   t = ncread(file,"TIME")
   x = getArrayAtTime(t0,t,ncread(file,"X"))
   xb = getArrayAtTime(t0,t,ncread(file,"XB"))

   beta_e = getValFrom2Darray(rho,t0,x,t,ncread(file,"BTE"))

   if (Zz == -999.9) && (Az == -999.9)
      include_impurity = false
   else
      include_impurity = true
   end

# Params and species:
# To return: beta_e, [Ti, Te, Tz, ni, ne, nz, + gradients]
# Corresponding: BTE, NI, NE, XDENS, TI, TE, TX

# Geometry:
# qinp, shat, [r_geo, shift], rmaj, akappa,akappri,tri,tripri,
# Q, SHAT, [RMJMP,gradient], RAXIS, [ELONG, TRIANG + gradients]


# Other relevant quantities:
# BDENS: beam ion density. How does this affect quasineutrality

   # TODO: Define nuii, nuee, nuzz

   specs=GS2species[]
   push!(specs,GS2species(dens=ni/ne,temp=1.0,z=Zi,mass=1.0,vnewk=nuii*a/vti,fprim=gradni*a,trpim=gradTi*a))
   if include_impurity
      push!(specs,GS2species(dens=nz/ne,temp=Tz/Ti,z=Zimp,mass=Aimp/Ai,vnewk=nuzz*a/vti,fprim=gradnz*a,trpim=gradTz*a))
   end
   push!(specs,GS2species(dens=1.0,temp=Te/Ti,z=-1.0,mass=me/(Ai*mp),vnewk=nuee*a/vti,fprim=gradne*a,trpim=gradTe*a))

   geom = GS2geometry(equilibrium_option="eik",rhoc=rho,irho=1,qinp=q,shat=shat, 
                       shift=Rprime,akappa=kappa,akappri=kappri,tri=tri,tripri=tripri)

   return beta_e, geom, specs
end
