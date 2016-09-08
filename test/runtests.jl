#!/usr/bin/env julia
include("../src/PlasmaTools.jl")
using PlasmaTools
using Base.Test


# Resolution for tests
Nv = 800

# Test conditions
n_cgs = 1.0e14
T_eV = 10.0e3

T = T_eV*el
n = n_cgs*1.0e6

deuterium = SpeciesData(density=0.5*n,temperature=T,charge=el,mass=2.0*mp)
tritium = SpeciesData(density=0.5*n,temperature=T,charge=el,mass=3.0*mp)
electron = SpeciesData(density=n,temperature=T,charge=-el,mass=me)
helium3 = SpeciesData(density=0.01*n,temperature=T,charge=2.0*el,mass=3.0*mp)

bulkspecs = [deuterium,tritium,electron]

vti = sqrt(2.0*T/mp)
vte = sqrt(2.0*T/me)
vmax = 5.0*vti
vgrid = createFVgrid(Nv,0.0,vmax,[1,0])

f0trace = helium3.density*(2.0*pi*helium3.temperature/helium3.mass)^(-1.5)*exp(-0.5*helium3.mass*vgrid.ctr.^2/helium3.temperature)

# Test value of lnLambda, from NRL plasma formulary
lnLambda_ee = 23.5 - log(sqrt(n_cgs)*T_eV^-1.25) - sqrt(1e-5 + (log(T_eV)-2)^2/16)
@test_approx_eq(lnLambda(electron,electron),lnLambda_ee)
lnLambda_ei = 24.0-log(sqrt(n_cgs)/T_eV)
@test_approx_eq(lnLambda(electron,helium3),lnLambda_ei)
@test_approx_eq(lnLambda(helium3,electron),lnLambda_ei)
lnLambda_ii = 23.0-log( (2/T_eV)*sqrt( (0.01*n_cgs*4/T_eV) + (0.5*n_cgs/T_eV)) )
@test_approx_eq(lnLambda(helium3,deuterium),lnLambda_ii)
@test_approx_eq(lnLambda(deuterium,helium3),lnLambda_ii)

# Test value of nuhat against deuterium
@test_approx_eq(nuhat(helium3,deuterium),10.411163139347321*lnLambda_ii)

# Test low-v limit of Chandrasekhar function
x = 1.0e-4
@test_approx_eq_eps(gfunc(x),2*x/(3.0*sqrt(pi)),1.0e-12)

# Test high-v limit of Chandrasekhar function
x = 8.0
@test_approx_eq(gfunc(x),0.5*x^-2)

# Test value of nu_s

# Test value of nu_par

# Test slowing-down time

# Test solver for simple problem

# Test that collision operator vanishes for Maxwellian
collop = collisionOperator(vgrid,helium3,bulkspecs)
testfunc = collop*f0trace
@test_approx_eq_eps(testfunc[1:Nv],zeros(Nv),1.0e-3)

# Test that solving collision operator gives Maxwellian

println("All tests passed.")
