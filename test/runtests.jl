#!/usr/bin/env julia
using PlasmaTools
using Base.Test

# Testing gs2ls
g = GS2geometry(equilibrium_option="eik",rhoc=0.6,irho=2,qinp=2.0,shat=1.0,shift=-0.1,akappa=0.2,rmaj=3.0,r_geo=3.1)

r = GS2resolution(nperiod=3,ky=0.2)

p = GS2params(runname="gs2ls_test",beta=0.003)

specs = GS2species[]
push!(specs, GS2species(dens=0.5,temp=1.0,mass=1.0,z=1.0,tprim=3.0,fprim=1.0))
push!(specs, GS2species(dens=0.5,temp=1.0,mass=1.5,z=1.0,tprim=3.0,fprim=1.0))
push!(specs, GS2species(dens=1.0,temp=1.1,mass=0.0002,z=-1.0,tprim=2.5,fprim=1.0))

createInputFile("gs2ls_test.in",p,r,g,specs)
info("Wrote a gs2 input file: gs2ls_test.in. Check for correctness.")


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

# Test solver for simple problem, advection only
f0anal = exp(-vgrid.ctr.^2/vti^2)
A0 = 5.0
D0 = 1.0
onefunc(x) = 1.0
Afunc(x) = A0
Dfunc(x) = D0
source = f0anal.*( (D0/vti^2)*(4.0*(vgrid.ctr.^2/vti^2) - 2.0) - 2.0*vgrid.ctr*A0/vti^2)

operator = createAdvecDiffOperator(vgrid,onefunc,Afunc,Dfunc,[0.0,0.0],source)
f0test = operator\source
println(f0test)
println(f0anal)
@test_approx_eq_eps(f0test,f0anal,1.0e-3)


# Test solver for simple problem, diffusion only

# Test that collision operator vanishes for Maxwellian
collop = collisionOperator(vgrid,helium3,bulkspecs)
testfunc = collop*f0trace
@test_approx_eq_eps(testfunc[1:Nv],zeros(Nv),1.0e-3)

# Test that solving collision operator gives Maxwellian

println("All tests passed.")

