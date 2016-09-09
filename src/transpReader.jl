#!/bin/julia
using NetCDF
using PyPlot

"""
plotSummary(file::ASCIIString)

Displays a set of time-dependent plots of relevant plasma quantities. Useful for determining the overall health of the run and ideal time samples.
"""
function plotSummary(file::ASCIIString)

   t = ncread(file,"TIME")
   Ti0 = ncread(file,"TI0")
   Te0 = ncread(file,"TE0")
   q95 = ncread(file,"Q95")
   ne0 = ncread(file,"NE")[1,:]'

   #TODO: Add some kind of global diagnostic of fast ions

end

"""
plotSummary(file::ASCIIString,time::Float64)

Displays a set of radius-dependent plots of relevant plasma quantities at the given time.
"""
function plotSummary(file::ASCIIString,time::Float64)

   # TODO: Make radial plots
   kappa = ncread(file,"ELONG")[1,:]
   tri = ncread(file,"TRIANG")[1,:]
   
end
