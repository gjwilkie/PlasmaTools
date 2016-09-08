module species

# TODO: Create custom constructor that creates v_t as part of type

export SpeciesData

"""
A container for species data. The following example creates a species of electrons in SI units.

 * `SpeicesData(mass=me,charge=-el,density=1.0e20,temp=1e3*el)`

SI units are not necessary, but is expected when passing to other methods in PlasmaTools.

"""
type SpeciesData
   mass::Float64
   charge::Float64
   density::Float64
   temperature::Float64
end

# Define constructor so that all keywords must be specified instead of ordered
SpeciesData(;mass=-1.0,charge=-1.0,density=-1.0,temperature=-1.0) =SpeciesData(mass,charge,density,temperature)

end
