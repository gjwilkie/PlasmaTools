module grids

export Grid,createFVgrid, createAdvecDiffOperator

type Grid
   lft::Array{Float64,1}
   ctr::Array{Float64,1}
   rgt::Array{Float64,1}
   BC::Array{Int64,1}
end

function createDifferentiationMatrix(vgrid::Array{Float64,1})
   Nv = length(vgrid)
   for i in 1:Nv
     ddv[i,2:end-1] = 1.0/(vgrid[3:end]-vgrid[1:end-2])
     ddv[i,1] = 1.0/(vgrid[2]-vgrid[1])
     ddv[i,end] = 1.0/(vgrid[end]-vgrid[end-1])
   end
end

function createFVgrid(Nv::Int64,v0::Float64,vf::Float64,BCs::Array{Int64,1})
   leftFaces = Float64[]
   centerPts = Float64[]
   rightFaces = Float64[]
   dv = (vf-v0)/Nv
   for i in 1:Nv
      push!(leftFaces, (i-1)*dv)
   end
   rightFaces = leftFaces[2:end]
   push!(rightFaces,leftFaces[end] + (leftFaces[end]-leftFaces[end-1]))
   
   centerPts = 0.5*( leftFaces[1:end] + rightFaces[1:end] )

   return Grid(leftFaces,centerPts,rightFaces,BCs)
end

function createAdvecDiffOperator(grid::Grid,jacobian::Function,advection::Function,diffusion::Function,boundaryVals::Array{Float64,1},source::Array{Float64,1})

   N = length(grid.ctr)
   a = zeros(Float64,N)
   b = zeros(Float64,N)
   c = zeros(Float64,N)
   operator = spzeros(N,N)

   # See ../docs/operator.pdf for description of this 
   # Interior points
   denom = jacobian(grid.ctr[2:end-1]).*(grid.rgt[2:end-1]-grid.lft[2:end-1])
   A_left = jacobian(grid.lft[2:end-1]).*advection(grid.lft[2:end-1])
   A_right = jacobian(grid.rgt[2:end-1]).*advection(grid.rgt[2:end-1])
   D_left = (jacobian(grid.lft[2:end-1]).*diffusion(grid.lft[2:end-1]))./(grid.ctr[2:end-1] - grid.ctr[1:end-2])
   D_right = (jacobian(grid.rgt[2:end-1]).*diffusion(grid.rgt[2:end-1]))./(grid.ctr[3:end] - grid.ctr[2:end-1])

   a[1:end-2] = (-0.5*A_left + D_left)./denom
   b[2:end-1] = (0.5*A_right - 0.5*A_left - D_right - D_left)./denom
   c[3:end] = (0.5*A_right + D_right)./denom

   # i=1 boundary condition
   if grid.BC[1] == 0
      # Dirichlet boundary conditions
      scaleFactor = b[2]
      b[1] = scaleFactor
      source[1] = boundaryVals[1]*scaleFactor
   elseif grid.BC[1] == 1
      # Flux-specified boundary
      denom = jacobian(grid.ctr[1])*(grid.rgt[1]-grid.lft[1])
      b[1] = (0.5*jacobian(grid.rgt[1]).*advection(grid.rgt[1]) - (jacobian(grid.rgt[1]).*diffusion(grid.rgt[1]))./(grid.ctr[2] - grid.ctr[1]) )/denom
      c[2] = (0.5*jacobian(grid.rgt[1]).*advection(grid.rgt[1]) + (jacobian(grid.rgt[1]).*diffusion(grid.rgt[1]))./(grid.ctr[2] - grid.ctr[1]) )/denom
      source[1] = source[1] + jacobian(grid.lft[1])*boundaryVals[1] / (jacobian(grid.ctr[1])*(grid.rgt[1]-grid.lft[1]))
   else
         error("Invalid boundary condition for i = 1")
   end

   # i=N boundary condition
   if grid.BC[2] == 0
      # Dirichlet boundary conditions
      scaleFactor = b[N-1]
      b[N] = scaleFactor
      source[N] = boundaryVals[2]*scaleFactor
   elseif grid.BC[2] == 1
      # Flux-specified boundary
      denom = jacobian(grid.ctr[N])*(grid.rgt[N]-grid.lft[N])
      a[N-1] = (-0.5*jacobian(grid.lft[N]).*advection(grid.lft[N]) + (jacobian(grid.lft[N]).*diffusion(grid.lft[N]))./(grid.ctr[N] - grid.ctr[N-1]) )/denom
      b[N] = (- 0.5*jacobian(grid.lft[N]).*advection(grid.lft[N]) - (jacobian(grid.lft[1]).*diffusion(grid.lft[1]))./(grid.ctr[N] - grid.ctr[N-1]) )/denom
      source[N] = source[N] - jacobian(grid.rgt[N])*boundaryVals[2] / (jacobian(grid.ctr[N])*(grid.rgt[N]-grid.lft[N]))
   else
         error("Invalid boundary condition for i = N")
   end

   for icol in 2:N-1
      operator[icol-1:icol+1,icol] = [a[icol],b[icol],c[icol]]
   end
   operator[1,1] = b[1]
   operator[1,2] = c[2]
   operator[N,N-1] = a[N-1]
   operator[N,N] = b[N]
   
   return operator

end

end
