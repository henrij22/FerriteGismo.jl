module FerriteGismo

# ==============================================================================
# DEPENDENCIES
# ==============================================================================

# Core Julia packages
using LinearAlgebra

# Utility packages
using Reexport
using ArgCheck

# Backend not exported
using TinyGismo: gsBasis, gsGeometry
# using Gismo: numElements, copyMatrix, Basis, Geometry, getElements, actives, compute, boundary

# Backend packages (re-exported)
@reexport using Ferrite
@reexport using TinyGismo

# ==============================================================================
# MODULE STRUCTURE
# ==============================================================================

include("utility.jl")

# Interpolations
include("interpolations.jl")

# Grid
include("grid/grid.jl")

# Dofs
include("dofs/dofhandler.jl")

# Iterators
include("iterators.jl")

# FEValues
include("fevalues/cellvalues.jl")
include("fevalues/geometry_mapping.jl")

# Export
# include("export/vtk.jl")

# ==============================================================================
# EXPORTS
# ==============================================================================

# Grid
export IGAGrid, numElements, parameterSpaceGrid, numElementsPerDirection

# Elements

# Dofs
export IGADofHandler, fixBoundary!

# Interpolations
export IGAInterpolation

# FE Values
export IGACellValues

# Utility
export interpolate

end # module FerriteGismo
