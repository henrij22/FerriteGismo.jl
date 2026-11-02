module FerriteGismo

# ==============================================================================
# DEPENDENCIES
# ==============================================================================

# Core Julia packages
using LinearAlgebra

# Utility packages
using Reexport
using ArgCheck
import WriteVTK

# Backend not exported
using TinyGismo: gsBasis, gsGeometry

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
include("dofs/constraints.jl")

# Iterators
include("iterators.jl")

# FEValues
include("fevalues/cellvalues.jl")
include("fevalues/facetvalues.jl")

# Export
include("export/exportmesh.jl")
include("export/l2projection.jl")
include("export/vtk.jl")
include("export/bezier.jl")

# ==============================================================================
# EXPORTS
# ==============================================================================

# Grid
export IGAGrid, numElements, parameterSpaceGrid, numElementsPerDirection
export hierarchicalSubdomains

# Export mesh / postprocessing
export exportGrid, exportPoints, evaluateAtExportNodes, breakpoints, isHierarchical
export VTKBezierFile

# Elements

# Dofs
export IGADofHandler, fixEdge!, prescribeEdge!, edgeControlPoints, fieldOffset

# Interpolations
export IGAInterpolation

# Utility
export interpolate

end # module FerriteGismo
