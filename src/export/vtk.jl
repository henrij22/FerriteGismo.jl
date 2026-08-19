# ==============================================================================
# VTK export
#
# The spline patch itself cannot be written to a VTK file: the formats know linear (and a
# few higher-order Lagrange) cells, not knot spans. FerriteGismo therefore hands the
# writers the *export mesh* — an ordinary Ferrite grid of the knot spans, see
# `exportmesh.jl` — and evaluates the spline fields at its nodes.
#
# Everything below is built on the public Ferrite API: the `VTKGridFile` /
# `VTKHDFGridFile` constructors for a plain grid, and `write_node_data`. What is added
# here are methods of *public* Ferrite functions dispatched on FerriteGismo's own types.
# ==============================================================================

"""
    VTKGridFile(filename, dh::IGADofHandler; subdivision::Int = 1, kwargs...)
    VTKHDFGridFile(filename, dh::IGADofHandler; subdivision::Int = 1, kwargs...)

Open a VTK (or VTKHDF) file for the grid of `dh`, written as the export mesh of the patch,
see [`exportGrid`](@ref). `subdivision` samples each knot span that many times per
direction; the same value then has to be passed to `write_solution` and `write_projection`,
since the data has to match the nodes of the file.

`VTKHDFGridFile` additionally requires
[VTKHDF.jl](https://github.com/Ferrite-FEM/VTKHDF.jl) to be loaded, as for a regular grid.

```julia
VTKGridFile("solution", dh) do vtk
    write_solution(vtk, dh, u)
end
```
"""
function Ferrite.VTKGridFile(filename::String, dh::IGADofHandler; subdivision::Int = 1, kwargs...)
    return VTKGridFile(filename, exportGrid(Ferrite.get_grid(dh); subdivision); kwargs...)
end

function Ferrite.VTKHDFGridFile(filename::String, dh::IGADofHandler; subdivision::Int = 1, kwargs...)
    return Ferrite.VTKHDFGridFile(filename, exportGrid(Ferrite.get_grid(dh); subdivision); kwargs...)
end

# `Ferrite.VTKHDFGridFile` always exists (Ferrite declares the type, its functionality lives
# in Ferrite's VTKHDF extension), so both file types are covered without depending on VTKHDF.
const VTKFile = Union{VTKGridFile, Ferrite.VTKHDFGridFile}

"""
    write_solution(vtk, dh::IGADofHandler, u::AbstractVector, suffix = ""; subdivision::Int = 1)

Write every field of the solution vector `u` to `vtk`, evaluated at the nodes of the export
mesh. Pass the same `subdivision` that the file was opened with.
"""
# One method per file type rather than one for their union: Ferrite's own method is
# `(::VTKGridFile, ::AbstractDofHandler, ...)`, and a union in the first argument would tie
# with it instead of winning on the more specific dof handler.
function Ferrite.write_solution(
        vtk::VTKGridFile, dh::IGADofHandler, u::AbstractVector, suffix = ""; subdivision::Int = 1
    )
    return _writeSolution!(vtk, dh, u, suffix, subdivision)
end

function Ferrite.write_solution(
        vtk::Ferrite.VTKHDFGridFile, dh::IGADofHandler, u::AbstractVector, suffix = ""; subdivision::Int = 1
    )
    return _writeSolution!(vtk, dh, u, suffix, subdivision)
end

function _writeSolution!(vtk, dh::IGADofHandler, u::AbstractVector, suffix, subdivision::Int)
    for name in Ferrite.getfieldnames(dh)
        data = evaluateAtExportNodes(dh, u, name; subdivision)
        Ferrite.write_node_data(vtk, data, string(name, suffix))
    end
    return vtk
end

"""
    write_projection(vtk, proj::L2ProjectorIGA, vals, name; subdivision::Int = 1)

Write projected quadrature-point data `vals` (the output of `project`) to `vtk` under
`name`, evaluated at the nodes of the export mesh. Pass the same `subdivision` that the
file was opened with.
"""
function Ferrite.write_projection(
        vtk::VTKFile, proj::L2ProjectorIGA, vals, name; subdivision::Int = 1
    )
    data = Ferrite.evaluate_at_grid_nodes(proj, vals; subdivision)
    Ferrite.write_node_data(vtk, data, name)
    return vtk
end

"""
    evaluate_at_grid_nodes(dh::IGADofHandler, u::AbstractVector, fieldname::Symbol)

Values of the field `fieldname` of `u` at the nodes of the export mesh, i.e. at the
knot-span corners. Use [`interpolate`](@ref) to evaluate at an arbitrary parametric point,
or [`evaluateAtExportNodes`](@ref) to sample a subdivided mesh.
"""
function Ferrite.evaluate_at_grid_nodes(dh::IGADofHandler, u::AbstractVector, fieldname::Symbol)
    return evaluateAtExportNodes(dh, u, fieldname)
end
