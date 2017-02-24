module CompScienceMeshes

using StaticArrays
import Base.getindex

export getcommonedge

export mesh, readmesh, writemesh
export dimension, universedimension, vertextype, celltype, coordtype
export meshsegment, meshrectangle, meshcircle, meshsphere, meshcuboid, meshwaveguidepost
export numvertices, vertices
export numcells, cells, cellvertices
export translate, translate!, rotate, rotate!, fliporientation!, fliporientation
export boundary, skeleton, vertextocellmap, connectivity, cellpairs
export barycentric_refinement
export bisecting_refinement

export simplex
export dimension, universedimension, vertextype, pointtype
export vertices, tangents, volume
export barytocart, carttobary, centroid
export euclidianbasis, defaultpointtype, point
export cartesian, jacobian, meshpoint
export intersection
export isinside
export isinclosure
export overlap

export meshpoint, meshpointtype
export meshpoints
export paramtype
export cartesian, parametric, barycentric
export jacobian, tangents, utangents, normal

export quadpoints
export trgauss, sqgauss, legendre

# TODO: remove this and use default mesh interface
export SegmentedAxis
export minmaxdist, rings, ring

typealias Pt{N,T} StaticArrays.SVector{N,T}

"""
    point(xs...)

Create point object of `CompScienceMeshes` default point class
and with coordinate type `Float64`
"""
point(xs...) = point(Float64, xs...)

"""
    point(type, xs...)

Create point object of `CompScienceMeshes` default point class
and with coordinate type `type`
"""
@generated function point(T::Type,xs...)
    D = length(xs)
    xp = :(SVector{$D,T}())
    for d in 1:D
        push!(xp.args, :(xs[$d]))
    end
    return xp
end

export index
"""
    index(i1, i2, ...) -> cell

Create a mesh cell by supplying the indices of the defining vertices in the
vertex buffer of the mesh.
"""
index(is...) = SVector{length(is),Int}(is...)


# """
#   defaultpointtype(T, dim) = SVector{T,dim}
#
# Returns the default point type used by package `CompScienceMeshes`
# """
# defaultpointtype(T, dim) = SVector{dim,T}



"""
  euclidian_basis(type, dim)

Returns an arrays of length (dim+1) containing the origin and the dim
Euclidian unit vectors with coordinate type `type`.
"""
function euclidianbasis(T::DataType, dim)
  P = SVector{dim,T}
  id = eye(dim)
  r = P[ P(id[:,i]...) for i in 1:dim ]
  z = P(zeros(T,dim)...)
  return unshift!(r, z)
end

include("combinatorics.jl")

# quadrature rules for segements, triangles, and squares
include("quadrature/SegmentGauss.jl")
include("quadrature/TriangleGauss.jl")
include("quadrature/SquareGauss.jl")

# simplices and related algorithms
include("patches.jl")
include("overlap.jl")
include("intersect.jl")
include("isinside.jl")
include("meshpoints.jl")
include("quadpoints.jl")

# mesh component
include("mesh.jl")
include("timeaxis.jl")

include("gmsh.jl")
include("gid.jl")
include("primitives.jl")
include("submesh.jl")
include("baryref.jl")
include("weld.jl")

# geometry API
include("geometry.jl")



include("utils.jl")

end # module
