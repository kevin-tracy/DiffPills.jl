__precompile__(true)
module DifferentialProximity

using LinearAlgebra
using StaticArrays
import ForwardDiff as FD
# import FiniteDiff
using Printf

include("utils/mrp.jl")
include("utils/quaternion.jl")
include("qp_solvers/interior_point.jl")
include("qp_solvers/active_set.jl")
include("capsules/proximity.jl")
include("capsules/proximity_jacobians.jl")
include("polygons/proximity.jl")
include("polygons/proximity_jacobians.jl")
include("capsule_v_polygon/proximity.jl")

end # module
