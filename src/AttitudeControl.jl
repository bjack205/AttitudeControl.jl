module AttitudeControl

using LinearAlgebra

import Rotations
import ForwardDiff

include("dynamics.jl")
include("integration.jl")

end # module AttitudeControl
