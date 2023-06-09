"""
`Inlet1D` is a struct that contains the inlet boundary condition for a 1D flame.
    - `T`: inlet temperature (K)
"""
struct Inlet1D
    T::Float64
    cs::Array{Float64,1}
end