"""
Tools for model edge analysis for automatic mechanism generation
"""

using Logging
using Sundials
using SparseArrays
using DiffEqBase: build_solution

abstract type AbstractTerminationCriterion end

struct TerminationTime <: AbstractTerminationCriterion
    time::Float64
end
struct TerminationConversion <: AbstractTerminationCriterion
    species::Species
    conversion::Float64
end
struct TerminationRateRatio <: AbstractTerminationCriterion
    ratio::Float64
end

export TerminationTime
export TerminationConversion
export TerminationRateRatio