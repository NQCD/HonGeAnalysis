using Statistics
using Distributions  # Import the Distributions package
"""
data_confident_interval:
This function computes the confidence interval of the data.

    data: Array{Float64,1}, the data to compute the confidence interval
    confidence_level: Float64, the confidence level of the interval
"""

function data_confident_interval(data::Array{Float64,1}; confidence_level::Float64=0.95)
    n = length(data)
    mean_data = mean(data)
    std_data = std(data)
    z = quantile(Normal(), (1 + confidence_level) / 2)
    margin_error = z * std_data / sqrt(n)
    return mean_data, margin_error
end

function binomial_standard_error(n::Int, p::Float64)
    return sqrt(p * (1 - p) / n)
end

