"""
This scripts have some functions to calculate the number of bins for a histogram plot.
"""
function FreedmanDiaconisEstimator(data::Array{Float64,1})
    n = length(data)
    q75, q25 = quantile(data, [0.75, 0.25])
    iqr = q75 - q25
    h = 2 * iqr / n^(1/3)
    nbins = ceil((maximum(data) - minimum(data)) / h)
    return nbins
end

function SturgesEstimator(data::Array{Float64,1})
    n = length(data)
    nbins = ceil(log2(n) + 1)
    return nbins
end

function ScottEstimator(data::Array{Float64,1})
    n = length(data)
    h = 3.5 * std(data) / n^(1/3)
    nbins = ceil((maximum(data) - minimum(data)) / h)
    return nbins
end

function RiceEstimator(data::Array{Float64,1})
    n = length(data)
    h = 2 * n^(-1/3) * std(data)
    nbins = ceil((maximum(data) - minimum(data)) / h)
    return nbins
end