using NQCModels
include("HGe_model/HGemodel.jl")
function unit_atomicToclassical(model::NQCModels.Morse; dict::Dict{Symbol, Any} = Dict{Symbol, Any}())
    """
    Morse potential: U_0(x) [:Dₑ, :x₀, :a]
    """
    morse_parameters = [:Dₑ, :x₀, :a]
    for key in morse_parameters
        try
            value = getfield(model, key)
        catch e
            error("Input model::Morse does not contain the field $key")
        end
        value = getfield(model, key)
        ## Convert to classical units and add it to the dictionary
        if key == :Dₑ
            dict[key] = auconvert(u"eV", value)
        elseif key == :x₀
            dict[key] = auconvert(u"Å", value)
        elseif key == :a
            dict[key] = auconvert(u"Å^-1", value)
        end
    end
    
    global add_prime = true
    return dict

end

function unit_atomicToclassical(model::Logistic; dict::Dict{Symbol, Any} = Dict{Symbol, Any}(), add_prime::Bool = false)
    """
    Logistic potential: h(x) [:L, :k, :x₀, :c, :a]

    If this function is called after unit_atomicToclassical(model::Morse), 
        some parameters ([:x₀′, :c′, :a′]) of the model need to add '
    """
    logistic_parameters = [:L, :k, :x₀, :c, :a]
    if add_prime
        dict_keys = [:L, :k, :x₀′, :c′, :a′]
    else
        dict_keys = logistic_parameters
    end

    for (i,key) in enumerate(logistic_parameters)
        try
            value = getfield(model, key)
        catch e
            error("Input model::Logistic does not contain the field $key")
        end
        value = getfield(model, key)
        ## Convert to classical units and add it to the dictionary
        if key == :L
            dict[dict_keys[i]] = auconvert(u"eV", value)
        elseif key == :k
            dict[dict_keys[i]] = auconvert(u"Å^-1", value)
        elseif key == :x₀
            dict[dict_keys[i]] = auconvert(u"Å", value)
        elseif key == :c
            dict[dict_keys[i]] = auconvert(u"eV", value)
        elseif key == :a
            dict[dict_keys[i]] = value
        end
    end
    return dict

end


function unit_atomicToclassical(model::Hokseon)
    """
    Hokseon model contains three parts:

        Morse potential (model.Morse): U_0(x) [:Dₑ, :x₀, :a] and its vertical shirft [:c]

        Logistic potential (model.logistic): h(x) [:L, :k, :x₀, :c, :a]

        Hyperbolic tangent coupling (model.parameters): A(x) [:q, :ã, :x̃, :Ā]
    
    """
    dict = Dict{Symbol, Any}(
        :morse => Dict{Symbol, Any}(),
        :logistic => Dict{Symbol, Any}(),
        :tanh => Dict{Symbol, Any}()
    )
    ## Morse potential
    try
        morse = model.morse
    catch e
        error("Input model does not contain the field model.morse")
    end
    morse = model.morse
    dict[:morse] = unit_atomicToclassical(morse)
    # Vertical shift of the Morse potential
    try
        c = model.c
    catch e
        error("Input model does not contain the field model.c")
    end
    c = model.c
    dict[:morse][:c] = auconvert(u"eV", c)


    ## Logistic potential
    try
        logistic = model.logistic
    catch e
        error("Input model does not contain the field model.logistic")
    end
    logistic = model.logistic
    dict[:logistic] = unit_atomicToclassical(logistic; add_prime)


    ## Hyperbolic tangent coupling
    tanh_parameters = [:q, :ã, :x̃, :Ā]
    for key in tanh_parameters
        try
            value = getfield(model, key)
        catch e
            error("Input model does not contain the field model.parameters.$key")
        end
        value = getfield(model, key)
        ## Convert to classical units and add it to the dictionary
        if key == :q
            dict[:tanh][key] = value
        elseif key == :ã
            dict[:tanh][key] = auconvert(u"Å", value)
        elseif key == :x̃
            dict[:tanh][key] = auconvert(u"Å", value)
        elseif key == :Ā
            dict[:tanh][key] = auconvert(u"eV", value)
        end
    end

    return dict
end