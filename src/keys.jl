"""
    keys_template(dict::Dict;thofkey::Int=1)

    This function is to extract the keys template from the dictionary.

        dict: Dict with the simulation results
        thofkey: Int, the index of the key to extract the template
        return: String, the keys template
    
    example
    ```julia
    dict = Dict("gap=0.1_incident_energy=0.1"=>[1,2,3], "gap=0.1_incident_energy=0.2"=>[1,2,3])
    keys_template(dict) # return "gap=0.1_incident_energy="
    ```
"""
function keys_template(dict::Dict; thofkey::Int=1)
    first_key = collect(keys(dict))[thofkey]
    parts = split(first_key, "=")
    keys_template = join(parts[1:2], "=") * "="
    return keys_template
end