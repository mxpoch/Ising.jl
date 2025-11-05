using JSON
using Plots
using Statistics

# Use the name of the attached file
data = JSON.parsefile("ising_pt_30.results.json", allownan=true)

# Dictionary to store data grouped by system size (Lx, Ly)
sizes = Dict{Tuple{Float64,Float64}, Dict{String,Vector{Float64}}}()

# Parse data from each entry in the JSON array
for entry in data
    # Extract system dimensions from parameters
    Lx = entry["parameters"]["Lx"]
    Ly = entry["parameters"]["Ly"]
    size_key = (Lx, Ly)
    
    # Initialize a dictionary for this system size if it's the first time we've seen it
    if !haskey(sizes, size_key)
        sizes[size_key] = Dict(
            "T" => Float64[],
            "E" => Float64[], "E_err" => Float64[],
            "M" => Float64[], "M_err" => Float64[],
            "Cv" => Float64[], "Cv_err" => Float64[],
            "Chi" => Float64[], "Chi_err" => Float64[],
            "U" => Float64[], "U_err" => Float64[]
        )
    end
    
    s = sizes[size_key]
    N = Lx * Ly
    
    # Extract the array of temperatures for this simulation
    temperatures = entry["parameters"]["parallel_tempering"]["values"]
    
    # Loop over each temperature and its corresponding result
    for (i, T) in enumerate(temperatures)
        push!(s["T"], T)
        
        # Energy per site
        push!(s["E"], entry["results"]["Energy"]["mean"][i])
        push!(s["E_err"], entry["results"]["Energy"]["error"][i] / N)
        
        # Absolute magnetization per site
        push!(s["M"], entry["results"]["AbsMagnetization"]["mean"][i])
        push!(s["M_err"], entry["results"]["AbsMagnetization"]["error"][i] / N)
        
        # Specific heat per site
        push!(s["Cv"], entry["results"]["SpecificHeat"]["mean"][i])
        push!(s["Cv_err"], entry["results"]["SpecificHeat"]["error"][i] / N)
        
        # Susceptibility per site
        push!(s["Chi"], entry["results"]["Susceptibility"]["mean"][i])
        push!(s["Chi_err"], entry["results"]["Susceptibility"]["error"][i] / N)
        
        # Binder ratio (dimensionless, not scaled by N)
        push!(s["U"], entry["results"]["BinderRatio"]["mean"][i])
        push!(s["U_err"], entry["results"]["BinderRatio"]["error"][i])
    end
end

# Sort the data for each system size by temperature for correct plotting
for (size_key, s) in sizes
    p = sortperm(s["T"])
    for key in keys(s)
        s[key] = s[key][p]
    end
end

# --- Plotting ---

# Define colors and markers for each system size to distinguish them on the plots
colors = [:blue, :red, :green, :purple, :orange]
markers = [:circle, :square, :diamond, :utriangle, :dtriangle]
size_keys = sort(collect(keys(sizes)))

# Create individual plots for each physical quantity
p1 = plot(xlabel="Temperature", ylabel="Energy per site",
         title="Energy vs Temperature", legend=:bottomright)
p2 = plot(xlabel="Temperature", ylabel="|Magnetization| per site",
         title="Magnetization vs Temperature", legend=:topright)
p3 = plot(xlabel="Temperature", ylabel="Specific Heat per site",
         title="Specific Heat vs Temperature", legend=:topright)
# The susceptibility plot is commented out as in the original script
# p4 = plot(xlabel="Temperature", ylabel="Susceptibility per site",
#          title="Susceptibility vs Temperature", legend=:topright)
p5 = plot(xlabel="Temperature", ylabel="Binder Ratio",
         title="Binder Ratio vs Temperature", legend=:bottomleft)

# Add data series to each plot for each system size
for (i, size_key) in enumerate(size_keys)
    s = sizes[size_key]
    Lx, Ly = size_key
    N = Lx * Ly
    label = "L=$(Int(Lx))×$(Int(Ly))"
    
    # Plot Energy (mean is scaled to per-site, error is already per-site)
    plot!(p1, s["T"], s["E"], yerror=s["E_err"],
          label=label, marker=markers[i], color=colors[i], markersize=4,
          markerstrokewidth=0)
    
    # Plot Magnetization (mean is scaled to per-site, error is already per-site)
    plot!(p2, s["T"], s["M"], yerror=s["M_err"],
          label=label, marker=markers[i], color=colors[i], markersize=4,
          markerstrokewidth=0)
    
    # Plot Specific Heat (mean is scaled to per-site, error is already per-site)
    plot!(p3, s["T"], s["Cv"], yerror=s["Cv_err"],
          label=label, marker=markers[i], color=colors[i], markersize=4,
          markerstrokewidth=0)
    
    # Plot Susceptibility (if uncommented)
    # plot!(p4, s["T"], s["Chi"] ./ N, yerror=s["Chi_err"],
    #       label=label, marker=markers[i], color=colors[i], markersize=4,
    #       markerstrokewidth=0)
    
    # Plot Binder Ratio (not scaled by N)
    plot!(p5, s["T"], s["U"], yerror=s["U_err"],
          label=label, marker=markers[i], color=colors[i], markersize=4,
          markerstrokewidth=0)
end

# Combine all plots into a single figure with a 2x2 layout
plot(p1, p2, p3, p5, layout=(2, 2), size=(1200, 1000))

# Save the final figure to a file
savefig("ising_results.png")

println("Plot saved to ising_results.png")
