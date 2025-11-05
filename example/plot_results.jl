using JSON
using Plots
using Statistics

# Read the results file
data = JSON.parsefile("example_job.results.json", allownan=true)

# Dictionary to store data by system size
sizes = Dict{Tuple{Float64,Float64}, Dict{String,Vector{Float64}}}()

# Parse data grouped by system size
for entry in data
    T = entry["parameters"]["T"]
    Lx = entry["parameters"]["Lx"]
    Ly = entry["parameters"]["Ly"]
    size_key = (Lx, Ly)
    
    # Initialize dictionary for this size if needed
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
    
    push!(s["T"], T)
    
    # Energy per site (scaled by N)
    push!(s["E"], entry["results"]["Energy"]["mean"])
    push!(s["E_err"], entry["results"]["Energy"]["error"] / N)
    
    # Absolute magnetization per site (scaled by N)
    push!(s["M"], entry["results"]["AbsMagnetization"]["mean"])
    push!(s["M_err"], entry["results"]["AbsMagnetization"]["error"] / N)
    
    # Specific heat per site (scaled by N)
    push!(s["Cv"], entry["results"]["SpecificHeat"]["mean"][1])
    push!(s["Cv_err"], entry["results"]["SpecificHeat"]["error"][1] / N)
    
    # Susceptibility per site (scaled by N)
    push!(s["Chi"], entry["results"]["Susceptibility"]["mean"][1])
    push!(s["Chi_err"], entry["results"]["Susceptibility"]["error"][1] / N)
    
    # Binder ratio (NOT scaled - dimensionless)
    push!(s["U"], entry["results"]["BinderRatio"]["mean"][1])
    push!(s["U_err"], entry["results"]["BinderRatio"]["error"][1])
end

# Sort each size by temperature
for (size_key, s) in sizes
    idx = sortperm(s["T"])
    for key in keys(s)
        s[key] = s[key][idx]
    end
end

# Define colors and markers for each size
colors = [:blue, :red, :green]
markers = [:circle, :square, :diamond]
size_keys = sort(collect(keys(sizes)))

# Create plots
# p1 = plot(xlabel="Temperature", ylabel="Energy per site",
#          title="Energy vs Temperature", legend=:bottomright)
# p2 = plot(xlabel="Temperature", ylabel="|Magnetization| per site",
#          title="Magnetization vs Temperature", legend=:topright)
# p3 = plot(xlabel="Temperature", ylabel="Specific Heat per site",
#          title="Specific Heat vs Temperature", legend=:topright)
# p4 = plot(xlabel="Temperature", ylabel="Susceptibility per site",
#          title="Susceptibility vs Temperature", legend=:topright)
p5 = plot(xlabel="Temperature", ylabel="Binder Ratio",
         title="Binder Ratio vs Temperature", legend=:bottomleft)

for (i, size_key) in enumerate(size_keys)
    s = sizes[size_key]
    Lx, Ly = size_key
    N = Lx*Ly
    label = "L=$(Int(Lx))×$(Int(Ly))"
    
#     plot!(p1, s["T"], s["E"]./N*N, yerror=s["E_err"],
#           label=label, marker=markers[i], color=colors[i], markersize=4,
#           markerstrokewidth=0)
    
#     plot!(p2, s["T"], s["M"]./N*N, yerror=s["M_err"],
#           label=label, marker=markers[i], color=colors[i], markersize=4,
#           markerstrokewidth=0)
    
    # plot!(p3, s["T"], s["Cv"]./N*N, yerror=s["Cv_err"],
    #       label=label, marker=markers[i], color=colors[i], markersize=4,
    #       markerstrokewidth=0)
    
#     # plot!(p4, s["T"], s["Chi"]./N*N, yerror=s["Chi_err"],
#     #       label=label, marker=markers[i], color=colors[i], markersize=4,
#     #       markerstrokewidth=0)
    
    plot!(p5, s["T"], s["U"], yerror=s["U_err"],
          label=label, marker=markers[i], color=colors[i], markersize=4,
          markerstrokewidth=0)
end

# Combine into a single figure
# plot(p1, p2, p3, p5, layout=(2, 2), size=(1200, 1200))
# plot(p3, size=(1200, 1200))
savefig("20x20Carlo.png")

println("Plot saved to ising_results.png")
