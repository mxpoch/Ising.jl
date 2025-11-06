using Plots
using DataFrames
using Statistics
using Carlo.ResultTools

df = DataFrame(ResultTools.dataframe("diluted_PT.results.json"))
df.T = (pt -> pt["values"]).(df.parallel_tempering)
df.Lx = Int.(df.Lx)

filtered_df = df[df.DP .≈ 0.2, :]

# Group the DataFrame by the 'grouping' column
grouped_df = groupby(filtered_df, :Lx)

# Calculate the mean of 'Energy' and 'T' for each group
averaged_df = combine(grouped_df, :BinderRatio => mean, :T => mean)

colors = [:blue, :red, :green]
markers = [:circle, :square, :diamond]
plot(averaged_df.T_mean, averaged_df.BinderRatio_mean, ma=0.5, markersize=3, markers=reshape(markers, 1, length(markers)), palette=colors, title="Carlo.jl Ising Model 20.0% Dilution - Binder Parameter", xlabel = "Temperature", ylabel="Binder Ratio", group=averaged_df.Lx, legendtitle="L")
savefig("CarloPTResults_0.2.png")