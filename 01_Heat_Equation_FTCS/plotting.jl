clearconsole()

using CSV
using Plots
font = Plots.font("Times New Roman", 18)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

final_field = CSV.read("field_final.csv")#, datarow = 2, type=Float64)

x = convert(Array,final_field[:,1])

u_e = convert(Array,final_field[:,2])
u_n = convert(Array,final_field[:,3])
u_error = convert(Array,final_field[:,4])

u = Array{Float64}(undef, length(u_e), 2)
u[:,1] = u_e
u[:,2] = u_n

for i = 1:Int64(length(u_error))
    u_error[i] = abs(u_error[i])
end

# plotting the exact and numerical solution
styles = [:solid, :dot]
styles = reshape(styles, 1, length(styles))
color=[:red :blue]

p1 = plot(x,u, line = (6,styles), color=color,
               xlabel="\$X\$", xlims=(minimum(x),maximum(x)),
               ylabel = "\$U\$",
               grid=(:none), label=["Exact solution" "FTCS scheme"],
               legend=:topright)

plot(p1, size = (800, 400))
savefig("ftcs.pdf")

p2 = plot(x,u_error, line = (3,styles), linestyle=:dot, color=:green,
                     xlabel="\$X\$", xlims=(minimum(x),maximum(x)),
                     ylabel="\$ϵ\$", grid=(:none),
                     label="Error", legend=:topright,
                     markershape = :circle, markercolor = :green,
                     markerstrokecolor = :black, markersize = 7)

plot(p2, size = (400, 400))
savefig("ftcs_error.pdf")

plot(p1,p2, size=(1400,500))
savefig("ftcs_all.pdf")
