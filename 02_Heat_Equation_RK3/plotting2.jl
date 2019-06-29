using CSV
using PyPlot
rc("font", family="Arial", size=16.0)


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

fig = figure("An example", figsize=(14,6));
ax1 = fig[:add_subplot](1,2,1);
ax2 = fig[:add_subplot](1,2,2);

ax1.plot(x, u_e, lw=4, ls = "-", color="b", label="Exact solution")
ax1.plot(x, u_n, lw=4, ls = "--", color="r", label="FTCS solution")
ax1.set_xlabel("\$x\$")
ax1.set_ylabel("\$v\$")
ax1.set_title("Solution field")
ax1.set_xlim(-1,1)
ax1.legend(fontsize=14, loc=0)

ax2.plot(x, u_error, marker = "o", markeredgecolor="k",
        markersize=8, color="g", lw=4)
ax2.set_xlabel("\$x\$")
ax2.set_ylabel("\$ϵ\$")
ax2.set_title("Discretization error")
ax2.set_xlim(-1,1)
#ax2.legend(fontsize=14, loc=0)

plt[:subplot](ax1);
plt[:subplot](ax2);

fig.tight_layout()
fig.savefig("rk3.pdf")
