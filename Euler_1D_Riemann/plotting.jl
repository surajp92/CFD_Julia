using DelimitedFiles

using CSV
using PyPlot

rc("font", family="Arial", size=16.0)
#pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)
ns = 20
tf = 0.01

solution = readdlm("solution_d.txt")#, datarow = 3, type=Float64)
xd =  solution[:,1]
ud = solution[:,2:ns+2]

solution = readdlm("solution_v.txt")#, datarow = 3, type=Float64)
xp =  solution[:,1]
up = solution[:,2:ns+2]

fig = figure("An example", figsize=(14,6));
ax1 = fig[:add_subplot](1,2,1);
ax2 = fig[:add_subplot](1,2,2);


ax1.plot(xd, ud[:,ns+1], lw=4, color="blue", label="t = 0.2")

ax1.set_xlabel("\$x\$")
ax1.set_ylabel("\$v\$")
ax1.set_title("Density")
ax1.set_xlim(0,1)
ax1.legend(fontsize=14, loc=1)


ax2.plot(xp, up[:,ns+1], lw=4, color="red", label="t = 0.2")

ax2.set_xlabel("\$x\$")
ax2.set_ylabel("\$v\$")
ax2.set_title("Velocity")
ax2.set_xlim(0,1)
ax2.legend(fontsize=14, loc=1)

plt[:subplot](ax1);
plt[:subplot](ax2);

fig.tight_layout()
fig.savefig("euler_1D.pdf")
