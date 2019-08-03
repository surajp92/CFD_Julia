using DelimitedFiles

using CSV
using PyPlot

rc("font", family="Arial", size=16.0)
#pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)
ns = 10
tf = 0.25

solution = readdlm("solution.txt")#, datarow = 3, type=Float64)
x =  solution[:,1]
u = solution[:,2:ns+1]

fig2 = figure("An example", figsize=(8,6));
ax3 = fig2[:add_subplot](1,1,1);

for i = 1:ns-1
    ax3.plot(x, u[:,i], lw=1, label="t = $(tf*i/ns)")
end
ax3.set_xlabel("\$x\$")
ax3.set_ylabel("\$u\$")
ax3.set_xlim(0,1)
ax3.legend(fontsize=14, loc=0, bbox_to_anchor=(0.3, 0.45, 0.5, 0.5))

fig2.tight_layout()
fig2.savefig("burgers_cds.pdf")
