using DelimitedFiles

using CSV
using PyPlot

rc("font", family="Arial", size=16.0)
#pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)
ns = 20
tf = 0.01
gamma = 1.4

solution = readdlm("solution_d.txt")#, datarow = 3, type=Float64)
xd =  solution[:,1]
ud = solution[:,2:ns+2]

solution = readdlm("solution_v.txt")#, datarow = 3, type=Float64)
xv =  solution[:,1]
uv = solution[:,2:ns+2]

solution = readdlm("solution_e.txt")#, datarow = 3, type=Float64)
xe =  solution[:,1]
ue = solution[:,2:ns+2]

nx = length(xv)
xp = xe
up = Array{Float64}(undef, nx,ns+1)

for i = 1:nx for j = 1:ns+1
        up[i,j] = (gamma-1.0)*(ue[i,j] - uv[i,j]*uv[i,j]/ud[i,j])
        uv[i,j] = uv[i,j]/ud[i,j]
        ue[i,j] = ue[i,j]/ud[i,j]
end end

solution = readdlm("solution_dF.txt")#, datarow = 3, type=Float64)
xdF =  solution[:,1]
udF = solution[:,2:ns+2]

solution = readdlm("solution_vF.txt")#, datarow = 3, type=Float64)
xvF =  solution[:,1]
uvF = solution[:,2:ns+2]

solution = readdlm("solution_eF.txt")#, datarow = 3, type=Float64)
xeF =  solution[:,1]
ueF = solution[:,2:ns+2]

nxF = length(xvF)
xpF = xeF
upF = Array{Float64}(undef, nxF,ns+1)

for i = 1:nxF for j = 1:ns+1
        upF[i,j] = (gamma-1.0)*(ueF[i,j] - uvF[i,j]*uvF[i,j]/udF[i,j])
        uvF[i,j] = uvF[i,j]/udF[i,j]
        ueF[i,j] = ueF[i,j]/udF[i,j]
end end

fig = figure("An example", figsize=(16,14));
ax1 = fig[:add_subplot](2,2,1);
ax2 = fig[:add_subplot](2,2,2);
ax3 = fig[:add_subplot](2,2,3);
ax4 = fig[:add_subplot](2,2,4);

ax1.plot(xdF, udF[:,ns+1], lw=4, color="black", label="True")
ax1.plot(xd, ud[:,ns+1], lw=2, color="blue", marker = "o",
        markersize=4,label="Low-resolution")
ax1.set_xlabel("\$x\$")
ax1.set_ylabel("\$ρ\$")
ax1.set_title("Density")
ax1.set_xlim(0,1)
ax1.legend(fontsize=14, loc=1)

ax2.plot(xvF, uvF[:,ns+1], lw=4, color="black", label="True")
ax2.plot(xv, uv[:,ns+1], lw=2, color="red", marker = "o",
        markersize=4,label="Low-resolution")
ax2.set_xlabel("\$x\$")
ax2.set_ylabel("\$v\$")
ax2.set_title("Velocity")
ax2.set_xlim(0,1)
ax2.legend(fontsize=14, loc=0)

ax3.plot(xeF, ueF[:,ns+1], lw=4, color="black", label="True")
ax3.plot(xe, ue[:,ns+1], lw=2, color="green", marker = "o",
        markersize=4,label="Low-resolution")
ax3.set_xlabel("\$x\$")
ax3.set_ylabel("\$e\$")
ax3.set_title("Energy")
ax3.set_xlim(0,1)
ax3.legend(fontsize=14, loc=0)

ax4.plot(xpF, upF[:,ns+1], lw=4, color="black", label="True")
ax4.plot(xp, up[:,ns+1], lw=2, color="orange", marker = "o",
        markersize=4,label="Low-resolution")
ax4.set_xlabel("\$x\$")
ax4.set_ylabel("\$p\$")
ax4.set_title("Pressure")
ax4.set_xlim(0,1)
ax4.legend(fontsize=14, loc=0)

plt[:subplot](ax1);
plt[:subplot](ax2);
plt[:subplot](ax3);
plt[:subplot](ax4);

fig.tight_layout()
fig.savefig("euler_1D.pdf")
