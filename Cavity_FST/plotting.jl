clearconsole()
using DelimitedFiles

using CSV
using PyPlot
#using GR
#using Plots
#plotly()
# font = plt.font("Times New Roman", 18)
# pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

nx = 64
ny = 64
residual_hist = readdlm("residual_plot.txt")#, datarow = 3, type=Float64)
iter_hist =  residual_hist[:,1]
res_hist = residual_hist[:,2]

fig, ax = plt.subplots()
ax.semilogy(iter_hist, res_hist[:,1], color="red",lw=2,label="rms")
ax.set_xlabel("Iteration count")
ax.set_ylabel("Residual \$L_2\$ norm")
ax.set_xlim(0,10000)
ax.legend()
fig.tight_layout()
fig.savefig("ldc_residual.pdf")

final_field = readdlm("field_final.txt")#, datarow = 3, type=Float64)

x = convert(Array,final_field[:,1])
y = convert(Array,final_field[:,2])

wn = convert(Array,final_field[:,3])
wn = reshape(wn, (nx+1, ny+1))

sn = convert(Array,final_field[:,4])
sn = reshape(sn, (nx+1, ny+1))

xx = x[1:nx+1]
yy = reshape(y, (nx+1, ny+1))[1,:]

XX = repeat(xx, 1, length(yy))
XX = convert(Matrix,transpose(XX))
YY = repeat(yy, 1, length(xx))

fig, ax = plt.subplots(figsize=(6,5))
cs = ax.contourf(xx, yy, transpose(wn),levels=60, cmap="jet")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Vorticity field")
cbar = fig.colorbar(cs)
fig.tight_layout()
fig.savefig("wn_contour.pdf")

fig, ax = plt.subplots(figsize=(6,5))
cs = ax.contourf(xx, yy, transpose(sn),levels=60, cmap="jet")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Streamfunction")
cbar = fig.colorbar(cs)
fig.tight_layout()
fig.savefig("sn_contour.pdf")
