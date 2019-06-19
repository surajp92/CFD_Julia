clearconsole()
using DelimitedFiles

using CSV
using PyPlot
#using GR
#using Plots
#plotly()
# font = plt.font("Times New Roman", 18)
# pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

nx = 128
ny = 128
residual_hist = readdlm("residual.txt")#, datarow = 3, type=Float64)
iter_hist =  residual_hist[:,1]
res_hist = convert(Matrix, residual_hist[:,2:3])

color=[:red :blue]

fig, ax = plt.subplots()
ax.semilogy(iter_hist, res_hist[:,1], color="red",lw=2,label="rms")
ax.semilogy(iter_hist, res_hist[:,2], color="blue",lw=2,label="rms/rms_0")
ax.set_xlim(0,30000)
ax.legend()
fig.savefig("gs_residual.pdf")

init_field = readdlm("field_initial.txt")#, type=Float64)
final_field = readdlm("field_final.txt")#, datarow = 3, type=Float64)

x = convert(Array,init_field[:,1])
y = convert(Array,init_field[:,2])

u_e = convert(Array,init_field[:,5])
u_e = reshape(u_e, (nx+1, ny+1))

u_n = convert(Array,final_field[:,4])
u_n = reshape(u_n, (nx+1, ny+1))

xx = x[1:nx+1]
yy = reshape(y, (nx+1, ny+1))[1,:]

XX = repeat(xx, 1, length(yy))
XX = convert(Matrix,transpose(XX))
YY = repeat(yy, 1, length(xx))

fig, ax = plt.subplots(figsize=(7,5))
cs = ax.contourf(xx, yy, transpose(u_e),levels=20, cmap="jet")
ax.set_xlabel("X")
ax.set_ylabel("Y")
cbar = fig.colorbar(cs)
fig.tight_layout()
fig.savefig("gs_contour.pdf")
