
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

final_field = readdlm("field_final.txt")#, datarow = 3, type=Float64)

x = convert(Array,final_field[:,1])
y = convert(Array,final_field[:,2])

wn_t20 = convert(Array,final_field[:,3])
wn_t20 = reshape(wn_t20, (nx+1, ny+1))

field = readdlm("vm0.txt")
wn_t0 = convert(Array,field[:,3])
wn_t0 = reshape(wn_t0, (nx+1, ny+1))

field = readdlm("vm5.txt")
wn_t10 = convert(Array,field[:,3])
wn_t10 = reshape(wn_t10, (nx+1, ny+1))

xx = x[1:nx+1]
yy = reshape(y, (nx+1, ny+1))[1,:]

fig = figure("An example", figsize=(16,6));
ax1 = fig[:add_subplot](1,3,1);
ax2 = fig[:add_subplot](1,3,2);
ax3 = fig[:add_subplot](1,3,3);

cs = ax1.contourf(xx, yy, transpose(wn_t0),levels=80, cmap="YlGn",
                interpolation="bilinear")

ax1.set_title("\$ t=0\$")
plt[:subplot](ax1); im

cs = ax2.contourf(xx, yy, transpose(wn_t10),levels=80, cmap="YlGn",
                interpolation="bilinear")
ax2.set_title("\$ t=10\$")
plt[:subplot](ax2); cs

cs = ax3.contourf(xx, yy, transpose(wn_t20),levels=80, cmap="YlGn",
                 interpolation="bilinear")
ax3.set_title("\$ t=20\$")
plt[:subplot](ax3); cs


# cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])
# cb = plt.barh(ax1, cax = cbaxes)
fig.tight_layout()

fig.subplots_adjust(bottom=0.2)
cbaxes = fig.add_axes([0.2, 0.05, 0.6, 0.05])
fig.colorbar(cs, orientation="horizontal", cax = cbaxes)
# fig.colorbar(cs, ax = ax2), orientation='horizontal', pad=0.2
# fig.colorbar(cs, ax = ax3, orientation='horizontal', pad=0.2)

fig.savefig("vm_contour.eps")
