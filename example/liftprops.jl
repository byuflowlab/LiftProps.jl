import XfoilSweep
import LiftProps
using Plots
pyplot()

airfoil_file = "naca2412.dat"

open(airfoil_file,"r") do f
  global x = []
  global y = []
  for line in eachline(f)
    x = append!(x,parse(split(chomp(line))[1]))
    y = append!(y,parse(split(chomp(line))[2]))
  end
end

aoa = linspace(-15,20,45)

cl,cd,cdp,cm,converged = Xfoil.xfoilsweep(x,y,aoa,400000.0,printdata=true,clminstop=true,clmaxstop=true)

aoa = aoa[find(converged)]
cl = cl[find(converged)]
cd = cd[find(converged)]
cdp = cdp[find(converged)]
cm = cm[find(converged)]

liftslope,zeroliftangle,aoafit,clfit = LiftProps.fitliftslope(aoa,cl)
_,clmax1 = LiftProps.findclmax(aoa,cl)
_,clmax2 = LiftProps.findclmax(aoa,cl,2)
_,clmax3 = LiftProps.findclmax(aoa,cl,1.0)
_,clmaxlinear1 = LiftProps.findclmaxlinear(aoa,cl,liftslope,zeroliftangle,interpolate=false)
_,clmaxlinear2 = LiftProps.findclmaxlinear(aoa,cl,liftslope,zeroliftangle,interpolate=true)
_,clmin1 = LiftProps.findclmin(aoa,cl)
_,clmin2 = LiftProps.findclmin(aoa,cl,2)
_,clmin3 = LiftProps.findclmin(aoa,cl,1.0)

# LIFT PLOT
plot(overwrite_figure = false)
myfont = font(12)
scatter!(aoa*180/pi,cl,label = "",markersize=5)
plot!(ylims = [min(minimum(cl),-0.1), maximum(cl)*1.1])
plot!(xlabel = "\$\\alpha\$ (degrees)",ylabel = "\$c_{l}\$")
plot!(legendfont = myfont,guidefont = myfont,tickfont = myfont)
# Plot XFOIL Data
scatter!(aoa*180/pi,cl,label = "XFOIL Data",markersize=5,markercolor=:red)
# Plot Linear Lift Slope Line
indclmin = indmin(abs.(cl-minimum(cl)))
liftslopelinex = linspace(min(zeroliftangle-2.5*pi/180,aoa[indclmin]),maximum(aoa),1000)
liftslopeliney = liftslope*(liftslopelinex-zeroliftangle)
minAoAidx = indmin(abs.(minimum(cl) - liftslopeliney))
maxAoAidx = indmin(abs.(maximum(cl) - liftslopeliney))
minAoA = liftslopelinex[minAoAidx]
maxAoA = liftslopelinex[maxAoAidx]
liftslopelinex = linspace(min(zeroliftangle-2.5*pi/180,minAoA),maxAoA,1000)
liftslopeliney = liftslope*(liftslopelinex-zeroliftangle)
plot!(liftslopelinex*180/pi,liftslopeliney,linewidth=2,
    label="Linear Lift Slope Line",linecolor=:blue)
# Plot Linear Lift Slope Line Data
scatter!(aoafit*180/pi,clfit,label = "Least Squares Data",
    markersize = 8,markercolor = :red,markershape = :x)
# Plot Zero Lift Angle
plot!([zeroliftangle*180/pi],seriestype = [:vline],label = "\$\\alpha_{0}\$",
    linestyle = :dashdot,linewidth = 1.5,linecolor = :brown)
# Plot clmax
plot!([clmax1],seriestype = [:hline],label = "\$c_{l,max,1}\$",
    linestyle = :dash,linewidth = 1,linecolor = :green)
plot!([clmax2],seriestype = [:hline],label = "\$c_{l,max,2}\$",
    linestyle = :dashdot,linewidth = 2,linecolor = :green)
plot!([clmax3],seriestype = [:hline],label = "\$c_{l,max,3}\$",
    linestyle = :dot,linewidth = 3,linecolor = :green)
plot!(show = true)
# Plot clmaxlinear
plot!([clmaxlinear1],seriestype = [:hline],label = "\$c_{l,max,linear,1}\$",
    linestyle = :dash,linewidth = 1,linecolor = :orange)
plot!([clmaxlinear2],seriestype = [:hline],label = "\$c_{l,max,linear,2}\$",
    linestyle = :dashdot,linewidth = 2,linecolor = :orange)
# Plot clmin
plot!([clmin1],seriestype = [:hline],label = "\$c_{l,min,1}\$",
    linestyle = :dash,linewidth = 1,linecolor = :purple)
plot!([clmin2],seriestype = [:hline],label = "\$c_{l,min,2}\$",
    linestyle = :dashdot,linewidth = 2,linecolor = :purple)
plot!([clmin3],seriestype = [:hline],label = "\$c_{l,min,3}\$",
    linestyle = :dot,linewidth = 3,linecolor = :purple)
plot!(show = true)
