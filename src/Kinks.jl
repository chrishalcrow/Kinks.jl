module Kinks

using Makie
using GLMakie
using CairoMakie
using DifferentialEquations
using DiffEqCallbacks
using LaTeXStrings

export setgrid
export makeinitial
export plot_field
export flow!, flowfix!, getXi
export plot_heat
export forceconstraint!
export f2makedata!, engplot
export f3makedata, bindplot

include("plotting_functions.jl")
include("diff_functions.jl")
include("initialise_functions.jl")
include("figtwo.jl")
include("figthree.jl")


function f3makedata(x,dx, N,muN,mN)
	
	#x = -40.0:0.1:40.0
	#dx=0.1
	
	dn = muN/100.0

	muL = [ dn*(i) for i=1:muN ]
	mL = [ dn*(i) for i=1:mN ]

	E11 = zeros(muN,mN)
	E22 = zeros(muN,mN)
	
	phi = zeros(2,N)

	ti = 100
	EList = zeros(ti)

	tspan = (0.0,ti+0.0)
	engtimes = [ 1.0*i for i=1:ti ]


	for a in 1:muN, b in 1:mN

	    #EList = zeros(ti)
	    p = (dx, N, EList, ti, muL[a], mL[b])

		# 2 kinks
	    phi[1,:] =  [ ( tanh(x[n]+4.0) + tanh(x[n]-4.0) )/2.0 for n=1:length(x)] ;
	    phi[2,:] =  [ mL[b]*0.5*( 1.0/cosh(x[n]+4.0)^2 + 1.0/cosh(x[n]-4.0)^2 ) for n=1:length(x)] ;
 
	    cb = PresetTimeCallback(engtimes, affect_checkeng!, save_positions=(false,false))
	    prob = ODEProblem(dEdu!,phi,tspan,p)
	    sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=false, callback=cb, tstops = engtimes )
		
		#export("plot(sol[end])
		#save("jfigs/2fig_"*string(a)*"_"*string(b)*".png", plot_field(sol[end],x));
		
		E22[a,b] = EList[end]
		
		# 1 kink
	    phi[1,:] =  [ ( tanh(x[n]) - 1.0 )/2.0 for n=1:length(x)] ;
	    phi[2,:] =  [ mL[b]*( tanh(x[n]) + 1.0 )*0.5 for n=1:length(x)] ;

	    cb = PresetTimeCallback(engtimes, affect_checkeng!, save_positions=(false,false))
	    prob = ODEProblem(dEdu!,phi,tspan,p)
	    sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=false, callback=cb, tstops = engtimes )

		#save("jfigs/1fig_"*string(a)*"_"*string(b)*".png", plot_field(sol[end],x));

		E11[a,b] = EList[end]
		
   end
   
   return (2.0*E11 .- E22 )./(2.0.*E11)
    
end

function bindplot(ED,muN,mN;
	scale = 1.0) 
	
	dn = muN/100.0
	
	print(dn)
	
	muL = [ dn*(i) for i=1:muN ]
	mL = [ dn*(i) for i=1:mN ]
	
	CairoMakie.activate!()
	
	#xtks = []
	#ytks = []
	xtks = [0.5,1.0,1.5,2.0]
	ytks = [0.2,0.4,0.6,0.8,1.0]
	
	#muL = [ 0.2*(i) for i=1:5 ]
	# mL = [ 0.2*(i) for i=1:10 ]
	
    spcolor = RGBf(0.7,0.7,0.7)

    f1 = Figure(backgroundcolor = :white, resolution = (5.0, 2.5) .* (72*scale), fontsize=10*scale)

    ax = Axis(f1[1,1],
		#title = titlename,
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale,
		xlabel = L"m",
		ylabel = L"\mu",
		xgridvisible= false,
		ygridvisible= false
		#xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)]),
		#yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
		#xticks  = ([-5.0, -2.5, 0.0, 2.5, 5.0],[ L"-5.0", L"-2.5", L"0.0", L"2.5", L"5.0"]),
		#yticks  = ([-1.0, -0.5, 0.0, 0.5, 1.0],[ L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"])
        #xlabel = L"x"
    )
	
	if length(xtks) != 0
		ax.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
	end
	if length(ytks) != 0
		ax.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
	end
	
	thelevels =
	
	contourf!(ax, mL, muL, ED',colormap=:rainbow1,levels=[i for i=0:0.02:0.2])
	
	contour!(ax, mL, muL, ED', grid=false, labels = true, levels=[i for i=0:0.04:0.3],colormap=:rainbow1,
    labelsize = 7,  labelcolor = :black, labelfont = "courier")# colormap=:plasma, levels=10, linewidth=scale*1.0)
	
	#levels=[i for i=0:0.04:0.3]

	ylims!(ax, 0.1, 1.1)
    xlims!(ax, 0.4, 2.0)

	#colb = Colorbar(f1[1,2],cp
	   # ,size = 5*scale
	   # vertical=false
	   # ,ticks = ([0.0, 0.05, 0.1],[L"0.0", L"0.05", L"0.1"])
	  #  ,tellheight = true
	  #  ,tellwidth = true
	#)
	
#    xp = x
#	ys = [ phi[1,:],phi[2,:] ]

#	ps = [ lines!(ax, xp, ys[a], linewidth=1.5*scale) for a in 1:2 ] 

 #   ylims!(ax, ylimB, ylimT)
 #   xlims!(ax, xlimL, xlimR)



#	Legend(f1[1,2],
#	    framevisible = false,
#	    ps,
#	    [L"\Phi_1(x)",L"\Phi_2(x)"],
 #       padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0),
#		linepoints = [Point2f(1.0-scale, 0.5), Point2f(1.0, 0.5)]
 #   )
	

    #ax.yticks  = ([-1.0, -0.5, 0.0, 0.5, 1.0],[ L"-1.0", L"-0.5", L"0.0", L"0.5", $1\.0$])
   #= #ax.xticks  = ([-20.0, -10.0, 0.0, 10.0, 20.0],[L"-20.0", L"-10.0", L"0.0", L"10.0", L"20.0"])
    =#
   
   	#colgap!(f1.layout,20*scale) 

    return f1
	#save("jfigs/bindplot.png", );

	
end


function setgrid(N,dx)
	
	println("grid has ", N, " points and goes from ", -0.5*dx*(N-1), " to ", 0.5*dx*(N-1))	
	return -0.5*dx*(N-1):dx:0.5*dx*(N-1)
	
end


function V(ϕ1,ϕ2,μ,m)
    0.5*ϕ1^2*(1.0 - ϕ1^2)^2 + 0.5*μ^2*(1 - ϕ1^2 - ϕ2/m)^2
end

function d1V(ϕ1,ϕ2,μ,m)
     ϕ1 - 4.0*ϕ1^3 + 3.0*ϕ1^5 + (2.0*μ^2*ϕ1*( m*(ϕ1^2-1.0) + ϕ2))/m 
end

function d2V(ϕ1,ϕ2,μ,m)
    ( m*(ϕ1^2 - 1.0) + ϕ2)*(μ^2/m^2)
end

function d11V(ϕ1,ϕ2,μ,m)
    1.0 - 12.0*ϕ1^2 + 15.0*ϕ1^4 + (2.0*μ^2*(m*(-1.0 + 3.0*ϕ1^2) + ϕ2))/m
end

function d12V(ϕ1,ϕ2,μ,m)
    (2.0*μ^2*ϕ1)/m
end

function d22V(ϕ1,ϕ2,μ,m)
    μ^2/m^2
end

function dxB(phi, a, i, ls)
     @fastmath @inbounds 0.5*(phi[a,i+1] - phi[a,i-1])/ls
end

function d2xB(phi, a, i, ls)
    @fastmath @inbounds (phi[a,i+1] - 2.0*phi[a,i] + phi[a,i-1])/ls^2
end

function dxD(phi, a, i,ls)
    (phi[a,i+2] - 8.0*phi[a,i+1] + 8.0*phi[a,i-1] - phi[a,i-2] )/(12.0*ls)
end

function d2xD(phi, a,i,ls)
    (-phi[a,i+2] + 16.0*phi[a,i+1] - 30.0*phi[a,i] + 16.0*phi[a,i-1] - phi[a,i-2] )/(12.0*ls^2)
end

function calcEng(u,lp,ls,μ,m)
    eng = 0.0
    for i in 3:lp-2
        eng += V(u[1,i],u[2,i],μ,m) + 0.5*(dxD(u,1,i,ls)^2 + dxD(u,2,i,ls)^2)
    end
    return eng*ls
end






end
