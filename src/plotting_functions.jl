

function plot_field(phi::Matrix{Float64}, x; 
	
	scale::Float64=1.0,
	xlimL = x[1],
	xlimR = x[end],
	ylimT = maximum(phi)*1.2,
	ylimB = minimum(phi)*1.2,
	titlename = "",
	xtks = [],
	ytks = []
	)
	
	println(xtks)

    spcolor = RGBf(0.7,0.7,0.7)

    f1 = Figure(backgroundcolor = :white, resolution = (6.161, 2.4) .* (72*scale), fontsize=10*scale)

    ax = Axis(f1[1,1],
		title = titlename,
		#title = titlename,
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale
    )
	
	if length(xtks) != 0
		ax.xticks = (xtks,[latexstring(xtks[a]) for a in 1:length(xtks)])
	end
	if length(ytks) != 0
		ax.yticks = (ytks,[latexstring(ytks[a]) for a in 1:length(ytks)])
	end

	
    xp = x
	ys = [ phi[1,:],phi[2,:] ]

	ps = [ lines!(ax, xp, ys[a], linewidth=1.5*scale) for a in 1:2 ] 

    ylims!(ax, ylimB, ylimT)
    xlims!(ax, xlimL, xlimR)



	Legend(f1[1,2],
	    framevisible = false,
	    ps,
	    [L"\Phi_1(x)",L"\Phi_2(x)"],
        padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0),
		linepoints = [Point2f(1.0-scale, 0.5), Point2f(1.0, 0.5)]
    )
	

    #ax.yticks  = ([-1.0, -0.5, 0.0, 0.5, 1.0],[ L"-1.0", L"-0.5", L"0.0", L"0.5", $1\.0$])
   #= #ax.xticks  = ([-20.0, -10.0, 0.0, 10.0, 20.0],[L"-20.0", L"-10.0", L"0.0", L"10.0", L"20.0"])
    =#
   
   	colgap!(f1.layout,20*scale) 

    return f1
	
end


function plot_heat(phi, N; scale=1.0)
	
	xc = -1.3:0.01:1.3
	yc =  -2.0:0.01:2.9
	VN = zeros(size(xc)[1],size(yc)[1])

	for a in 1:size(xc)[1], b in 1:size(yc)[1]
	    VN[a,b] = min(0.11,V(xc[a],yc[b],0.25,1.0))
	end

	f = Figure(backgroundcolor = :white, resolution = (0.5*6.161, 0.5*6.161) .* (72*scale), fontsize=10*scale)
	ax = Axis(f[1,1],
	    title = L"V(\Psi_1, \Psi_2)",
	    xlabel = L"\Phi_1",
	    ylabel = L"\Phi_2"
	    #,protrusions = 0
	)
	
	minL = 0.002
	maxL = 0.122
	nlev = 16
	
	hm = contour!(ax,xc,yc,VN,
	    colormap = Reverse(:rainbow1),
	    levels= minL:(maxL-minL)/nlev:maxL,
	    linewidth = 0.5*scale
	)
		
	colb = Colorbar(f[1,2],hm
	    ,size = 5*scale
	   # vertical=false
	    ,ticks = ([0.0, 0.05, 0.1],[L"0.0", L"0.05", L"0.1"])
	    ,tellheight = true
	    ,tellwidth = true
	)
	
	p1 = lines!(ax,[ Point(phi[1,n],phi[2,n]) for n in 1:N ],
	    linewidth = 1*scale,
	    color=:black
	)
	
	scatter!(0.0,phi[2,Int(round(N/2))], marker = :rtriangle, color=:black, markersize=9*scale)
	
		
		return f
		
		

	#=
    


	p1 = lines!(ax,[ Point(ϕN[1,n],ϕN[2,n]) for n in 1:400 ],
	    linewidth = 1,
	    color=:black
	)
    
	#ax.yticks  = ([-1.0, -0.5, 0.0, 0.5, 1.0],[L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"])
	#ax.xticks  = ([-20.0, -10.0, 0.0, 10.0, 20.0],[L"-20.0", L"-10.0", L"0.0", L"10.0", L"20.0"])

	ax.yticks = ([-2.0, 0.0, 2.0],[L"-2.0", L"0.0", L"2.0"])
	ax.xticks  = ([-1.0, -0.5, 0.0, 0.5, 1.0],[L"-1.0", L"-0.5", L"0.0", L"0.5", L"1.0"])
    

    
	scatter!(0.0,0.49, marker = :rtriangle, color=:black, markersize=6)

    
	
	colgap!(f.layout, 5)=#
	
end




