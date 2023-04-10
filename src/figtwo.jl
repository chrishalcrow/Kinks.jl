function f2makedata!(EL,x, dx, N, tanhdiff, sechlumps)
	
	thestart = 12.0
	thejump = -0.2
	theend = 1.0
	
	
	frames = Int((thestart - theend)/abs(thejump) + 1)
	println("there are ", frames, " frames.")
	
	phiL = zeros(frames,2,N)
	 
	a = 1
	for X0 = thestart:thejump:theend
	
		ti =100
		EList = zeros(ti)
		
		phiL[a,:,:] = makeinitial([tanhdiff,sechlumps],x, pars=true, X=X0)
		forceconstraint!(phiL,a, X0, [-0.5,0.5] ,dx,N)
		phiL[a,:,:] = flowfix!(phiL[a,:,:], (dx, N, EList, ti, 0.25, 1.0, getXi(X0,dx,N) ) )
		EL[a,1] = X0
		EL[a,2] = EList[ti]
		
		a += 1
		
		
	end
	
	return phiL 
		
end

function engplot(EL,Kpts,f2colors,phiL,x,therange;
	scale=3.0
	)

	xtks1 = [2.5,5.0,7.5,10.0]
	ytks1 = [0.58,0.59,0.60,0.61]
	
	xtks2=[-30,-15,0,15,30]
	ytks2=[-1.0,-0.5,0.0,0.5,1.0]
	
	xtks3=[-30,-15,0,15,30]
	ytks3=[0.0,0.2,0.4,0.6,0.8,1.0]

    spcolor = RGBf(0.7,0.7,0.7)

    f1 = Figure(backgroundcolor = :white, resolution = (6.161, 2.0) .* (72*scale), fontsize=10*scale)

    ax = Axis(f1[1,1],
		#title = L"\text{Energy}",
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale
    )
	
	if length(xtks1) != 0
		ax.xticks = (xtks1,[latexstring(xtks1[a]) for a in 1:length(xtks1)])
	end
	if length(ytks1) != 0
		ax.yticks = (ytks1,[latexstring(ytks1[a]) for a in 1:length(ytks1)])
	end
	

    xp = EL[:,1]
	ys = EL[:,2]

	lines!(ax, xp, ys, linewidth=1.5*scale)

	println(length(Kpts))

	for a in 1:length(Kpts)
		scatter!(EL[Kpts[a],1],EL[Kpts[a],2], marker = :circle, color=f2colors[a], markersize=6*scale)
	end

    ylims!(ax, 0.57, 0.62)
    xlims!(ax, 0.8,12.0)
	
    ax2 = Axis(f1[1,2],
		#title = L"\text{Energy}",
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale
    )
	
    xp = x

	ys = [ phiL[Kpts[a],1,:] for a in 1:length(Kpts) ] 

	ps = [ lines!(ax2, xp, ys[a], color=f2colors[a],linewidth=1.0*scale) for a in 1:length(Kpts) ] 
	
	if length(xtks2) != 0
		ax2.xticks = (xtks2,[latexstring(xtks2[a]) for a in 1:length(xtks2)])
	end
	if length(ytks2) != 0
		ax2.yticks = (ytks2,[latexstring(ytks2[a]) for a in 1:length(ytks2)])
	end
	
    ylims!(ax2, -1.1, 1.1)
    xlims!(ax2, -30,30)
	
    ax3 = Axis(f1[1,3],
		#title = L"\text{Energy}",
        bottomspinecolor = spcolor,
        topspinecolor = spcolor,
        leftspinecolor = spcolor,
        rightspinecolor = spcolor,
        xtickcolor = spcolor,
        ytickcolor = spcolor,
        titlesize = 12*scale,
    )
	
    xp = x

	ys = [ phiL[Kpts[a],2,:] for a in 1:length(Kpts) ] 

	ps = [ lines!(ax3, xp, ys[a], color=f2colors[a], linewidth=1.0*scale, label = latexstring(therange[Kpts[a]])) for a in 1:length(Kpts) ] 
	
	if length(xtks3) != 0
		ax3.xticks = (xtks3,[latexstring(xtks3[a]) for a in 1:length(xtks3)])
	end
	if length(ytks3) != 0
		ax3.yticks = (ytks3,[latexstring(ytks3[a]) for a in 1:length(ytks3)])
	end
	
    ylims!(ax3, -0.1,1.1)
    xlims!(ax3, -30,30)
	
	f1[1,4] = Legend(f1,
		ax3,
		L"X",
		    framevisible = false,	    
		    #[latexstring(therange[Kpts[a]]) for a=1:5],
	        padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0),
			labelhalign = :center
	    )
	
	#=
	 Legend(f1[1,4],
	ps,
	"hello",
	    framevisible = false,	    
	    [latexstring(therange[Kpts[a]]) for a=1:5],
        padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0),
		labelhalign = :center
    )
	
	=#
   
   	#colgap!(f1.layout,20*scale) 

    return f1
	
end