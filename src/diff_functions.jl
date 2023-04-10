


function dEdu!(du,u,p,t)
    @inbounds for i in 3:p[2]-2
        du[1,i] = d2xD(u, 1, i, p[1]) - d1V(u[1,i],u[2,i],p[5],p[6])
        du[2,i] = d2xD(u, 2, i, p[1]) - d2V(u[1,i],u[2,i],p[5],p[6])
    end
end

function dEdufix!(du,u,p,t)
    @inbounds for i in 3:p[2]-2
        du[1,i] = d2xD(u, 1, i, p[1]) - d1V(u[1,i],u[2,i],p[5],p[6])
        du[2,i] = d2xD(u, 2, i, p[1]) - d2V(u[1,i],u[2,i],p[5],p[6])
    end
	du[1, Int(p[7][1])] = 0.0
	du[1, Int(p[7][2])] = 0.0
end


function affect_checkeng!(integrator)
    
    n = floor(Int,round(integrator.t))
    
    integrator.p[3][n] = calcEng(integrator.u,integrator.p[2],integrator.p[1],integrator.p[5],integrator.p[6])
    
    if n > 3 && abs(integrator.p[3][n] - integrator.p[3][n-1]) < 10.0^(-10)
        terminate!(integrator)
        integrator.p[3][(n+1):integrator.p[4]] .= integrator.p[3][n]
    end

end



function flow!(phi,p)
    
	ti = p[4]
	
    EList = zeros(ti)

    #(dx, N, EList, ti, μ, m) = p


    tspan = (0.0,ti+0.0)
    engtimes = [ 1.0*i for i=1:ti ]

    cb = PresetTimeCallback(engtimes, affect_checkeng!, save_positions=(false,false))

    prob = ODEProblem(dEdu!,phi,tspan,p)
	sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=false, callback=cb )
    #sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=true, callback=cb, tstops = engtimes )

    if sol.t[end] == ti 
        println("AGGGG");
    else
        println(sol.t[end])
    end

    return sol[end]
    
end

function getXi(X,dx,N)
	x = -(N-1)*dx*0.5:dx:(N-1)*dx*0.5
	return [ findall(alp->alp==-X, x)[1], findall(alp->alp==X, x)[1] ]
end
	
function forceconstraint!(phi::Matrix{Float64},X,theconstraint,dx,N)
	phi[1,getXi(X,dx,N)] = theconstraint
end
function forceconstraint!(phi::Array{Float64,3},a,X,theconstraint,dx,N)
	phi[a,1,getXi(X,dx,N)] = theconstraint
end
	
## assuming that X is positive
function flowfix!(phi,p)
    
	ti = p[4]
	
    EList = zeros(ti)
	
    tspan = (0.0,ti+0.0)
    engtimes = [ 1.0*i for i=1:ti ]

    cb = PresetTimeCallback(engtimes, affect_checkeng!, save_positions=(false,false))

    prob = ODEProblem(dEdufix!,phi,tspan,p)
	sol = solve(prob, Tsit5(),reltol=1e-8, abstol=1e-8, save_everystep=false, callback=cb )

    if sol.t[end] == ti 
        println("AGGGG");
    else
        println(sol.t[end])
    end
	
	println("The values ", sol[end][1, p[7][1]], " and ", sol[end][1, p[7][1]])

    return sol[end]
    
end

