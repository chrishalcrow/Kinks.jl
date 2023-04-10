function makeinitial(fl, x; pars = false, X)
	if pars == false
		return [ fl[a](x[n]) for a=1:2, n=1:size(x)[1] ]
	else
		return [ fl[a](x[n],X) for a=1:2, n=1:size(x)[1] ]
	end
	#f(x[n]) for n=1:size(x)[1] ]
end

#function makeinitial(fl, x)
#	return [ fl[a](x[n]) for a=1:2, n=1:size(x)[1] ]
	#f(x[n]) for n=1:size(x)[1] ]
	#end