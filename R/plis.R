plis = function (LIS, fdr = 0.05, adjust = F) 
{
	n = length(LIS)
	s.LIS = sort(LIS)
	for (i in 1:n) {
		if (mean(s.LIS[1:i]) > fdr) 
		break
	}
	nNonNull = i - 1
	states = rep(0, n)
	if (nNonNull > 0) 
		states[LIS <= s.LIS[nNonNull]] = 1
	
	if (adjust) {
		aLIS = sapply(LIS, function(cut) mean(LIS[which(LIS <= cut)]))
		return(list(state = states, aLIS = aLIS))
	} else {
		return(states)
	}
}
