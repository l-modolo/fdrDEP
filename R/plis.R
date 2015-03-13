plis = function (LIS, fdr = 0.05, adjust = F) 
{
	m = length(LIS)
	s.LIS = sort(LIS)
	for (i in 1:m) {
		if (mean(s.LIS[1:i]) > fdr) 
		break
	}
	nNonNull = i - 1
	states = rep(0, m)
	if (nNonNull > 0) 
		states[LIS <= s.LIS[nNonNull]] = 1
	
	if (adjust) {
		aLIS = sapply(LIS, function(cut) mean(LIS[which(LIS <= cut)]))
		return(list(state = states, aLIS = aLIS))
	} else {
		return(states)
	}
}

bh = function(pvalues, fdr = 0.05, mu = 1, pi_0 = 1)
{
	m = length(pvalues)
	n = length(pvalues[pvalues >= mu])
	if(mu < 1)
	{
		exess = n - pi_0 * m *(1-mu)
		pi_0 = 1
	}
	else
		exess = 0
	s.pvalues = sort(pvalues)
	for (i in 1:m) {
		if ( s.pvalues[i] > (i/(m-exess)) * (fdr/pi_0))
		break
	}
	nNonNull = i - 1
	states = rep(0, m)
	if (nNonNull > 0) 
		states[pvalues <= s.pvalues[nNonNull]] = 1
	return(states)
}