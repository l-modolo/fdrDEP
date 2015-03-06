LPO = function(x, p=1, thread_number=1, max_alternative = 0.2, max_grid_size = 2187, min_previous_central_column = 0.5, v = F, real_density = NULL)
{
	# void crossvalidation(double** x, int* m, int* p, double* lambda, double* mu, double** score, int** D)
	m = length(x)
	lambda = 0
	mu = 1
	pi_0 = 0.9
	breaks = 0
	lambda_real = 0
	mu_real = 1
	pi_0_real = 0.9
	breaks_real = 0
	N_max = m/log(m)
	if(N_max > max_grid_size)
		N_max = max_grid_size
	steps = unlist(lapply(c(1:(log(N_max)/log(3))), function(x){3^x}))
	steps = ( steps*(steps+1) )/ 2
	steps = sum(steps)
	risk = rep(1000000, min(max_grid_size,floor(steps)))
	real_risk = rep(1000000, min(max_grid_size,floor(steps)))
	real_risk_computation = 0
	if(!is.null(real_density))
	{
		real_density = real_density[order(x)]
		real_risk_computation = 1
	}
	x = x[order(x)]
	res = .C('crossvalidation',
			as.double(x),
			as.integer(length(x)),
			as.integer(p),
			as.integer(thread_number),
			as.double(max_alternative),
			as.integer(max_grid_size),
			as.double(min_previous_central_column),
			as.double(real_density),
			as.integer(real_risk_computation),
			real_risk = as.double(real_risk),
			lambda_real = as.double(lambda_real),
			mu_real = as.double(mu_real),
			breaks_real = as.double(breaks_real),
			pi_0_real = as.double(pi_0_real),
			lambda = as.double(lambda),
			mu = as.double(mu),
			pi_0 = as.double(pi_0),
			breaks = as.integer(breaks),
			risk = as.double(risk),
			as.integer(floor(steps)))
	if(v)
	{
		cat("lambda = ",res$lambda," mu = ",res$mu," pi_0 = ",res$pi_0, "\n")
		hist(x, breaks=res$breaks, freq=F, main=paste("LPO with p=", p))
		abline(v=res$lambda, col='green', lwd=1)
		abline(v=res$mu, col='red', lwd=1)
		abline(h=res$pi_0, col="blue", lwd=1)
		plot(c(1:length(res$risk[res$risk < 1])), res$risk[res$risk < 1], type ="l", lwd=1, xlab="D", ylab="risk")
		abline(v=c(1:length(res$risk[res$risk < 1]))[ res$risk[res$risk < 1] == min(res$risk[res$risk < 1])], col="red")
	}
	if(real_risk_computation == 1)
	{
		return(list(lambda = res$lambda, mu = res$mu, pi_0 = res$pi_0, risk = res$risk[res$risk < 1000000], breaks=res$breaks, real_risk = res$real_risk[res$real_risk < 1000000], lambda_real = res$lambda_real, mu_real = res$mu_real, pi_0_real = res$pi_0_real, breaks_real=res$breaks_real))
	}
	else
		return(list(lambda = res$lambda, mu = res$mu, pi_0 = res$pi_0, risk = res$risk[res$risk < 1], breaks=res$breaks))
}