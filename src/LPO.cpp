#include <R.h>
#include <Rdefines.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include <iostream>
#include <ostream>


using namespace std;

extern "C" {

double pi_0_computation(double* x, 
						int m, 
						double lambda, 
						double mu)
{
	int i = 0;
	for(int j = 0; j < m; j++)
	{
		if(lambda <= x[j] && x[j] <= mu)
			i++;
	}
	return(((double)i)/((double)m * (mu - lambda)));
}

void irregular_histogram(double* x, // the ordered vector of p-values
						int m, // size of x
						int breaks, // the number of breaks of the fine grid
						int* D, // the dimention of the model (number of breaks).
						double lambda, // lambda the left limit of the uniform region
						double mu,  // the right limit of the uniform region
						int* m_k, // the number of x in interval k
						double* omega_k, // the width of interval k
						double* segment_k) // the right limit of interval k
{
	for(int i = 0; i < breaks; i++)
	{
		segment_k[i] = 0.0;
		omega_k[i] = 0.0;
		m_k[i] = 0.0;
	}
	*D = 0;
	double segment = 0.0;
	double segment_size = 1.0 / ((double)breaks);
	int j = 0;
	for(int i = 1; i <= breaks; i++)
	{
		segment = segment_size * ((double)i);
		omega_k[(*D)] += 1.0 / ((double)breaks);
		if(segment <= lambda || mu <= segment)
		{
			segment_k[(*D)] = segment;
			while(x[j] <= segment_k[(*D)] && j < m)
			{
				m_k[(*D)] = m_k[(*D)] + 1;
				j++;
			}
			*D = *D + 1;
		}
	}
}

double risk(int m_tmp, 
			int D_tmp, 
			int p_tmp, 
			int* m_k, 
			double* omega_k)
{
	// score = ((2*m-p)/((m-1)*(m-p))) * sum( (m_k/(m*omega_k)) )
	// score = score - ( ( (m*(m-p+1))/((m-1)*(m-p)) ) * sum( (1/omega_k)*((m_k/m)^2) ) )
	double m = (double)m_tmp;
	double D = (double)D_tmp;
	double p = (double)p_tmp;
	// for q given p we compute the risk of a given irregular histogram
	double score = (2.0*m-p)/((m-1.0)*(m-p));
	double indice1 = 0.0;
	double indice2 = 0.0;
	for(int i = 0; i < D; i++)
		indice1 += ((double)m_k[i])/(m*omega_k[i]);
	score = score * indice1;
	indice1 = (m*(m-p+1.0))/((m-1.0)*(m-p));
	for(int i = 0; i < D; i++)
		indice2 += (1.0/omega_k[i])*((((double)m_k[i])/m)*(((double)m_k[i])/m));
	return(score - (indice1 * indice2));
}

double rrisk(double* x, // the ordered vector of p-values
				int m, // size of x
				int D, // the dimention of the model (number of breaks).
				int* m_k, // the number of x in interval k
				double* omega_k, // the width of interval k
				double* segment_k, // the right limit of interval k
				double* real_density)
{
	int i = 0;
	double real_risk = 0.0;
	for(int j = 0; j < D; j++)
	{
		while(i < m && x[i] < segment_k[j])
		{
			// cout << real_density[i] << endl;
			// cout << x[i] << "/" << segment_k[j] << endl;
			real_risk += (( ( ((double)m_k[j])/(omega_k[j] * (double)m) ) - real_density[i] ) * ( ( ((double)m_k[j])/(omega_k[j] * (double)m) ) - real_density[i] ));
			i++;
		}
	}
	return(-sqrt(real_risk));
}

void crossvalidation(double* x, 
					int* m, 
					int* p, 
					int* thread_mumber, 
					double* max_alternative, 
					int* max_grid_size, 
					double* min_previous_central_column, 
					double* real_density, 
					int* real_risk_computation,
					double* real_risk, 
					double* lambda_real,
					double* mu_real,
					double* breaks_real,
					double* pi_0_real,
					double* lambda, 
					double* mu, 
					double* pi_0, 
					int* breaks, 
					double* score, 
					int score_size)
{
	// the results variables
	double best_risk = 1000000.0;
	double best_lambda = 0.0;
	double best_mu = 0.0;
	int best_grid = 0;

	double best_risk_real = 1000000.0;
	double best_lambda_real = 0.0;
	double best_mu_real = 0.0;
	int best_grid_real = 0;

	// we compute the grids list to explore
	int grid_max = (int)floor( (double)*m / log((double)*m) );
	if(grid_max > *max_grid_size)
		grid_max = *max_grid_size;
	int grids_size = (int)floor( log((double)grid_max) / log(3.0) - 1);
	int* grids = new int[grids_size];
	for(int i = 0; i < grids_size - 1; i++)
		grids[i] = pow(3.0, ((double)i)+3.0);
	grids[grids_size-1] = grid_max;

	int* D = NULL;
	double r = 0.0;
	double r_real = 0.0;
	int* m_k = NULL;
	double* omega_k = NULL;
	double* segment_k = NULL;
	double new_lambda = 0.0;
	double new_mu = 1.0;
	int min_central = 2;
	int max_j = (int)round(*max_alternative * (double)grids[0]);
	omp_set_num_threads(*thread_mumber);

	#pragma omp parallel shared(x, p, m, max_alternative, max_j, min_previous_central_column, best_risk, best_lambda, best_mu, best_grid, score, grids, grids_size, min_central) private(D, r, m_k, omega_k, segment_k, new_lambda, new_mu)
	{
		D = new int;
		*D = 0;
		r = 0.0;
		r_real = 0.0;
		new_lambda = 0.0;
		new_mu = 1.0;
		m_k = new int[grid_max];
		omega_k = new double[grid_max];
		segment_k = new double[grid_max];
		// given the maximal size of grid we know the size of the following vector
		for(int i = 0; i < grids_size; i++)
		{
			#pragma omp for schedule(dynamic, 20)
			for(int j = 2; j < max_j - min_central - 1; j++)
			{
				new_lambda = 1.0/((double)grids[i]) * (double)j;
				for(int k = j + min_central; k < grids[i] - 2; k++)
				{
					new_mu = 1.0/((double)grids[i]) * (double)k;
					irregular_histogram(x, *m, grids[i], D, new_lambda, new_mu, m_k, omega_k, segment_k);
					r = risk(*m, *D, *p, m_k, omega_k);
					#pragma omp critical(dataupdate)
					{
						if(score[(*D)] > r)
 							score[(*D)] = r;
						if(score[(*D)] < best_risk)
						{
							best_risk = score[(*D)];
							best_lambda = new_lambda;
							best_mu = new_mu;
							best_grid = grids[i];
						}
					}
				}
			}
			#pragma omp barrier
			#pragma omp master
			{
				if(i < grids_size)
					min_central = (int)round((best_mu - best_lambda) * (double)grids[i+1] * *min_previous_central_column);
				if(min_central < 1)
					min_central = 2;
				max_j = (int)round(*max_alternative * (double)grids[i+1]);
			}
			#pragma omp barrier
		}
		delete D;
		delete[] m_k;
		delete[] omega_k;
		delete[] segment_k;
	}
	if(*real_risk_computation == 1)
	{
		min_central = 2;
		max_j = (int)round(*max_alternative * (double)grids[0]);
		#pragma omp parallel shared(x, p, m, max_alternative, max_j, min_previous_central_column, best_risk_real, best_lambda_real, best_mu_real, best_grid_real, real_risk, grids, grids_size, min_central, real_risk_computation) private(D, r, m_k, omega_k, segment_k, new_lambda, new_mu)
		{
			D = new int;
			*D = 0;
			r = 0.0;
			r_real = 0.0;
			new_lambda = 0.0;
			new_mu = 1.0;
			m_k = new int[grid_max];
			omega_k = new double[grid_max];
			segment_k = new double[grid_max];
			// given the maximal size of grid we know the size of the following vector
			for(int i = 0; i < grids_size; i++)
			{
				#pragma omp for schedule(dynamic, 20)
				for(int j = 2; j < max_j - min_central - 1; j++)
				{
					new_lambda = 1.0/((double)grids[i]) * (double)j;
					for(int k = j + min_central; k < grids[i] - 2; k++)
					{
						new_mu = 1.0/((double)grids[i]) * (double)k;
						irregular_histogram(x, *m, grids[i], D, new_lambda, new_mu, m_k, omega_k, segment_k);
						r_real = rrisk(x, *m, *D, m_k, omega_k, segment_k, real_density);
						#pragma omp critical(dataupdate)
						{
							if(real_risk[(*D)] > r_real)
								real_risk[(*D)] = r_real;
							if(real_risk[(*D)] < best_risk_real)
							{
								best_risk_real = real_risk[(*D)];
								best_lambda_real = new_lambda;
								best_mu_real = new_mu;
								best_grid_real = grids[i];
							}
						}
					}
				}
				#pragma omp barrier
				#pragma omp master
				{
					if(i < grids_size)
						min_central = (int)round((best_mu_real - best_lambda_real) * (double)grids[i+1] * *min_previous_central_column);
					if(min_central < 1)
						min_central = 2;
					max_j = (int)round(*max_alternative * (double)grids[i+1]);
				}
				#pragma omp barrier
			}
			delete D;
			delete[] m_k;
			delete[] omega_k;
			delete[] segment_k;
		}
	}

	delete[] grids;
	*lambda = best_lambda;
	*mu = best_mu;
	*breaks = best_grid;
	*pi_0 = pi_0_computation(x, *m, *lambda, *mu);
	if(*real_risk_computation == 1)
	{
		*lambda_real = best_lambda_real;
		*mu_real = best_mu_real;
		*breaks_real = best_grid_real;
		*pi_0_real = pi_0_computation(x, *m, *lambda_real, *mu_real);
	}
}

}