// Model 2 - Multimodal prior model
//
// ---------------------------------------------------------------------------------------------------------------
//
// Pain experience depends not only on the pain stimulation, but also on the information in the cue, 
// where that information is represented as a multi-modal ‘prior’ distribution.
// The distribution consists of a weighted average of normal distributions, modulated by the parameter . The
// normals are centered at each possible magnitude (v_1ij … v_mij) of the incoming pain stimulus and are weighted by
// their corresponding probabilities of occurrence, (w_1ij … w_mij).
// 
// P(R|X,Z) ~ Noraml(X, beta) * SUM [ w_kij * Normal(v_kij, rho)
// 
// ---------------------------------------------------------------------------------------------------------------

//Functions which calculate the means and standard deviation of the product of Normal Distriubtions
functions {
  
  //Product of 2 Normal Distributions
  real get_mean(real m1, real m2, real s1, real s2){
    return (m1*s2*s2 + m2*s1*s1)/ (s1*s1 + s2*s2);}
  real get_std(real m1, real m2, real s1, real s2){
    return sqrt(s1*s1*s2*s2/(s1*s1+s2*s2));}
    
  //Product of 3 Normal Distributions
  real get_mean3(real m1, real m2, real m3, real s1, real s2, real s3){
    return get_mean(m1, get_mean(m2,m3,s2,s3), s1, get_std(m2,m3,s2,s3));}
  real get_std3(real m1, real m2, real m3, real s1, real s2, real s3){
    return get_std(m1, get_mean(m2,m3,s2,s3), s1, get_std(m2,m3,s2,s3));}
}

data {
  int<lower=0> N;
  int subject[N];                   //Subject ID
  vector[N] q;                      //Mean of Presented Que
  matrix[N, 4] v;                   //Presented Magnitudes
  matrix[N, 4] w;                   //Magnitude Weights
  vector[N] std;                    //Standard Deviation of Presented Que
  vector[N] X;                      //Administerd Stimulus
  vector<lower=0, upper=100>[N] r;  //Subject Rating
} 

transformed data{
  int<lower=0> S = max(subject);
}

parameters {
  real<lower=0, upper=100> beta;
  real<lower=0, upper=100> rho;
}

transformed parameters {
  
  //We assign the appropriate parameters to each trial
  vector<lower=0>[N] beta_n;
  vector<lower=0>[N] rho_n;
  
  //mean_p and std_p are the Mean and Standard Deviation of the resulting 'prior distribution'
  vector[N] mean_p;
  vector[N] std_p;
  
  for (n in 1:N){
    int s = subject[n];
    beta_n[n] = beta;
    rho_n[n] = rho;
    
    //the 'prior distribution' results from the product of two Normal terms: A stimulus term and a que term
    mean_p[n] = get_mean(X[n], q[n], beta_n[n], rho_n[n]);
    std_p[n] = get_std(X[n], q[n], beta_n[n], rho_n[n]);
  }
}


model {
  
  //Weak priors on the parameters make the simulations more stable
  //beta ~ cauchy(0,100);
  //rho ~ cauchy(0,100);
  
  //r ~ normal(X, beta_n);
  //r ~ normal(q, rho_n);
  
  for (n in 1:N){
    
    real m = 0.0;
    real prob = 0.0;
    real prior = 0.0;
    real prob_integral = 0.0;
    real prior_integral = 0.0;
    
    for (i in 1:4){
      m = m + v[n,i]*w[n,i];
    }
    
    //Calculate Normalization Constant
    for (theta_i in 1:100){
      prior_integral = 0.0;
      for (i in 1:4){
        prior_integral = prior_integral + w[n,i] * exp(normal_lpdf(theta_i | v[n,i], rho_n[n]));
      }
      prob_integral = prob_integral + exp(normal_lpdf(theta_i | X[n], beta_n[n])) * prior_integral;
    }
    
    prior = 0;
    for (i in 1:4){
      prior = prior + w[n,i] * exp(normal_lpdf(r[n] | v[n,i], rho_n[n]));
    }
    
    prob = exp(normal_lpdf(r[n] | X[n], beta_n[n])) * prior;
    prob = prob / prob_integral;
    
    target += log(prob);
  } 
  

}

generated quantities {
  real log_lik_t;         //Log Likelihood of the experiment with given parameters 
  vector[S] log_lik_s;    //Log Likelihood of each subject
  vector[N] log_lik;      //Log Likelihood of each trial
  
  {
  log_lik_t = 0.0;
  for (s in 1:S){log_lik_s[s] = 0.0;}
  
  for (n in 1:N){
    int s = subject[n];
    real prob = 0.0;
    real prior = 0.0;
    real prob_integral = 0.0;
    real prior_integral = 0.0;

    //Calculate Normalization Constant
    for (theta_i in 1:100){
      prior_integral = 0.0;
      for (i in 1:4){
        prior_integral = prior_integral + w[n,i] * exp(normal_lpdf(theta_i | v[n,i], rho_n[n]));
      }
      prob_integral = prob_integral + exp(normal_lpdf(theta_i | X[n], beta_n[n])) * prior_integral;
    }
    
    prior = 0.0;
    for (i in 1:4){
      prior = prior + w[n,i] * exp(normal_lpdf(r[n] | v[n,i], rho_n[n]));
    }
    
    prob = exp(normal_lpdf(r[n] | X[n], beta_n[n])) * prior;
    prob = prob / prob_integral;
    
    log_lik[n] = log(prob);
    log_lik_s[s] = log_lik_s[s] + log_lik[n];
    log_lik_t = log_lik_t + log_lik[n];
  }  
  }
}
