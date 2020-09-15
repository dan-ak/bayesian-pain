// 
// Model 6 - Variable full Bayesian model - MLE
//
// ---------------------------------------------------------------------------------------------------------------
//
// In this model pain experience is determined jointly by the pain stimulation and two tiers of expectations.
// In this model parameters are allowed to vary by participant. 
// 
// P(R|X,Z) ~ Noraml(X, beta_i) * Normal(q, rho_i + eta_i*sd) * Normal(mu_i, nu_i)
// 
// ---------------------------------------------------------------------------------------------------------------
//
// The MLE version takes in the MLE hyperparameters and returns likelihood samples conditioned on the hyperparameters
// The samples can be used as an estimate of the marginalizle MLE likelihood. 
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
  vector[N] std;                    //Standard Deviation of Presented Que
  vector[N] X;                      //Administerd Stimulus
  vector<lower=0, upper=100>[N] r;  //Subject Rating

  real<lower=0, upper=80> mu_u;
  real<lower=0, upper=100> rho_u;
  real<lower=0, upper=100> nu_u;
  real<lower=0, upper=100> beta_u;
  real<lower=-10, upper=10> eta_u;
  
  real<lower=0, upper=100> mu_s;
  real<lower=0, upper=100> rho_s;
  real<lower=0, upper=100> nu_s;
  real<lower=0, upper=100> beta_s;
  real<lower=0, upper=100> eta_s;
} 

transformed data{
  int<lower=0> S = max(subject);
}

parameters {
  
  vector<lower=0, upper=100>[S] beta;
  vector<lower=0, upper=100>[S] rho;
  vector<lower=-10, upper=10>[S] eta;
  vector<lower=0, upper=100>[S] mu;
  vector<lower=0, upper=100>[S] nu;
}

transformed parameters {
  
  //We assign the appropriate parameters to each trial
  vector[N] beta_n;
  vector[N] rho_n;
  vector[N] eta_n;
  vector[N] mu_n;
  vector[N] nu_n;
  
  //mean_p and std_p are the Mean and Standard Deviation of the resulting 'prior distribution' on each trial
  vector[N] mean_p;
  vector[N] std_p;
  
  for (n in 1:N){
    int s = subject[n];
    beta_n[n] = beta[s];
    rho_n[n] = rho[s];
    eta_n[n] = eta[s];
    mu_n[n] = mu[s];
    nu_n[n] = nu[s];
    
    mean_p[n] = get_mean3(X[n], q[n], mu_n[n], beta_n[n], rho_n[n] + eta_n[n]*std[n], nu_n[n]);
    std_p[n] = get_std3(X[n], q[n], mu_n[n], beta_n[n], rho_n[n] + eta_n[n]*std[n], nu_n[n]);
  }
}
model {
  
  beta ~ normal(beta_u, beta_s);
  rho ~ normal(rho_u, rho_s);
  eta ~ normal(eta_u, eta_s);
  mu ~ normal(mu_u, mu_s);
  nu ~ normal(nu_u, nu_s);
  
  r ~ normal(mean_p, std_p);

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

    log_lik[n] = normal_lpdf(r[n] | mean_p[n], std_p[n]);
    log_lik_s[s] = log_lik_s[s] + log_lik[n];
    log_lik_t = log_lik_t + log_lik[n];
  }  
  }
}

