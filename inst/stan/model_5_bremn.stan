// 
// Model 5 - Full Bayesian model
//
// ---------------------------------------------------------------------------------------------------------------
//
// In this model pain experience is determined jointly by the pain stimulation and two tiers of expectations: 
// 1. Stable cue-independent effects akin to optimism/pessimism trait and effects informed by the cue. 
// 2. Cue-independent effects are represented by a probability distribution with mean mu and variance nu, 
// which is integrated with the information provided by the cue to form the final prior distribution. 
// 
// P(R|X,Z) ~ Noraml(X, beta) * Normal(q, rho + eta*sd) * Normal(mu, nu)
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
} 

transformed data{
  int<lower=0> S = max(subject);
}

parameters {
  real<lower=0, upper=100> beta;
  real<lower=0, upper=100> rho;
  real<lower=-1, upper=10> eta;
  real<lower=0, upper=100> mu;
  real<lower=0, upper=100> nu;
}

transformed parameters {
  
  //We assign the appropriate parameters to each trial
  vector<lower=0>[N] beta_n;
  vector<lower=0>[N] rho_n;
  vector<lower=0>[N] eta_n;
  vector<lower=0>[N] mu_n;
  vector<lower=0>[N] nu_n;
  
  //mean_p and std_p are the Mean and Standard Deviation of the resulting 'prior distribution' on each trial
  vector[N] mean_p;
  vector[N] std_p;
  
  for (n in 1:N){
    int s = subject[n];
    beta_n[n] = beta;
    rho_n[n] = rho;
    eta_n[n] = eta;
    mu_n[n] = mu;
    nu_n[n] = nu;
    
    //the 'prior distribution' results from the product of two Normal terms: A stimulus term and a que term
    mean_p[n] = get_mean3(X[n], q[n], mu_n[n], beta_n[n], rho_n[n] + eta_n[n]*std[n], nu_n[n]);
    std_p[n] = get_std3(X[n], q[n], mu_n[n], beta_n[n], rho_n[n] + eta_n[n]*std[n], nu_n[n]);
  }
}
model {
  
  //Weak priors on the parameters make the simulations more stable
  //beta ~ cauchy(0,100);
  //rho ~ cauchy(0,100);
  //eta ~ cauchy(0,100);
  //mu ~ cauchy(0,100);
  //nu ~ cauchy(0,100);
  
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

