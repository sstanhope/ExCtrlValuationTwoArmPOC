

# Two sample economic model, comparison of means - Analytic approach

ext_control_means_analytic <- function(d_tgt,				# Targeted difference in means, expressed in terms of (unit) standard deviation
		       		       d_true,				# True difference in means, expressed in terms of (unit) standard deviation
				       t_obs,				# Observation time for trial participants (in time units)
                       		       alpha_poc,			# Type I error for POC analysis
		       		       power_poc,			# Power for POC analysis
		       		       n_ext,				# Number of external controls for POC analysis
		       		       alpha_conf,			# Type I error for confirmatory analysis
		       		       power_conf,			# Power for confirmatory analysis
		       		       cost_subject,			# Cost of a trial subject
		       		       enrollment_rate,			# Enrolled patients per unit time
				       compound_value,			# Present value of the successful compound today
		       		       discount_rate			# Discount rate per unit time
		       		       )
{

# POC design - Sample size per arm for POC analysis

n_poc <- ceiling(power.t.test(delta=d_tgt,sd=1,sig.level=alpha_poc,power=power_poc)$n)	
cost_poc <- (2*n_poc*cost_subject) / (1+discount_rate)^( ceiling(2*n_poc/enrollment_rate) + t_obs)

# True power of POC test

z_alpha_poc <- qnorm(1 - alpha_poc/2)
pwr_poc <- 1 - pnorm(z_alpha_poc,mean=d_true/sqrt(1/n_poc + 1/(n_poc+n_ext)))

# Expected delta conditional on POC success

f.cond.norm <- integrate(dnorm,lower=z_alpha_poc*sqrt(1/n_poc + 1/(n_poc+n_ext)),upper=Inf,mean=d_true,sd=sqrt(1/n_poc + 1/(n_poc+n_ext)))$value

f.cond <- function(x,mean=mean,sd=sd)  {
        ret <- x * dnorm(x,mean=mean,sd=sd)
        ret
        }
d_condpoc <- (1/f.cond.norm) * integrate(f.cond,lower=z_alpha_poc*sqrt(1/n_poc + 1/(n_poc+n_ext)),upper=Inf,mean=d_true,sd=sqrt(1/n_poc + 1/(n_poc+n_ext)))$value

# Confirmatory design

n_conf <- ceiling(power.t.test(delta=d_condpoc,sd=1,sig.level=alpha_conf,power=power_conf)$n)	# Sample size per arm for confirmatory analysis
cost_conf <- (2*n_conf*cost_subject) / (1+discount_rate)^( ceiling(2*n_poc/enrollment_rate) + t_obs + ceiling(2*n_conf/enrollment_rate) + t_obs )

# Power of confirmatory test

z_alpha_conf <- qnorm(1 - alpha_conf/2)
pwr_conf <- 1 - pnorm(z_alpha_conf,mean=d_true/sqrt(1/n_conf + 1/n_conf))

# Experimental outputs

n_total <- 2*n_poc + pwr_poc*(2*n_conf)
pwr_total <- pwr_poc*pwr_conf
cost_total <- cost_poc + pwr_poc*cost_conf
value_total <- -cost_poc + -cost_conf*pwr_poc + (compound_value / (1+discount_rate)^( ceiling(2*n_poc/enrollment_rate) + t_obs + ceiling(2*n_conf/enrollment_rate) + t_obs ))*pwr_total

# Conclude

dat.out <- data.frame(n_poc=n_poc,cost_poc=cost_poc,pwr_poc=pwr_poc,d_condpoc=d_condpoc,
		      n_conf=n_conf,cost_conf=cost_conf,pwr_conf=pwr_conf,
		      n_total=n_total,pwr_total=pwr_total,cost_total=cost_total,value_total=value_total)
dat.out

}



# Two sample economic model, comparison of proportions - Analytic approach

ext_control_proportions_analytic <- function(p_trt,				# Targeted proportion of responders in treatment group
		       		       	     p_trt_true,			# True proportion of responders in treatment group
				             p_con,				# Targeted proportion of responders in treatment group
		       		             p_con_true,			# True proportion of responders in treatment group
					     t_obs,				# Observation time for trial participants (in time units)
                       		             alpha_poc,			        # Type I error for POC analysis
		       		             power_poc,			        # Power for POC analysis
		       		             n_ext,				# Number of external controls for POC analysis
		       		             alpha_conf,			# Type I error for confirmatory analysis
		       		             power_conf,			# Power for confirmatory analysis
		       		             cost_subject,			# Cost of a trial subject
		       		             enrollment_rate,			# Enrolled patients per unit time
				             compound_value,			# Present value of the successful compound today
		       		             discount_rate			# Discount rate per unit time
		       		             )
{

# POC design - Sample size per arm for POC analysis

z_alpha_poc <- qnorm(1 - alpha_poc/2)
z_beta_poc <- qnorm(power_poc)

n_poc <- ceiling( ((z_alpha_poc+z_beta_poc)^2)*(p_trt*(1-p_trt) + p_con*(1-p_con)) / (p_trt-p_con)^2 )

cost_poc <- (2*n_poc*cost_subject) / (1+discount_rate)^( ceiling(2*n_poc/enrollment_rate) + t_obs )

# True power of POC test

pwr_poc <- 1 - pnorm(z_alpha_poc,mean=(p_trt_true - p_con_true)/sqrt(p_trt_true*(1-p_trt_true)/n_poc + p_con_true*(1-p_con_true)/(n_poc+n_ext)))

# Expected proportions conditional on POC success

trt.mc <- rnorm(100000,mean=p_trt_true,sd=sqrt(p_trt_true*(1-p_trt_true)/n_poc))
con.mc <- rnorm(100000,mean=p_con_true,sd=sqrt(p_con_true*(1-p_con_true)/(n_poc+n_ext)))
ind <- trt.mc<0 | con.mc<0 | trt.mc>1 | con.mc>1
trt.mc <- trt.mc[!ind]
con.mc <- con.mc[!ind]
ind <- ( (trt.mc - con.mc) / sqrt(trt.mc*(1-trt.mc)/n_poc + con.mc*(1-con.mc)/(n_poc+n_ext)) ) > z_alpha_poc
p_trt_condpoc <- mean(trt.mc[ind])
p_con_condpoc <- mean(con.mc[ind])

# Confirmatory design

z_alpha_conf <- qnorm(1 - alpha_conf/2)
z_beta_conf <- qnorm(power_conf)

n_conf <- ceiling( ((z_alpha_conf+z_beta_conf)^2)*(p_trt_condpoc*(1-p_trt_condpoc) + p_con_condpoc*(1-p_con_condpoc)) / (p_trt_condpoc-p_con_condpoc)^2 )

cost_conf <- (2*n_conf*cost_subject) / (1+discount_rate)^( ceiling(2*n_poc/enrollment_rate) + t_obs + ceiling(2*n_conf/enrollment_rate) + t_obs )

# Power of confirmatory test

pwr_conf <- 1 - pnorm(z_alpha_conf,mean=(p_trt_true - p_con_true)/sqrt(p_trt_true*(1-p_trt_true)/n_conf + p_con_true*(1-p_con_true)/n_conf))

# Experimental outputs

n_total <- 2*n_poc + pwr_poc*(2*n_conf)
pwr_total <- pwr_poc*pwr_conf
cost_total <- cost_poc + pwr_poc*cost_conf
value_total <- -cost_poc + -cost_conf*pwr_poc + (compound_value / (1+discount_rate)^( ceiling(2*n_poc/enrollment_rate) + t_obs + ceiling(2*n_conf/enrollment_rate) + t_obs ))*pwr_total

# Conclude

dat.out <- data.frame(n_poc=n_poc,cost_poc=cost_poc,pwr_poc=pwr_poc,
		      p_trt_condpoc=p_trt_condpoc,p_con_condpoc=p_con_condpoc,
		      n_conf=n_conf,cost_conf=cost_conf,pwr_conf=pwr_conf,
		      n_total=n_total,pwr_total=pwr_total,cost_total=cost_total,value_total=value_total)
dat.out

}



# Two sample economic model, comparison of survival and hazard ratio - Analytic approach

ext_control_surv_hazard_analytic <- function(s_trt,				# Targeted proportion of survivors at time t_obs in treatment group
		       		       	     s_trt_true,			# True proportion of survivors at time t_obs in treatment group
				             s_con,				# Targeted proportion of survivors at time t_obs in control group
		       		             s_con_true,			# True proportion of survivors at time t_obs in control group
					     t_obs,				# Observation time for trial participants
                       		             alpha_poc,			        # Type I error for POC analysis
		       		             power_poc,			        # Power for POC analysis
		       		             n_ext,				# Number of external controls for POC analysis
		       		             alpha_conf,			# Type I error for confirmatory analysis
		       		             power_conf,			# Power for confirmatory analysis
		       		             cost_subject,			# Cost of a trial subject
		       		             enrollment_rate,			# Enrolled patients per unit time
				             compound_value,			# Present value of the successful compound today
		       		             discount_rate			# Discount rate per unit time
		       		             )
{

# POC design - Sample size per arm for POC analysis

z_alpha_poc <- qnorm(1 - alpha_poc/2)
z_beta_poc <- qnorm(power_poc)

irr_poc <- log(s_trt) / log(s_con)
n_poc <- ceiling( ((irr_poc + 1) / (irr_poc - 1))^2 * (z_alpha_poc + z_beta_poc)^2 / (1 - s_trt + 1 - s_con) )

cost_poc <- (2*n_poc*cost_subject) / (1+discount_rate)^( (2*n_poc/enrollment_rate) + t_obs )

# True power of POC test

irr_true <- log(s_trt_true) / log(s_con_true)
pwr_poc <- pnorm( sqrt( n_poc*(1-s_trt_true) + (n_poc + n_ext)*(1-s_con_true) )*abs(irr_true - 1)/(irr_true + 1) - z_alpha_poc)

# Expected proportions conditional on POC success - note: calculation assumes no censoring

trt.mc <- rnorm(100000,mean=s_trt_true,sd=sqrt(s_trt_true*(1-s_trt_true)/n_poc))
con.mc <- rnorm(100000,mean=s_con_true,sd=sqrt(s_con_true*(1-s_con_true)/(n_poc+n_ext)))
ind <- trt.mc<0 | con.mc<0 | trt.mc>1 | con.mc>1
trt.mc <- trt.mc[!ind]
con.mc <- con.mc[!ind]
irr.mc <- log(trt.mc) / log(con.mc)
ind <- abs( sqrt( n_poc*(1-trt.mc) + (n_poc + n_ext)*(1-con.mc) )*abs(irr.mc - 1)/(irr.mc + 1) ) > z_alpha_poc
s_trt_condpoc <- mean(trt.mc[ind])
s_con_condpoc <- mean(con.mc[ind])

# Confirmatory design

z_alpha_conf <- qnorm(1 - alpha_conf/2)
z_beta_conf <- qnorm(power_conf)

irr_conf <- log(s_trt_condpoc) / log(s_con_condpoc)
n_conf <- ceiling( ((irr_conf + 1) / (irr_conf - 1))^2 * (z_alpha_conf + z_beta_conf)^2 / (1 - s_trt_condpoc + 1 - s_con_condpoc) )

cost_conf <- (2*n_conf*cost_subject) / (1+discount_rate)^( ceiling(2*n_conf/enrollment_rate) + t_obs + ceiling(2*n_conf/enrollment_rate) + t_obs )

# Power of confirmatory test

pwr_conf <- pnorm( sqrt( n_conf*(1-s_trt_true) + n_conf*(1-s_con_true) )*abs(irr_true - 1)/(irr_true + 1) - z_alpha_conf)

# Experimental outputs

n_total <- 2*n_poc + pwr_poc*(2*n_conf)
pwr_total <- pwr_poc*pwr_conf
cost_total <- cost_poc + pwr_poc*cost_conf
value_total <- -cost_poc + -cost_conf*pwr_poc + (compound_value / (1+discount_rate)^( ceiling(2*n_poc/enrollment_rate) + t_obs + ceiling(2*n_conf/enrollment_rate) + t_obs ))*pwr_total

# Conclude

dat.out <- data.frame(n_poc=n_poc,cost_poc=cost_poc,pwr_poc=pwr_poc,
		      s_trt_condpoc=s_trt_condpoc,s_con_condpoc=s_con_condpoc,
		      n_conf=n_conf,cost_conf=cost_conf,pwr_conf=pwr_conf,
		      n_total=n_total,pwr_total=pwr_total,cost_total=cost_total,value_total=value_total)
dat.out

}

