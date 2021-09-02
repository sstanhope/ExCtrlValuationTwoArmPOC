rm(list=ls())

options(width=200)

source('aca_valuation_library_20210428.R')


c.v.0 <- 1070000000
c.v.1 <- 0

ta.mult <- 1


# Two sample comparison of proportions analysis - fixed example calibrated to belimumab

tspa.0 <- ext_control_proportions_analytic(p_trt=0.57,p_trt_true=0.57,p_con=0.35,p_con_true=0.35,t_obs=12,
                  	             	   alpha_poc=0.05,power_poc=0.80,n_ext=0,
                  	             	   alpha_conf=0.05,power_conf=0.90,
                  	             	   cost_subject=42000,enrollment_rate=10,compound_value=c.v.0,
                  	             	   discount_rate=0.11/12)
print(tspa.0)

tspa.ext <- ext_control_proportions_analytic(p_trt=0.57,p_trt_true=0.57,p_con=0.35,p_con_true=0.35,t_obs=12,
                  	             	   alpha_poc=0.05,power_poc=0.80,n_ext=ta.mult*tspa.0$n_poc,
                  	             	   alpha_conf=0.05,power_conf=0.90,
                  	             	   cost_subject=42000,enrollment_rate=10,compound_value=c.v.0,
                  	             	   discount_rate=0.11/12)
print(tspa.ext)

tspa.ext.match <- ext_control_proportions_analytic(p_trt=0.57,p_trt_true=0.57,p_con=0.35,p_con_true=0.35,t_obs=12,
                  	             	           alpha_poc=0.05,power_poc=tspa.ext$pwr_poc,n_ext=0,
                  	             	           alpha_conf=0.05,power_conf=0.90,
                  	             	           cost_subject=42000,enrollment_rate=10,compound_value=c.v.0,
                  	             	           discount_rate=0.11/12)
print(tspa.ext.match)

print('TSPA - Baseline, Ext, Ext Match, Ext vs Ext Match')
print(paste('Total Power:',tspa.0$pwr_total,tspa.ext$pwr_total,tspa.ext.match$pwr_total,tspa.ext$pwr_total - tspa.ext.match$pwr_total))
print(paste('Total Cost:',tspa.0$cost_total,tspa.ext$cost_total,tspa.ext.match$cost_total,tspa.ext$cost_total - tspa.ext.match$cost_total))
print(paste('Total Value:',tspa.0$value_total,tspa.ext$value_total,tspa.ext.match$value_total,tspa.ext$value_total - tspa.ext.match$value_total))


# Two sample comparison of proportions analysis - simulation example

mc.its <- 1000

n_poc.0 <- rep(NA,mc.its)
cost_poc.0 <- rep(NA,mc.its)
pwr_poc.0 <- rep(NA,mc.its)
p_trt_condpoc.0 <- rep(NA,mc.its)
p_con_condpoc.0 <- rep(NA,mc.its)
n_conf.0 <- rep(NA,mc.its)
cost_conf.0 <- rep(NA,mc.its)
pwr_conf.0 <- rep(NA,mc.its)
n_total.0 <- rep(NA,mc.its)
pwr_total.0 <- rep(NA,mc.its)
cost_total.0 <- rep(NA,mc.its)
value_total.0 <- rep(NA,mc.its)

n_poc.ext <- rep(NA,mc.its)
cost_poc.ext <- rep(NA,mc.its)
pwr_poc.ext <- rep(NA,mc.its)
p_trt_condpoc.ext <- rep(NA,mc.its)
p_con_condpoc.ext <- rep(NA,mc.its)
n_conf.ext <- rep(NA,mc.its)
cost_conf.ext <- rep(NA,mc.its)
pwr_conf.ext <- rep(NA,mc.its)
n_total.ext <- rep(NA,mc.its)
pwr_total.ext <- rep(NA,mc.its)
cost_total.ext <- rep(NA,mc.its)
value_total.ext <- rep(NA,mc.its)

n_poc.ext.match <- rep(NA,mc.its)
cost_poc.ext.match <- rep(NA,mc.its)
pwr_poc.ext.match <- rep(NA,mc.its)
p_trt_condpoc.ext.match <- rep(NA,mc.its)
p_con_condpoc.ext.match <- rep(NA,mc.its)
n_conf.ext.match <- rep(NA,mc.its)
cost_conf.ext.match <- rep(NA,mc.its)
pwr_conf.ext.match <- rep(NA,mc.its)
n_total.ext.match <- rep(NA,mc.its)
pwr_total.ext.match <- rep(NA,mc.its)
cost_total.ext.match <- rep(NA,mc.its)
value_total.ext.match <- rep(NA,mc.its)

for (i in 1:mc.its)  {
	d_true <- sample(c(0.35,0.46,0.57),size=1,prob=c(0.5,0.35,0.15))

	tspa.mc.0 <- ext_control_proportions_analytic(p_trt=0.57,p_trt_true=d_true,p_con=0.35,p_con_true=0.35,t_obs=12,
	                  	             	      alpha_poc=0.05,power_poc=0.80,n_ext=0,
	                  	             	      alpha_conf=0.05,power_conf=0.90,
	                  	             	      cost_subject=42000,enrollment_rate=10,compound_value=c.v.0,
	                  	             	      discount_rate=0.11/12)

	tspa.mc.ext <- ext_control_proportions_analytic(p_trt=0.57,p_trt_true=d_true,p_con=0.35,p_con_true=0.35,t_obs=12,
	                  	             	     	alpha_poc=0.05,power_poc=0.80,n_ext=ta.mult*tspa.0$n_poc,
	                  	             	   	alpha_conf=0.05,power_conf=0.90,
	                  	             	   	cost_subject=42000,enrollment_rate=10,compound_value=c.v.0,
	                  	             	   	discount_rate=0.11/12)

	tspa.mc.ext.match <- ext_control_proportions_analytic(p_trt=0.57,p_trt_true=d_true,p_con=0.35,p_con_true=0.35,t_obs=12,
                  	             	           	      alpha_poc=0.05,power_poc=tspa.ext$pwr_poc,n_ext=0,
                  	             	           	      alpha_conf=0.05,power_conf=0.90,
                  	             	           	      cost_subject=42000,enrollment_rate=10,compound_value=c.v.0,
                  	             	           	      discount_rate=0.11/12)

	n_poc.0[i] <- tspa.mc.0$n_poc
	cost_poc.0[i] <- tspa.mc.0$cost_poc
	pwr_poc.0[i] <- tspa.mc.0$pwr_poc
	p_trt_condpoc.0[i] <- tspa.mc.0$p_trt_condpoc
	p_con_condpoc.0[i] <- tspa.mc.0$p_con_condpoc
	n_conf.0[i] <- tspa.mc.0$n_conf
	cost_conf.0[i] <- tspa.mc.0$cost_conf
	pwr_conf.0[i] <- tspa.mc.0$pwr_conf
	n_total.0[i] <- tspa.mc.0$n_total
	pwr_total.0[i] <- tspa.mc.0$pwr_total
	cost_total.0[i] <- tspa.mc.0$cost_total
	value_total.0[i] <- tspa.mc.0$value_total

	n_poc.ext[i] <- tspa.mc.ext$n_poc
	cost_poc.ext[i] <- tspa.mc.ext$cost_poc
	pwr_poc.ext[i] <- tspa.mc.ext$pwr_poc
	p_trt_condpoc.ext[i] <- tspa.mc.ext$p_trt_condpoc
	p_con_condpoc.ext[i] <- tspa.mc.ext$p_con_condpoc
	n_conf.ext[i] <- tspa.mc.ext$n_conf
	cost_conf.ext[i] <- tspa.mc.ext$cost_conf
	pwr_conf.ext[i] <- tspa.mc.ext$pwr_conf
	n_total.ext[i] <- tspa.mc.ext$n_total
	pwr_total.ext[i] <- tspa.mc.ext$pwr_total
	cost_total.ext[i] <- tspa.mc.ext$cost_total
	value_total.ext[i] <- tspa.mc.ext$value_total

	n_poc.ext.match[i] <- tspa.mc.ext.match$n_poc
	cost_poc.ext.match[i] <- tspa.mc.ext.match$cost_poc
	pwr_poc.ext.match[i] <- tspa.mc.ext.match$pwr_poc
	p_trt_condpoc.ext.match[i] <- tspa.mc.ext.match$p_trt_condpoc
	p_con_condpoc.ext.match[i] <- tspa.mc.ext.match$p_con_condpoc
	n_conf.ext.match[i] <- tspa.mc.ext.match$n_conf
	cost_conf.ext.match[i] <- tspa.mc.ext.match$cost_conf
	pwr_conf.ext.match[i] <- tspa.mc.ext.match$pwr_conf
	n_total.ext.match[i] <- tspa.mc.ext.match$n_total
	pwr_total.ext.match[i] <- tspa.mc.ext.match$pwr_total
	cost_total.ext.match[i] <- tspa.mc.ext.match$cost_total
	value_total.ext.match[i] <- tspa.mc.ext.match$value_total
	}

tspa.mc.out <- data.frame(n_poc.0,cost_poc.0,pwr_poc.0,p_trt_condpoc.0,p_con_condpoc.0,n_conf.0,cost_conf.0,pwr_conf.0,n_total.0,pwr_total.0,cost_total.0,value_total.0,
		          n_poc.ext,cost_poc.ext,pwr_poc.ext,p_trt_condpoc.ext,p_con_condpoc.ext,n_conf.ext,cost_conf.ext,pwr_conf.ext,n_total.ext,pwr_total.ext,cost_total.ext,value_total.ext,
		          n_poc.ext.match,cost_poc.ext.match,pwr_poc.ext.match,p_trt_condpoc.ext.match,p_con_condpoc.ext.match,n_conf.ext.match,cost_conf.ext.match,pwr_conf.ext.match,n_total.ext.match,pwr_total.ext.match,cost_total.ext.match,value_total.ext.match)

print(paste(mean(n_poc.0),mean(cost_poc.0),mean(pwr_poc.0),mean(p_trt_condpoc.0),mean(p_con_condpoc.0),mean(n_conf.0),mean(cost_conf.0),mean(pwr_conf.0),mean(n_total.0),mean(pwr_total.0),mean(cost_total.0),mean(value_total.0)))
print(paste(mean(n_poc.ext),mean(cost_poc.ext),mean(pwr_poc.ext),mean(p_trt_condpoc.ext),mean(p_con_condpoc.ext),mean(n_conf.ext),mean(cost_conf.ext),mean(pwr_conf.ext),mean(n_total.ext),mean(pwr_total.ext),mean(cost_total.ext),mean(value_total.ext)))
print(paste(mean(n_poc.ext.match),mean(cost_poc.ext.match),mean(pwr_poc.ext.match),mean(p_trt_condpoc.ext.match),mean(p_con_condpoc.ext.match),mean(n_conf.ext.match),mean(cost_conf.ext.match),mean(pwr_conf.ext.match),mean(n_total.ext.match),mean(pwr_total.ext.match),mean(cost_total.ext.match),mean(value_total.ext.match)))

print('TSPA MC Results - Baseline, Ext, Ext Match, Ext vs Ext Match')
print(paste('Total Power:',mean(tspa.mc.out$pwr_total.0),mean(tspa.mc.out$pwr_total.ext),mean(tspa.mc.out$pwr_total.ext.match),mean(tspa.mc.out$pwr_total.ext - tspa.mc.out$pwr_total.ext.match)))
print(paste('Total Cost:',mean(tspa.mc.out$cost_total.0),mean(tspa.mc.out$cost_total.ext),mean(tspa.mc.out$cost_total.ext.match),mean(tspa.mc.out$cost_total.ext - tspa.mc.out$cost_total.ext.match)))
print(paste('Total Value:',mean(tspa.mc.out$value_total.0),mean(tspa.mc.out$value_total.ext),mean(tspa.mc.out$value_total.ext.match),mean(tspa.mc.out$value_total.ext - tspa.mc.out$value_total.ext.match)))


