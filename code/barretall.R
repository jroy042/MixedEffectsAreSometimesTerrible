#### replicate Barr et al (2013)

compare_lmer_stan_barretal <- function(B=10000,threads=8){
	# B is the number of models to run for each of the 8 model types
	# threads is the number of threads to split this task by; set to 1 to not parallelize
	#	in which case the individual stan models will be parallelized
	
	loop_lmer_stan_barretal <- function(pardat){
		require(afex)
		require(rstan)
		
		{### create data frames without response

		## within-item, 24 items
		S <- 24 # number of subjects
		I <- 24 # number of items

		wi24 <- expand.grid(subject=1:S,item=1:I)
		wi24$x <- 0.5
		wi24[wi24$subject %in% 1:(S/2) & wi24$item %in% 1:(I/2),"x"] <- -0.5
		wi24[wi24$subject %in% (S/2+1):S & wi24$item %in% (I/2+1):I,"x"] <- -0.5

		## within-item, 12 items
		S <- 24 # number of subjects
		I <- 12 # number of items

		wi12 <- expand.grid(subject=1:S,item=1:I)
		wi12$x <- 0.5
		wi12[wi12$subject %in% 1:(S/2) & wi12$item %in% 1:(I/2),"x"] <- -0.5
		wi12[wi12$subject %in% (S/2+1):S & wi12$item %in% (I/2+1):I,"x"] <- -0.5

		## between-item, 24 items
		S <- 24 # number of subjects
		I <- 24 # number of items

		bi24 <- expand.grid(subject=1:S,item=1:I)
		bi24$x <- 0.5
		bi24[bi24$item %in% 1:(I/2),"x"] <- -0.5

		## between-item, 12 items
		S <- 24 # number of subjects
		I <- 12 # number of items

		bi12 <- expand.grid(subject=1:S,item=1:I)
		bi12$x <- 0.5
		bi12[bi12$item %in% 1:(I/2),"x"] <- -0.5
		}

		### Generate data frame with y and missing data
		generate_linear_data <- function(btwn,I,H0){
			# I should be 12 or 24
			# btwn should be TRUE for between-item, FALSE for within-item
			# H0 should be TRUE if the treatment effect is 0 and FALSE if it is 0.8
			
			require(MASS)
			
			S <- 24
			
			if(btwn){
				if(I==12){
					dat <- bi12
				} else dat <- bi24
			} else if(I==12){
				dat <- wi12
			} else dat <- wi24
			N <- nrow(dat)
			
			# Barr et al (2013) Table 2
			b0 <- runif(1,-3,3)
			b1 <- ifelse(H0,0,0.8)
			sigma_subj <- sqrt(runif(2,0,3))
			rho_subj <- runif(1,-0.8,0.8)
			if(!btwn){
				sigma_item <- sqrt(runif(2,0,3))
				rho_item <- runif(1,-0.8,0.8)
			} else {
				sigma_item <- sqrt(runif(1,0,3))
				rho_item <- NA
			}
			sigma_res <- sqrt(runif(1,0,3))
			p_missing <- runif(1,0,0.05)
			
			gamma_subj <- mvrnorm(n = S, mu = c(0,0),
				Sigma = diag(sigma_subj) %*% matrix(c(1,rho_subj,rho_subj,1),2,2) %*% diag(sigma_subj))
			
			if(btwn){
				gamma_item <- cbind(rnorm(I,0,sigma_item),0)
			} else {
				gamma_item <- mvrnorm(n = I, mu = c(0,0),
					Sigma = diag(sigma_item) %*% matrix(c(1,rho_item,rho_item,1),2,2) %*% diag(sigma_item))
			}
			
			dat$y <- 0
			for(n in 1:N){
				dat[n,"y"] <- sum(c(1,dat[n,"x"]) * 
					(c(b0,b1) + gamma_subj[dat[n,"subject"],] + gamma_item[dat[n,"item"],]))
			}
			dat$y <- dat$y + rnorm(N,0,sigma_res)
			dat <- dat[-sample(1:N,ceiling(N*p_missing),replace=FALSE),]
			
			return(list(dat = dat, pars = data.frame(b0, b1, s0=sigma_subj[1], s1=sigma_subj[2], rho_subj,
				i0 = sigma_item[1], i1 = ifelse(btwn,NA,sigma_item[2]), rho_item, sigma_res, p_missing)))
		}
		
		### Fit single model with both lmer and stan
		model_lmer_stan_barretal <- function(btwn,I,H0){
			dat <- generate_linear_data(btwn,I,H0)
			data <- dat$dat
			dat <- dat$pars
			y_mu <- mean(data$y)
			y_sd <- sd(data$y)
			stan_data <- list(N = nrow(data), S = max(data$subject), I = max(data$item),
				y = data$y, y_mu = y_mu, y_sd = y_sd,
				x = data$x, subj = data$subject, item = data$item)

			stan_pars <- c("b0","b1","sigma_e","sigma_subj","sigma_item",
				"gamma_subj","gamma_item","rho_subj","y_hat")
				
			if(btwn){
				formula <- y ~ x + (1 + x | subject) + (1 | item)
				stan_model <- stan_model_bi
			} else {
				formula <- y ~ x + (1 + x | subject) + (1 + x | item)
				stan_model <- stan_model_wi
				stan_pars <- c(stan_pars,"rho_item")
			}
			
			data$y <- scale(data$y)[,1]
			
			m_lmer <- mixed(formula, data)
			p_data_H0_lmer <- m_lmer$anova$P
			conv_rest <- summary(m_lmer$restricted$x)$optinfo$conv$lme4
			m_lmer <- summary(m_lmer$full)
			conv_lmer <- m_lmer$optinfo$conv$lme4
			t_lmer <- m_lmer$coef[2,3]
			if(!is.null(names(conv_rest))){
				conv_rest <- paste(conv_rest$messages,collapse="; ")
				conv_rest <- str_replace_all(conv_rest,",","")
				conv_rest <- str_replace_all(conv_rest,"\n","")
			} else conv_rest <- ""
			if(!is.null(names(conv_lmer))){
				conv_lmer <- paste(conv_lmer$messages,collapse="; ")
				conv_lmer <- str_replace_all(conv_lmer,",","")
				conv_lmer <- str_replace_all(conv_lmer,"\n","")
			} else conv_lmer <- ""
			
			m_stan <- sampling(object = stan_model, data = stan_data, pars = stan_pars,
				chains = 3, iter = 13000, warmup = 3000, thin = 10)
			samp <- data.frame(do.call(rbind,args=get_sampler_params(m_stan,inc_warmup=FALSE)))
			nmtd <- sum(samp$treedepth__>10)
			ndiv <- sum(samp$divergent__)
			p_HA_data_stan <- mean(extract(m_stan,pars="b1")$b1>0)
			m_stan <- summary(m_stan)$summary
			ci_0_stan <- m_stan["b1",4]<0 & m_stan["b1",8]>0
			rhat <- quantile(m_stan[,"Rhat"],c(0,.25,.5,.75,1))
			neff <- quantile(m_stan[,"n_eff"],c(0,.25,.5,.75,1))
			m_stan <- m_stan[,1]
			
			vc <- data.frame(m_lmer$varcor)
			if(btwn){
				i1_lmer <- ri_lmer <- i1_stan <- ri_stan <- NA
				i0_stan <- m_stan["sigma_item"]
				res_lmer <- vc[5,5] * y_sd
			} else {
				i1_lmer <- vc[5,5] * y_sd
				ri_lmer <- vc[6,5]
				res_lmer <- vc[7,5] * y_sd
				i0_stan <- m_stan["sigma_item[1]"]
				i1_stan <- m_stan["sigma_item[2]"]
				ri_stan <- m_stan["rho_item"]
			}
			
			return(data.frame(
				btwn, I, H0, pmiss = dat$p_missing,
				b0 = dat$b0, b0_lmer = m_lmer$coef[1,1]*y_sd+y_mu, b0_stan = m_stan["b0"],
				b1 = dat$b1, b1_lmer = m_lmer$coef[2,1]*y_sd, b1_stan = m_stan["b1"],
				s0 = dat$s0, s0_lmer = vc[1,5]*y_sd, s0_stan = m_stan["sigma_subj[1]"],
				s1 = dat$s1, s1_lmer = vc[2,5]*y_sd, s1_stan = m_stan["sigma_subj[2]"],
				rs = dat$rho_subj, rs_lmer = vc[3,5], rs_stan = m_stan["rho_subj"],
				i0 = dat$i0, i0_lmer = vc[4,5]*y_sd, i0_stan,
				i1 = dat$i1, i1_lmer, i1_stan,
				ri = dat$rho_item, ri_lmer, ri_stan,
				res = dat$sigma_res, res_lmer, res_stan = m_stan["sigma_e"],
				t_lmer, p_data_H0_lmer, p_HA_data_stan, ci_0_stan,
				conv_lmer, conv_rest, nmtd, ndiv,
				rhat_min = rhat[1], rhat_q25 = rhat[2], rhat_med = rhat[3], rhat_q75 = rhat[4], rhat_max = rhat[5],
				neff_min = neff[1], neff_q25 = neff[2], neff_med = neff[3], neff_q75 = neff[4], neff_max = neff[5]))
		}

		## unpack arguments
		B <- pardat[[1]]
		thread <- pardat[[2]]
		stan_model_wi <- pardat[[3]]
		stan_model_bi <- pardat[[4]]

		### run B iterations of each of the 8 models and output progress to thread-specific file
		barretal_comparison <- NULL
		tfile <- paste(getwd(),"/","barretal_thread_",thread,".progress",sep="")
		tcsvfile <- paste(getwd(),"/","barretal_thread_",thread,".csv",sep="")
		cat(paste(Sys.time(),"Starting Thread",thread),file=tfile,sep="\n",append=FALSE)
		for(b in 1:B){  # Each iteration has 8 calls to model_lmer_stan_barretal, one for each combination
			cat("\n",paste(Sys.time(),"Iteration",b,"/",B),file=tfile,sep="",append=TRUE)
			barretal_comparison <- rbind(barretal_comparison,model_lmer_stan_barretal(TRUE,12,TRUE))
			cat(" . ",file=tfile,sep="",append=TRUE)
			barretal_comparison <- rbind(barretal_comparison,model_lmer_stan_barretal(TRUE,12,FALSE))
			cat(" . ",file=tfile,sep="",append=TRUE)
			barretal_comparison <- rbind(barretal_comparison,model_lmer_stan_barretal(FALSE,12,TRUE))
			cat(" . ",file=tfile,sep="",append=TRUE)
			barretal_comparison <- rbind(barretal_comparison,model_lmer_stan_barretal(FALSE,12,FALSE))
			cat(" . ",file=tfile,sep="",append=TRUE)
			barretal_comparison <- rbind(barretal_comparison,model_lmer_stan_barretal(TRUE,24,TRUE))
			cat(" . ",file=tfile,sep="",append=TRUE)
			barretal_comparison <- rbind(barretal_comparison,model_lmer_stan_barretal(TRUE,24,FALSE))
			cat(" . ",file=tfile,sep="",append=TRUE)
			barretal_comparison <- rbind(barretal_comparison,model_lmer_stan_barretal(FALSE,24,TRUE))
			cat(" . ",file=tfile,sep="",append=TRUE)
			barretal_comparison <- rbind(barretal_comparison,model_lmer_stan_barretal(FALSE,24,FALSE))
			cat(" . ",file=tfile,sep="",append=TRUE)
			write.csv(barretal_comparison,file=tcsvfile)
		}
		cat("\n",paste(Sys.time(),"Finished Thread",thread),file=tfile,sep="",append=TRUE)
		
		### return results
		return(barretal_comparison)
	}

	{### Stan code

{## code for within-item
stan_code_wi <- "// stan code for within-item Barr et al (2013) replication
data {
  int<lower=2> N;
  int<lower=2> S;
  int<lower=2> I;

  vector[N] y;
  real y_mu;
  real<lower=0> y_sd;

  vector[N] x;
  int<lower=1,upper=S> subj[N];
  int<lower=1,upper=I> item[N];
}

parameters {
  real b0_raw;
  real b1_raw;
  real<lower=0> sigma_e_raw;

  vector<lower=0>[2] sigma_subj_raw;
  cholesky_factor_corr[2] L_subj;
  vector[2] gamma_subj_raw[S];

  vector<lower=0>[2] sigma_item_raw;
  cholesky_factor_corr[2] L_item;
  vector[2] gamma_item_raw[I];
}

transformed parameters {
  real b0;
  real b1;
  real<lower=0> sigma_e;

  vector<lower=0>[2] sigma_subj;
  vector[2] gamma_subj[S];

  vector<lower=0>[2] sigma_item;
  vector[2] gamma_item[I];

  vector[N] y_hat;

  b0 = y_mu + 0.25 * y_sd * b0_raw;
  b1 = y_sd * b1_raw;
  sigma_e = 0.5 * y_sd * sigma_e_raw;

  sigma_subj = y_sd * sigma_subj_raw;
  for(s in 1:S){
    gamma_subj[s] = sigma_subj .* (L_subj * gamma_subj_raw[s]);
  }

  sigma_item = y_sd * sigma_item_raw;
  for(i in 1:I){
    gamma_item[i] = sigma_item .* (L_item * gamma_item_raw[i]);
  }

  for(n in 1:N){
    y_hat[n] = b0 + gamma_subj[subj[n],1] + gamma_item[item[n],1]
        + x[n] * (b1 + gamma_subj[subj[n],2] + gamma_item[item[n],2]);
  }
}

model {
  b0_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  sigma_e_raw ~ normal(0,1);

  sigma_subj_raw ~ normal(0,1);
  L_subj ~ lkj_corr_cholesky(2);
  for(s in 1:S) gamma_subj_raw[s] ~ normal(0,1);

  sigma_item_raw ~ normal(0,1);
  L_item ~ lkj_corr_cholesky(2);
  for(i in 1:I) gamma_item_raw[i] ~ normal(0,1);

  y ~ normal(y_hat,sigma_e);
}

generated quantities {
  real<lower=-1,upper=1> rho_subj;
  real<lower=-1,upper=1> rho_item;
  {
    matrix[2,2] omega_subj;
    matrix[2,2] omega_item;

    omega_subj = tcrossprod(L_subj);
    rho_subj = omega_subj[2,1];
    omega_item = tcrossprod(L_item);
    rho_item = omega_item[2,1];
  }
}
"
}

{## code for between-item
stan_code_bi <- "// stan code for between-item Barr et al (2013) replication
data {
  int<lower=2> N;
  int<lower=2> S;
  int<lower=2> I;

  vector[N] y;
  real y_mu;
  real<lower=0> y_sd;

  vector[N] x;
  int<lower=1,upper=S> subj[N];
  int<lower=1,upper=I> item[N];
}

parameters {
  real b0_raw;
  real b1_raw;
  real<lower=0> sigma_e_raw;

  vector<lower=0>[2] sigma_subj_raw;
  cholesky_factor_corr[2] L_subj;
  vector[2] gamma_subj_raw[S];

  real<lower=0> sigma_item_raw;
  vector[I] gamma_item_raw;
}

transformed parameters {
  real b0;
  real b1;
  real<lower=0> sigma_e;

  vector<lower=0>[2] sigma_subj;
  vector[2] gamma_subj[S];

  real<lower=0> sigma_item;
  vector[I] gamma_item;

  vector[N] y_hat;

  b0 = y_mu + 0.25 * y_sd * b0_raw;
  b1 = y_sd * b1_raw;
  sigma_e = 0.5 * y_sd * sigma_e_raw;

  sigma_subj = y_sd * sigma_subj_raw;
  for(s in 1:S){
	gamma_subj[s] = sigma_subj .* (L_subj * gamma_subj_raw[s]);
  }

  sigma_item = y_sd * sigma_item_raw;
  gamma_item = sigma_item * gamma_item_raw;

  for(n in 1:N){
    y_hat[n] = b0 + gamma_subj[subj[n],1] + gamma_item[item[n]]
        + x[n] * (b1 + gamma_subj[subj[n],2]);
  }
}

model {
  b0_raw ~ normal(0,1);
  b1_raw ~ normal(0,1);
  sigma_e_raw ~ normal(0,1);

  sigma_subj_raw ~ normal(0,1);
  L_subj ~ lkj_corr_cholesky(2);
  for(s in 1:S) gamma_subj_raw[s] ~ normal(0,1);

  sigma_item_raw ~ normal(0,1);
  gamma_item_raw ~ normal(0,1);

  y ~ normal(y_hat,sigma_e);
}

generated quantities {
  real<lower=-1,upper=1> rho_subj;
  {
    matrix[2,2] omega_subj;

    omega_subj = tcrossprod(L_subj);
    rho_subj = omega_subj[2,1];
  }
}
"
}
	}
	
	require(rstan)
	stan_model_wi <- stan_model(model_name = "barretal_within_item", model_code = stan_code_wi, save_dso = TRUE)
	stan_model_bi <- stan_model(model_name = "barretal_between_item", model_code = stan_code_bi, save_dso = TRUE)

	if(threads==1){
		options(mc.cores = parallel::detectCores())
		barretal_comparison <- loop_lmer_stan_barretal(list(B,threads,stan_model_wi,stan_model_bi))
	} else {
		require(parallel)
		cl = makeCluster(rep("localhost", threads))
		
		# simulations per thread and thread number
		left <- B %% threads
		B <- rep(floor(B/threads),threads)
		if(left>0) B[1:left] <- B[1:left] + 1
		pardat <- list()
		for(i in 1:threads) pardat[[i]] <- list(B[i],i,stan_model_wi,stan_model_bi)
		
		barretal_comparison <- parLapply(cl,pardat,loop_lmer_stan_barretal)
		barretal_comparison <- do.call(rbind,args=barretal_comparison)
	}
	
	return(barretal_comparison)
}
