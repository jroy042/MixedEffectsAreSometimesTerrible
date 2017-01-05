#### compare_lme4_stan_imbalanced
#### imbalanced logistic and linear models

compare_lme4_stan_imbalanced <- function(B=10000,threads=16){
	# B is the number of models to run for each of the 2 model types
	# threads is the number of threads to split this task by; set to 1 to not parallelize
	# in which case the individual stan models will be parallelized

	single_thread <- function(pardat){
		require(lme4)
		require(rstan)
		require(stringr)
		require(gtools)
		require(clusterGeneration)
		options(contrasts=c("contr.sum","contr.poly"))

		### create data frame without response
		keep <- c("coef","s0","s1","s2","s3",
			"r01","r02","r03","r12","r13","r23")
		formula <- y ~ x1 + x2 + (1 + x1 + x2 | subject)

		### Generate model
		generate_model <- function(linear){
			balance <- function(dat){
				nobs <- nrow(dat)
				counts <- as.data.frame(xtabs(~.,dat))
				cells <- nrow(counts)
				non_zero_counts <- subset(counts, Freq > 0)
				non_empty_cells <- nrow(non_zero_counts)
				pct_non_empty_cells <- non_empty_cells / cells
				expected <- nobs / non_empty_cells
				return(pct_non_empty_cells * (
					1 - sum((non_zero_counts$Freq - expected)^2)/nobs^2))
			}
			balance_ratio <- function(dat){
				nobs <- nrow(dat)
				counts <- as.data.frame(xtabs(~.,dat))
				cells <- nrow(counts)
				non_zero_counts <- subset(counts, Freq > 0)
				non_empty_cells <- nrow(non_zero_counts)
				pct_non_empty_cells <- non_empty_cells / cells
				expected <- nobs / non_empty_cells
				freq <- non_zero_counts$Freq
				ratio <- freq / expected
				ratio[ratio < 1] <- 1 / ratio[ratio < 1]
				bal <- mean(ratio) / pct_non_empty_cells
				bal <- 2 * (1 - bal / (1 + bal))
				return(bal)
			}

			mod <- list(linear = linear)

			# S is the total number of subjects
			# L is the mean observations per subject
			# NS is a vector of the number of observations for each subject
			mod$S <- S <- sample(30:60,1)
			mod$L <- L <- runif(1,20,30)
			NS <- rpois(S,L)

			# Create data frame and check balance
			subject <- NULL
			for(s in 1:S) subject <- c(subject,rep(s,NS[s]))
			x1 = sample(c("a","b"), sum(NS),
				prob=rdirichlet(1,c(1,1)), replace=T)
			x2 = sample(c("a","b","c"), sum(NS),
				prob=rdirichlet(1,c(1,1,1)), replace=T)
			x1[sample(1:sum(NS),2,replace=F)] <- c("a","b")
			x2[sample(1:sum(NS),3,replace=F)] <- c("a","b","c")
			mod$frame <- data.frame(x1,x2,subject)
			mod$N <- N <- nrow(mod$frame)
			mod$balance <- balance(mod$frame)
			mod$balance_ratio <- balance_ratio(mod$frame)
			mod$frame$y <- 0
			
			# Create model matrix
			x <- model.matrix(nobars(formula),mod$frame)
			z <- mkReTrms(findbars(formula),mod$frame)
			mod$mat <- cbind(x,t(z$Ztlist[[1]]))
			sparse <- extract_sparse_parts(mod$mat)
			mod$stan <- list(N = N, S = S, P = 4, QS = 4,
				nz = length(sparse$w), x_w = sparse$w, x_v = sparse$v,
				x_u = sparse$u)
			
			# Generate model
			mod$b <- b <- c(runif(1,-2,2),runif(3,-1,1))
			mod$ss <- ss <- c(runif(1,0,1),runif(3,0,0.5))
			mod$os <- os <- rcorrmatrix(4)
			Ls <- t(chol(os))
			gs <- Ls %*% matrix(rnorm(S*4),4,S)
			for(i in 1:4) gs[i,] <- gs[i,] * ss[i]

			y <- as.vector(mod$mat %*% c(b,as.vector(gs)))

			if(linear){
				mod$se <- se <- runif(1,0,1)
				mod$pext <- NA
				y <- y + rnorm(N,0,se)
				mod$frame$y <- scale(y)[,1]
				mod$stan$y <- y
				mod$mu <- mod$stan$mu <- mean(y)
				mod$sdev <- mod$stan$sdev <- sd(y)
			} else {
				mod$se <- NA
				mod$mu <- 0
				mod$sdev <- 1
				mod$pext <- mean(abs(y)>5)
				mod$frame$y <- mod$stan$y <- rbinom(N,1,plogis(y))
			}

			return(mod)
		}

		### Fit single model with both lme4 and stan
		single_model <- function(linear){
			possfail_genmod <- tryCatch(mod <- generate_model(linear), error = function(e) e)
			while(inherits(possfail_genmod, "error")){
				possfail_genmod <- tryCatch(mod <- generate_model(linear), error = function(e) e)
			}

			if(linear){
				possfail_lme4 <- tryCatch(m_lme4 <- summary(lmer(formula, mod$frame)), error = function(e) e)
				possfail_stan <- tryCatch(m_stan <- sampling(object = gaustan, data = mod$stan, chains = 3,
					pars = c(keep,"res"), control = list(adapt_delta = 0.99)), error = function(e) e)
			} else {
				possfail_lme4 <- tryCatch(m_lme4 <- summary(glmer(formula, mod$frame, family = binomial)), error = function(e) e)
				possfail_stan <- tryCatch(m_stan <- sampling(object = binstan, data = mod$stan, chains = 3,
					pars = keep, control = list(adapt_delta = 0.99)), error = function(e) e)
			}
			
			if(inherits(possfail_lme4, "error") | inherits(possfail_stan, "error")){
				save(mod, possfail_lme4, possfail_stan, file = paste(getwd(),"/unconverged_thread_",
					thread,"_iteration_",b,"_linear_",linear,".rda",sep=""))
				return(NULL)
			} else {
				mess_lme4 <- m_lme4$optinfo$conv$lme4
				if(!is.null(names(mess_lme4))){
					mess_lme4 <- paste(mess_lme4$messages,collapse="; ")
					mess_lme4 <- str_replace_all(mess_lme4,",","")
					mess_lme4 <- str_replace_all(mess_lme4,"\n","")
				} else {
					mess_lme4 <- ""
				}
				vc <- data.frame(m_lme4$varcor)
				theta <- m_lme4$optinfo$val

				samp <- data.frame(do.call(rbind,
					args=get_sampler_params(m_stan,inc_warmup=FALSE)))
				nmtd <- sum(samp$treedepth__ > 10)
				ndiv <- sum(samp$divergent__)
				m_stan <- summary(m_stan,probs=c(.025,.975))$summary
				rhat <- quantile(m_stan[,"Rhat"],c(0,.25,.5,.75,1))
				neff <- quantile(m_stan[,"n_eff"],c(0,.25,.5,.75,1))
				m_stan <- m_stan[,1]
				
				if(mess_lme4 != "" | ndiv > 0 | as.numeric(rhat[5]) >= 1.1){
					save(mod, file = paste(getwd(),"/unconverged_thread_",
						thread,"_iteration_",b,"_linear_",linear,".rda",sep=""))
				}

				resb <- data.frame(
					thread = thread, iteration = b, linear = linear,
					pext = mod$pext, N = mod$N, S = mod$S, L = mod$L,
					balance = mod$balance, balance_ratio = mod$balance_ratio, mu = mod$mu, sdev = mod$sdev,
					b0 = mod$b[1], b0_lme4 = m_lme4$coef[1,1] * mod$sdev + mod$mu, b0_stan = m_stan["coef[1]"],
					b1 = mod$b[2], b1_lme4 = m_lme4$coef[2,1] * mod$sdev, b1_stan = m_stan["coef[2]"],
					b2 = mod$b[3], b2_lme4 = m_lme4$coef[3,1] * mod$sdev, b2_stan = m_stan["coef[3]"],
					b3 = mod$b[4], b3_lme4 = m_lme4$coef[4,1] * mod$sdev, b3_stan = m_stan["coef[4]"],
					s0 = mod$ss[1], s0_lme4 = vc[1,5] * mod$sdev, s0_stan = m_stan["s0"],
					s1 = mod$ss[2], s1_lme4 = vc[2,5] * mod$sdev, s1_stan = m_stan["s1"],
					s2 = mod$ss[3], s2_lme4 = vc[3,5] * mod$sdev, s2_stan = m_stan["s2"],
					s3 = mod$ss[4], s3_lme4 = vc[4,5] * mod$sdev, s3_stan = m_stan["s3"],
					r01 = mod$os[2,1], r01_lme4 = vc[5,5], r01_stan = m_stan["r01"],
					r02 = mod$os[3,1], r02_lme4 = vc[6,5], r02_stan = m_stan["r02"],
					r03 = mod$os[4,1], r03_lme4 = vc[7,5], r03_stan = m_stan["r03"],
					r12 = mod$os[3,2], r12_lme4 = vc[8,5], r12_stan = m_stan["r12"],
					r13 = mod$os[4,2], r13_lme4 = vc[9,5], r13_stan = m_stan["r13"],
					r23 = mod$os[4,3], r23_lme4 = vc[10,5], r23_stan = m_stan["r23"],
					se = mod$se, se_lme4 = ifelse(linear,vc[11,5]*mod$sdev,NA), se_stan = ifelse(linear,m_stan["res"],NA),
					mess_lme4, nmtd, ndiv,
					rhat_min = rhat[1], rhat_q25 = rhat[2], rhat_med = rhat[3], rhat_q75 = rhat[4], rhat_max = rhat[5],
					neff_min = neff[1], neff_q25 = neff[2], neff_med = neff[3], neff_q75 = neff[4], neff_max = neff[5],
					stringsAsFactors = FALSE)
				for(i in 1:10) resb[,paste("theta",i,sep="_")] <- theta[i]
				return(resb)
			}
		}

		## unpack arguments
		B <- pardat[[1]]
		thread <- pardat[[2]]
		gaucode <- pardat[[3]]
		bincode <- pardat[[4]]

		gaustan <- stan_model(model_name = "linear_sim", model_code = gaucode, save_dso = TRUE)
		binstan <- stan_model(model_name = "logistic_sim", model_code = bincode, save_dso = TRUE)

		### run B iterations of each of the 2 models based on H0 and output progress to thread-specific file
		results <- NULL
		progfile <- paste(getwd(),"/","imbalanced_sim_thread_",thread,".progress",sep="")
		csvfile <- paste(getwd(),"/","imbalanced_sim_thread_",thread,".csv",sep="")
		cat(paste(Sys.time(),"Starting Thread",thread),file=progfile,sep="\n",append=FALSE)
		
		for(b in 1:B){  # Each iteration has 2 calls to single_model, one for each H0 condition
			cat("\n",paste(Sys.time(),"Iteration",b,"/",B),file=progfile,sep="",append=TRUE)
			
			temp <- single_model(TRUE)
			if(is.null(temp)){
				results[nrow(results)+1,1:2] <- c(thread,b)
				results[nrow(results),3] <- TRUE
				cat(" ERROR ",file=progfile,sep="",append=TRUE)
			} else {
				results <- rbind(results, temp)
				cat(" . ",file=progfile,sep="",append=TRUE)
			}
			
			temp <- single_model(FALSE)
			if(is.null(temp)){
				cat(" ERROR ",file=progfile,sep="",append=TRUE)
				results[nrow(results)+1,1:2] <- c(thread,b)
				results[nrow(results),3] <- FALSE
			} else {
				results <- rbind(results, temp)
				cat(" . ",file=progfile,sep="",append=TRUE)
			}
			write.csv(results,row.names=F,quote=F,file=csvfile)
		}
		cat("\n",paste(Sys.time(),"Finished Thread",thread),file=progfile,sep="",append=TRUE)
		
		### return results
		return(results)
	}

	{### Stan code

{	bincode <- "// Stan code for logistic regression simulation
// adapted from Kimball, Shantz, Eager, and Roy (2016)

data {
  int<lower=2> N;  // number of observations
  int<lower=2> S;  // number of subjects

  int<lower=1> P;  // number of fixed effects
  int<lower=1,upper=P> QS;  // number of subject effects

  // sparse model matrix (CSR)
  int<lower=1> nz;  // number of non-zero elements in x
  vector[nz] x_w;  // non-zero elements in x
  int x_v[nz];  // column indices for x_w
  int x_u[N+1];  // row-start indices for x

  int<lower=0,upper=1> y[N];  // binary response
}

transformed data {
  int K;  // number of columns in x
  int SF;  // first subject effect column in x
  int SL;  // last subject effect column in x

  K = P + S * QS;
  SF = P + 1;
  SL = P + S * QS;
}

parameters {
  vector[P] beta;

  matrix[QS,S] gamma_subj_raw;
  vector<lower=0>[QS] sigma_subj;  // subject effect SDs
  cholesky_factor_corr[QS] omega_subj_raw;
}

transformed parameters {
  vector[K] coef;  // all coefficients
  vector[N] y_hat;  // predicted log-odds

  // transform fixed effects
  coef[1:P] = 2 * beta;

  // transform subject effects
  coef[SF:SL]
    = to_vector(rep_matrix(sigma_subj,S)
      .* (omega_subj_raw * gamma_subj_raw));

  // y_hat = x * coef
  y_hat = csr_matrix_times_vector(N,K,x_w,x_v,x_u,coef);
}

model {
  beta ~ normal(0,1);

  to_vector(gamma_subj_raw) ~ normal(0,1);
  sigma_subj ~ normal(0,1);
  omega_subj_raw ~ lkj_corr_cholesky(2);

  y ~ bernoulli_logit(y_hat);  // logistic model defined
}

generated quantities {
  real<lower=0> s0;
  real<lower=0> s1;
  real<lower=0> s2;
  real<lower=0> s3;
  real<lower=-1,upper=1> r01;
  real<lower=-1,upper=1> r02;
  real<lower=-1,upper=1> r03;
  real<lower=-1,upper=1> r12;
  real<lower=-1,upper=1> r13;
  real<lower=-1,upper=1> r23;

  s0 = sigma_subj[1];
  s1 = sigma_subj[2];
  s2 = sigma_subj[3];
  s3 = sigma_subj[4];
  {
    matrix[QS,QS] omega_subj;  // correlation in subject effects
    omega_subj = tcrossprod(omega_subj_raw);
    r01 = omega_subj[2,1];
    r02 = omega_subj[3,1];
    r03 = omega_subj[4,1];
    r12 = omega_subj[3,2];
    r13 = omega_subj[4,2];
    r23 = omega_subj[4,3];
  }
}
"
}

{	gaucode <- "// Stan code for linear regression simulation
// adapted from Kimball, Shantz, Eager, and Roy (2016)

data {
  int<lower=2> N;  // number of observations
  int<lower=2> S;  // number of subjects

  int<lower=1> P;  // number of fixed effects
  int<lower=1,upper=P> QS;  // number of subject effects

  // sparse model matrix (CSR)
  int<lower=1> nz;  // number of non-zero elements in x
  vector[nz] x_w;  // non-zero elements in x
  int x_v[nz];  // column indices for x_w
  int x_u[N+1];  // row-start indices for x

  vector[N] y;  // continuous response
  real mu;  // mean(y)
  real<lower=0> sdev;  // sd(y)
}

transformed data {
  int K;  // number of columns in x
  int SF;  // first subject effect column in x
  int SL;  // last subject effect column in x

  K = P + S * QS;
  SF = P + 1;
  SL = P + S * QS;
}

parameters {
  vector[P] beta_raw;
  real<lower=0> res_raw;

  matrix[QS,S] gamma_subj_raw;
  vector<lower=0>[QS] sigma_subj_raw;
  cholesky_factor_corr[QS] omega_subj_raw;
}

transformed parameters {
  vector[K] coef;  // all coefficients
  real<lower=0> res;  // residual standard error
  vector[N] y_hat;  // predicted log-odds

  // transform fixed effects
  coef[1:P] = beta_raw * 2 * sdev;
  coef[1] = coef[1] + mu;
  res = 0.5 * sdev * res_raw;

  // transform subject effects
  coef[SF:SL]
    = to_vector((rep_matrix(sigma_subj_raw,S) * sdev)
      .* (omega_subj_raw * gamma_subj_raw));

  // y_hat = x * coef
  y_hat = csr_matrix_times_vector(N,K,x_w,x_v,x_u,coef);
}

model {
  beta_raw ~ normal(0,1);
  res_raw ~ normal(0,1);

  to_vector(gamma_subj_raw) ~ normal(0,1);
  sigma_subj_raw ~ normal(0,1);
  omega_subj_raw ~ lkj_corr_cholesky(2);

  y ~ normal(y_hat,res);  // linear model defined
}

generated quantities {
  real<lower=0> s0;
  real<lower=0> s1;
  real<lower=0> s2;
  real<lower=0> s3;
  real<lower=-1,upper=1> r01;
  real<lower=-1,upper=1> r02;
  real<lower=-1,upper=1> r03;
  real<lower=-1,upper=1> r12;
  real<lower=-1,upper=1> r13;
  real<lower=-1,upper=1> r23;

  s0 = sigma_subj_raw[1] * sdev;
  s1 = sigma_subj_raw[2] * sdev;
  s2 = sigma_subj_raw[3] * sdev;
  s3 = sigma_subj_raw[4] * sdev;
  {
    matrix[QS,QS] omega_subj;  // correlation in subject effects
    omega_subj = tcrossprod(omega_subj_raw);
    r01 = omega_subj[2,1];
    r02 = omega_subj[3,1];
    r03 = omega_subj[4,1];
    r12 = omega_subj[3,2];
    r13 = omega_subj[4,2];
    r23 = omega_subj[4,3];
  }
}
"
}
	}

	if(threads==1){
		options(mc.cores = parallel::detectCores())
		results <- single_thread(list(B,threads,gaucode,bincode))
	} else {
		require(parallel)
		cl = makeCluster(rep("localhost", threads))
		
		# simulations per thread and thread number
		left <- B %% threads
		B <- rep(floor(B/threads),threads)
		if(left>0) B[1:left] <- B[1:left] + 1
		pardat <- list()
		for(i in 1:threads) pardat[[i]] <- list(B[i],i,gaucode,bincode)
		
		results <- parLapply(cl,pardat,single_thread)
		
		stopCluster(cl)
		
		results <- do.call(rbind,args=results)
	}
	
	return(results)
}
