#### compare_lme4_stan
#### simple logistic model

compare_lme4_stan <- function(B=10000,threads=16){
	# B is the number of models to run for each of the 8 model types
	# threads is the number of threads to split this task by; set to 1 to not parallelize
	# in which case the individual stan models will be parallelized

	single_thread <- function(pardat){
		require(lme4)
		require(rstan)
		
		### create data frame without response
		S <- 50
		I <- 30
		N <- S * I
		slist <- rep(c(1,2),length.out=S)
		ilist <- list(
			rep(c(1,-1),length.out=I),
			rep(c(-1,1),length.out=I))
		data <- data.frame(
			y = rep(0,N),
			subject = sort(rep(1:S,I)),
			item = rep(1:I,S),
			condition = unlist(ilist[slist]))
		formula <- y ~ condition + (1 + condition | subject) + (1 + condition | item)
		x <- model.matrix(nobars(formula),data)
		z <- mkReTrms(findbars(formula),data)
		xz <- cbind(x,t(z$Ztlist[[1]]),t(z$Ztlist[[2]]))
		xzsparse <- extract_sparse_parts(xz)
		standata <- list(N = N, S = S, I = I, P = 2, QS = 2, QI = 2, y = 0,
			nz = length(xzsparse$w), x_w = xzsparse$w, x_v = xzsparse$v,
			x_u = xzsparse$u)
		keep <- c("s0","s1","rs","i0","i1","ri","coef")

		### Generate model
		generate_model <- function(H0){
			b0 <- runif(1,-1,1)
			b1 <- ifelse(H0,0,0.8)

			ss <- runif(1,0,1.5)
			ss[2] <- runif(1,0,0.75)

			si <- runif(1,0,1)
			si[2] <- runif(1,0,0.5)

			rs <- runif(1,-0.9,0.9)
			Ls <- t(chol(matrix(c(1,rs,rs,1),2,2)))

			ri <- runif(1,-0.9,0.9)
			Li <- t(chol(matrix(c(1,ri,ri,1),2,2)))

			gs <- Ls %*% matrix(rnorm(S*2),2,S)
			gs[1,] <- gs[1,] * ss[1]
			gs[2,] <- gs[2,] * ss[2]

			gi <- Li %*% matrix(rnorm(I*2),2,I)
			gi[1,] <- gi[1,] * si[1]
			gi[2,] <- gi[2,] * si[2]

			bg <- c(b0,b1,as.vector(gs),as.vector(gi))
			lo <- as.vector(xz %*% bg)
			pext <- mean(abs(lo)>5)
			y <- rbinom(N,1,plogis(lo))

			return(list(y = y, s0 = ss[1], s1 = ss[2],
				rs = rs, i0 = si[1], i1 = si[2], ri = ri,
				pext = pext, b0 = b0, b1 = b1))
		}
		
		### Fit single model with both lme4 and stan
		single_model <- function(H0){
			mod <- generate_model(H0)
			y <- data$y <- standata$y <- mod$y
			
			m_lme4 <- summary(glmer(formula, data, family = binomial))
			p_H0_lme4 <- m_lme4$coef[2,4]
			mess_lme4 <- m_lme4$optinfo$conv$lme4
			if(!is.null(names(mess_lme4))){
				mess_lme4 <- paste(mess_lme4$messages,collapse="; ")
				mess_lme4 <- str_replace_all(mess_lme4,",","")
				mess_lme4 <- str_replace_all(mess_lme4,"\n","")
			} else {
				mess_lme4 <- ""
			}
			vc <- data.frame(m_lme4$varcor)
			
			m_stan <- sampling(object = stanmod, data = standata, pars = keep, chains = 3)
			samp <- data.frame(do.call(rbind,args=get_sampler_params(m_stan,inc_warmup=FALSE)))
			nmtd <- sum(samp$treedepth__ > 10)
			ndiv <- sum(samp$divergent__)
			p_HA_stan <- mean(extract(m_stan,pars="coef")$coef[,2]>0)
			m_stan <- summary(m_stan,probs=c(.025,.975))$summary
			ci_0_stan <- m_stan["coef[2]",4]<0 & m_stan["coef[2]",5]>0
			rhat <- quantile(m_stan[,"Rhat"],c(0,.25,.5,.75,1))
			neff <- quantile(m_stan[,"n_eff"],c(0,.25,.5,.75,1))
			m_stan <- m_stan[,1]
			
			if(mess_lme4 != "" | ndiv > 0 | as.numeric(rhat[5]) >= 1.1){
				save(y, file = paste(getwd(),"/unconverged_thread_",thread,"_iteration_",b,"_H0_",H0,".rda",sep=""))
			}

			return(data.frame(
				H0, pext = mod$pext,
				b0 = mod$b0, b0_lme4 = m_lme4$coef[1,1], b0_stan = m_stan["coef[1]"],
				b1 = mod$b1, b1_lme4 = m_lme4$coef[2,1], b1_stan = m_stan["coef[2]"],
				s0 = mod$s0, s0_lme4 = vc[1,5], s0_stan = m_stan["s0"],
				s1 = mod$s1, s1_lme4 = vc[2,5], s1_stan = m_stan["s1"],
				rs = mod$rs, rs_lme4 = vc[3,5], rs_stan = m_stan["rs"],
				i0 = mod$i0, i0_lme4 = vc[4,5], i0_stan = m_stan["i0"],
				i1 = mod$i1, i1_lme4 = vc[5,5], i1_stan = m_stan["i1"],
				ri = mod$ri, ri_lme4 = vc[6,5], ri_stan = m_stan["ri"],
				p_H0_lme4, p_HA_stan, ci_0_stan,
				mess_lme4, nmtd, ndiv,
				rhat_min = rhat[1], rhat_q25 = rhat[2], rhat_med = rhat[3], rhat_q75 = rhat[4], rhat_max = rhat[5],
				neff_min = neff[1], neff_q25 = neff[2], neff_med = neff[3], neff_q75 = neff[4], neff_max = neff[5]))
		}

		## unpack arguments
		B <- pardat[[1]]
		thread <- pardat[[2]]
		stanmod <- pardat[[3]]

		### run B iterations of each of the 2 models based on H0 and output progress to thread-specific file
		results <- NULL
		progfile <- paste(getwd(),"/","simple_logistic_sim_thread_",thread,".progress",sep="")
		csvfile <- paste(getwd(),"/","simple_logistic_sim_thread_",thread,".csv",sep="")
		cat(paste(Sys.time(),"Starting Thread",thread),file=progfile,sep="\n",append=FALSE)
		
		for(b in 1:B){  # Each iteration has 2 calls to single_model, one for each H0 condition
			cat("\n",paste(Sys.time(),"Iteration",b,"/",B),file=progfile,sep="",append=TRUE)
			results <- rbind(results,single_model(TRUE))
			cat(" . ",file=progfile,sep="",append=TRUE)
			results <- rbind(results,single_model(FALSE))
			cat(" . ",file=progfile,sep="",append=TRUE)
			write.csv(results,row.names=F,quote=F,file=csvfile)
		}
		cat("\n",paste(Sys.time(),"Finished Thread",thread),file=progfile,sep="",append=TRUE)
		
		### return results
		return(results)
	}

	{### Stan code
	stancode <- "// Stan code for simple logistic regression simulation
// adapted from Kimball, Shantz, Eager, and Roy (2016)

data {
  int<lower=2> N;  // number of observations
  int<lower=2> S;  // number of subjects
  int<lower=2> I;  // number of items

  int<lower=1> P;  // number of fixed effects
  int<lower=1,upper=P> QS;  // number of subject effects
  int<lower=1,upper=P> QI;  // number of item effects

  int<lower=0,upper=1> y[N];  // binary response

  // sparse model matrix (CSR)
  int<lower=1> nz;  // number of non-zero elements in x
  vector[nz] x_w;  // non-zero elements in x
  int x_v[nz];  // column indices for x_w
  int x_u[N+1];  // row-start indices for x
}

transformed data {
  int K;  // number of columns in x
  int SF;  // first subject effect column in x
  int SL;  // last subject effect column in x
  int IF;  // first item effect column in x
  int IL;  // last item effect column in x

  K = P + S * QS + I * QI;
  SF = P + 1;
  SL = P + S * QS;
  IF = SL + 1;
  IL = SL + I * QI;
}

parameters {
  vector[P] beta;

  matrix[QS,S] gamma_subj_raw;
  vector<lower=0>[QS] sigma_subj;  // subject effect SDs
  cholesky_factor_corr[QS] omega_subj_raw;

  matrix[QI,I] gamma_item_raw;
  vector<lower=0>[QI] sigma_item;  // item effect SDs
  cholesky_factor_corr[QI] omega_item_raw;
}

transformed parameters {
  vector[K] coef;  // all coefficients
  vector[N] y_hat;  // predicted log-odds

  // transform fixed effects
  coef[1:P] = beta;

  // transform subject effects
  coef[SF:SL]
    = to_vector(rep_matrix(sigma_subj,S)
      .* (omega_subj_raw * gamma_subj_raw));

  // transform item effects
  coef[IF:IL]
    = to_vector(rep_matrix(sigma_item,I)
      .* (omega_item_raw * gamma_item_raw));

  // y_hat = x * coef
  y_hat = csr_matrix_times_vector(N,K,x_w,x_v,x_u,coef);
}

model {
  beta ~ normal(0,1);

  to_vector(gamma_subj_raw) ~ normal(0,1);
  sigma_subj ~ normal(0,1);
  omega_subj_raw ~ lkj_corr_cholesky(2);

  to_vector(gamma_item_raw) ~ normal(0,1);
  sigma_item ~ normal(0,1);
  omega_item_raw ~ lkj_corr_cholesky(2);

  y ~ bernoulli_logit(y_hat);  // logistic model defined
}

generated quantities {
  real<lower=0> s0;
  real<lower=0> s1;
  real<lower=0> i0;
  real<lower=0> i1;
  real<lower=-1,upper=1> rs;
  real<lower=-1,upper=1> ri;

  s0 = sigma_subj[1];
  s1 = sigma_subj[2];
  i0 = sigma_item[1];
  i1 = sigma_item[2];
  {
    matrix[QS,QS] omega_subj;  // correlation in subject effects
    matrix[QI,QI] omega_item;  // correlation in item effects
    omega_subj = tcrossprod(omega_subj_raw);
    omega_item = tcrossprod(omega_item_raw);
    rs = omega_subj[1,2];
    ri = omega_item[1,2];
  }
}
"
	}

	require(rstan)
	stanmod <- stan_model(model_name = "simple_logistic_sim", model_code = stancode, save_dso = TRUE)

	if(threads==1){
		options(mc.cores = parallel::detectCores())
		results <- single_thread(list(B,threads,stanmod))
	} else {
		require(parallel)
		cl = makeCluster(rep("localhost", threads))
		
		# simulations per thread and thread number
		left <- B %% threads
		B <- rep(floor(B/threads),threads)
		if(left>0) B[1:left] <- B[1:left] + 1
		pardat <- list()
		for(i in 1:threads) pardat[[i]] <- list(B[i],i,stanmod)
		
		results <- parLapply(cl,pardat,single_thread)
		
		stopCluster(cl)
		
		results <- do.call(rbind,args=results)
	}
	
	return(results)
}

results <- compare_lme4_stan()
