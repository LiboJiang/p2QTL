###--- written on 2014.12.12 by yemeixia ---###
library("mvtnorm");
get_Parent_RIL <- function(RILs_file){
		parent <- c("B73",
								"B97","CML103","CML228","CML247","CML277","CML322","CML333","CML52","CML69","HP301",
								"Il14H","Ki11","Ki3","Ky21","M162W","M37W","Mo18W","MS71","NC350","NC358",
								"Oh43","OH7B","P39","Tx303","Tzi8"
								);
		
		RILs <- as.matrix( read.table(RILs_file, head=F) );
		n <- dim(RILs)[1];
		m <- dim(RILs)[2];
		RILs <- as.vector( RILs[n,n:m] );
		
		all <- c(parent, RILs);
		return(all);
}

get_marker_name <- function(map_file){
		map <- as.matrix( read.table(map_file, head=T) );
		marker <- map[,1];
		
		return(marker);
}

process_pheno <- function(pheno_file, ind_name, trait_name, env_name){
		pheno <- as.matrix( read.table(pheno_file, head=T) );
		
		m_id <- which( pheno[1,] == trait_name );
		pheno <- pheno[ ,c(1,m_id)];
		final_pheno <- pheno[1:2,];
		
		for( i in 1:length(ind_name) ){
				cur_id <- ind_name[i];
				n_id <- which(pheno[,1] == cur_id);
				final_pheno <- rbind( final_pheno, pheno[n_id, ] );
		}
		
		final_pheno1 <- NULL;
		for( i in 1:length(env_name) ){
				cur_id <- env_name[i];
				n_id <- which(final_pheno[2,] == cur_id);
				final_pheno1 <- cbind(final_pheno1, final_pheno[,n_id]);
		}
		
		rownames(final_pheno1) <- final_pheno[,1];
		final_pheno1 <- final_pheno1[-c(1,2), ];
		rname <- rownames(final_pheno1);
	
		final_pheno1 <- matrix(as.numeric(final_pheno1), nrow=nrow(final_pheno1));
		final_pheno1 <- abs(final_pheno1);
		rownames(final_pheno1) <- rname;
		
		
		return(final_pheno1);
}

process_geno <- function(ind, mk, raw_snp_file, imputed_snp_file){
		########------- extract specific snp bases for parents at all markers ---###
		snp_base <- as.matrix( read.table(raw_snp_file, head=T) );
		rownames(snp_base) <- snp_base[,1];
		snp_base <- snp_base[,-1];
		
		parent_snp <- NULL;
		for( i in 1:length(mk) ){
				mk_name <- mk[i];
				mk_name <- strsplit(mk_name, "/")[[1]][1];
				
				mk_id <- which( rownames(snp_base)==mk_name );
				parent_snp <- rbind(parent_snp, snp_base[mk_id,]);
				rownames(parent_snp)[i] <- mk_name;
		}
		
		pr_snp <- NULL;
		parent_ind <- ind[1:26];
		for( i in 1:length(parent_ind) ){
				pr_name <- parent_ind[i];
				pr_id <- which( colnames(parent_snp) == pr_name );
				pr_snp <- cbind(pr_snp, parent_snp[ ,pr_id]);
				colnames(pr_snp)[i] <- pr_name;
		}
		pr_snp <- t(pr_snp);
		
		########------- transform bases of parents to 0,1,2 ---###
		pr_snp_num <- matrix(0, dim(pr_snp)[1], dim(pr_snp)[2]);
		rownames(pr_snp_num) <- rownames(pr_snp);
		colnames(pr_snp_num) <- colnames(pr_snp);
		
		pr_snp_num[1,] <- rep(0, times=dim(pr_snp)[2]);
		for(j in 1:dim(pr_snp)[2]){
			for(i in 2:dim(pr_snp)[1]){
					if(pr_snp[i,j] == pr_snp[1,j]){ pr_snp_num[i,j] <- 0; }
					else{
							bases <- strsplit(pr_snp[i,j],"/")[[1]];
							if(bases[1] != bases[2]){ pr_snp_num[i,j] <- 1; }
							if(bases[1] == "-"){ pr_snp_num[i,j] <- 0.3; }
							else{ pr_snp_num[i,j] <- 2; }
					}
			}
		}
		
		########------- transform bases of parents to 0,1,2 ---###
		impute_RIL0 <- as.matrix( read.table(imputed_snp_file, head=F) );
		rownames(impute_RIL0) <- impute_RIL0[,1];
		impute_RIL0 <- impute_RIL0[,-1];
		
		impute_RIL2 <- impute_RIL0[4:dim(impute_RIL0)[1],];
		impute_RIL21 <- matrix(as.numeric(impute_RIL2), nrow=nrow(impute_RIL2));
		rownames(impute_RIL21) <- rownames(impute_RIL2);
	
		impute_RILs <- rbind(pr_snp_num, impute_RIL21);
		return(impute_RILs);
}

process_map <- function(imputed_snp_file){
		impute_RIL0 <- as.matrix( read.table(imputed_snp_file, head=F) );
		rownames(impute_RIL0) <- impute_RIL0[,1];
		impute_RIL0 <- impute_RIL0[,-1];
		
		final <- impute_RIL0[1:3,];
		return(final);
}

############--- get data --------#####################
trait_name <- c("PlantHeight");
#env_name <- c("06A","07A","06CL1","07CL1","26M3","27M3","06FL1","07FL1","06PR","065","07U");
env_name <- c("06A","07A","06CL1","07CL1","06FL1","06PR","065","07U");

all_ind_name <- get_Parent_RIL("RILs_for_NAM_map.txt");
all_mk_name  <- get_marker_name("NAM_map.txt");
all_pheno <- process_pheno("Trait_maize282NAM.txt", all_ind_name, trait_name, env_name);
all_geno <- process_geno(all_ind_name, all_mk_name, "NAM_genos_raw.txt", "NAM_genos_imputed.txt");
maps <- process_map("NAM_genos_imputed.txt");

############------ calculate EI -------################
sums <- apply(all_pheno, 1, sum);
valid_id <- which(sums != "NaN");
new_all_pheno <- all_pheno[valid_id, ];

mean_each <- apply(new_all_pheno, 2, mean);
mean_all  <- sum(new_all_pheno)/(dim(new_all_pheno)[1] * dim(new_all_pheno)[2]);
dfrt <- mean_each - mean_all;

xuhao <- order(dfrt);
env_name <- env_name[xuhao];
env_indice <- sort(dfrt);

renew_all_pheno <- new_all_pheno[ ,xuhao];
renew_all_geno  <-  all_geno[valid_id, ];
for(i in 1:dim(renew_all_geno)[1]){
		for(j in 1:dim(renew_all_geno)[2]){
				if(renew_all_geno[i,j] != 0 && renew_all_geno[i,j] != 2){ renew_all_geno[i,j] <- -1; }
		}
}

#####---- test how many individuals had missing genotypes at flanking markers --#############
#impute <- rep(0, times=1105);
#for(j in 1:1105){
#	for(i in 1:3502){
#		if(renew_all_geno[i,j] == 0 && renew_all_geno[i,j+1] == 0){ impute[j] <- impute[j] + 0; }
#		else if(renew_all_geno[i,j] == 2 && renew_all_geno[i,j+1] == 2){ impute[j] <- impute[j] + 0; }
#		else if(renew_all_geno[i,j] == 0 && renew_all_geno[i,j+1] == 2){ impute[j] <- impute[j] + 0; }
#		else if(renew_all_geno[i,j] == 2 && renew_all_geno[i,j+1] == 0){ impute[j] <- impute[j] + 0; }
#		else{ impute[j] <- impute[j] + 1; }
#	}
#}
#fivenum(impute); 94, 252, 295, 363, 887

#########---- define functions ----------#################
get_mu <- function(a, b, xj){
		mu <- a * (b^xj);
		return(mu);
}

get_var <- function(rho, sig, xj){
		n <- length(xj);
		sigma <- matrix(0, n, n);
		for(i in 1:n){
			for(j in 1:n){
					if(i==j){ sigma[i,j] <- 1; }
					else{ sigma[i,j] <- rho^( abs(xj[i]-xj[j]) ); }
			}
		}
		
		sigma <- sig * sigma;
		return(sigma);
}

get_loglike1 <- function(parin, y, markers, xj){
		a0 <- parin[1];	b0 <- parin[2];			mu0 <- get_mu(a0, b0, xj);
		a2 <- parin[3];	b2 <- parin[4];			mu2 <- get_mu(a2, b2, xj);
		rho <- parin[5];sig <- parin[6];	sigma <- get_var(rho, sig, xj);
		
		if(rho < 0.1 || rho > 0.98){ return(NaN); }
		if(sig < 0){ return(NaN); }
		
		id2 <- which(markers == 2);
		id0 <- which(markers == 0);
		y2  <- y[id2,];
		y0  <- y[id0,];
		
		fy2 <- dmvnorm(y2, mu2, sigma, log=TRUE);
		fy0 <- dmvnorm(y0, mu0, sigma, log=TRUE);
		
		loglike <- 0 - (sum(fy2) + sum(fy0));
		return( loglike );
}

get_loglike1_a <- function(pars, parin, y, markers, xj){
		a0 <- pars[1];			a2 <- pars[2];	
		
		b0 <- parin[1];			b2 <- parin[2];	
		rho <- parin[3];		sig <- parin[4];	
		
		if(rho < 0.1 || rho > 0.98){ return(NaN); }
		if(sig < 0){ return(NaN); }
		
		mu0 <- get_mu(a0, b0, xj);
		mu2 <- get_mu(a2, b2, xj);
		sigma <- get_var(rho, sig, xj);
		
		id2 <- which(markers == 2);
		id0 <- which(markers == 0);
		y2  <- y[id2,];
		y0  <- y[id0,];
		
		fy2 <- dmvnorm(y2, mu2, sigma, log=TRUE);
		fy0 <- dmvnorm(y0, mu0, sigma, log=TRUE);
		
		loglike <- 0 - (sum(fy2) + sum(fy0));
		return( loglike );
}

get_loglike1_var <- function(par_var, parin, y, markers, xj){
		rho <- par_var[1];	sig <- par_var[2];		sigma <- get_var(rho, sig, xj);
		a0 <- parin[1];	b0 <- parin[2];			mu0 <- get_mu(a0, b0, xj);
		a2 <- parin[3];	b2 <- parin[4];			mu2 <- get_mu(a2, b2, xj);
		
		if(rho < 0.1 || rho > 0.98){ return(NaN); }
		if(sig < 0){ return(NaN); }
		
		id2 <- which(markers == 2);
		id0 <- which(markers == 0);
		y2  <- y[id2,];
		y0  <- y[id0,];
		
		fy2 <- dmvnorm(y2, mu2, sigma, log=TRUE);
		fy0 <- dmvnorm(y0, mu0, sigma, log=TRUE);
		
		loglike <- 0 - (sum(fy2) + sum(fy0));
		return( loglike );
}

get_loglike0 <- function(parin, y, xj){
		a   <- parin[1];		b <- parin[2];		mu <- get_mu(a, b, xj);
		rho <- parin[3];		sig <- parin[4];	sigma <- get_var(rho, sig, xj);
		
		if(rho < 0.1 || rho > 0.98){ return(NaN); }
		if(sig < 0){ return(NaN); }
		#cat("parin",parin,"\n");
		
		fy <- dmvnorm(y, mu, sigma, log=TRUE);
		loglike <- 0 - sum(fy);
		
		#cat("loglike0",loglike,"\n");
		return( loglike );
}

get_loglike0_var <- function(par_var, parin, y, xj){
		rho <- par_var[1];	sig <- par_var[2];	sigma <- get_var(rho, sig, xj);
		a   <- parin[1];		b <- parin[2];		mu <- get_mu(a, b, xj);
		
		if(rho < 0.1 || rho > 0.98){ return(NaN); }
		if(sig < 0){ return(NaN); }
		
		fy <- dmvnorm(y, mu, sigma, log=TRUE);
		loglike <- 0 - sum(fy);
		
		return( loglike );
}

qtlmodel.get_est_LR2 <- function(dat, id_x, parin, xj){
	y <- as.matrix( dat$phenos );
	snp_x <- dat$genos[, id_x]; 
	
	##########-------- remove missing data --------###############
	missing_x <- which( snp_x == -1 );
	
	missing <- sort( unique(c(missing_x)) );
	if( length(missing) > 0){
		snp_x <- snp_x[ -(missing) ];
		y <- y[ -(missing), ];
	}
	
	cat("N of individuals:",dim(y)[1],"\n");
	#markers <- cbind( snp_x, snp_y );
	markers <- snp_x;
	##########-------- logL for H0 hypothesis --------#########
	min.par <- c();
	min.val <- Inf;
	loop <- 1;
	while(loop <= 5){
			h0 <- try( optim( parin, get_loglike0, y=y, xj=xj, method="BFGS" ), TRUE);
			if( class(h0)=="try-error" ){
					if(length(min.par)==0) { parin <- c( parin[1:2] * runif(2, 0.99, 1.01), 0.8, 500 ); }
          else { parin <- min.par * runif(4, 0.99, 1.01); }
          next;		
			}
			
			if( length(min.par) == 0 ){	min.par <- h0$par;	min.val <- h0$value;	}
			else{
					if( h0$value <= min.val ){	min.par <- h0$par;	min.val <- h0$value; }
			}
			
			parin <- min.par * runif(4, 0.98, 1.02);
			loop <- loop + 1;
			
			cat(min.val, " ", min.par, "\n");
	}
	par_h0 <- min.par;
	val_h0 <- min.val;
	
	par_var <- par_h0[3:4];
	parin   <- par_h0[1:2];
	h0 <- try( optim( par_var, get_loglike0_var, parin=parin, y=y, xj=xj, method="L-BFGS-B",lower=c(0.5,400), upper=c(0.98,700) ), TRUE);
	if( h0$value <= min.val ){	min.par <- c(par_h0[1:2], h0$par);	min.val <- h0$value; }
	par_h0 <- min.par;
	val_h0 <- min.val;
	cat(min.val, " ", min.par, "\n");
	
	##########-------- logL for H1 hypothesis --------#########
	parin1 <- c( par_h0[1:2]*runif(2, 1.00, 1.00), 
							 par_h0[1:2]*runif(2, 1.00, 1.00), 
							 par_h0[3:4] 
						 );

	min.par <- c();
	min.val <- Inf;
	loop <- 1;
	while(loop <= 5){
			if(loop %% 2 == 0){ h1 <- try( optim( parin1, get_loglike1, y=y, markers=markers, xj=xj, method="Nelder-Mead" ), TRUE); }
			if(loop %% 2 == 1){ h1 <- try( optim( parin1, get_loglike1, y=y, markers=markers, xj=xj, method="L-BFGS-B", lower=c(145,0.99,145,0.99,0.1,300), upper=c(160,1.1,160,1.1,0.98,700) ), TRUE); }
			if( class(h1)=="try-error" ){
					if(length(min.par)==0) { parin <- c( parin[1:4] * runif(4), 0.7, 500 ); }
          else { parin <- min.par * runif(6, 0.98, 1.02); }
          next;		
			}
			
			if( length(min.par) == 0 ){	min.par <- h1$par;	min.val <- h1$value;	}
			else{
					if( h1$value <= min.val ){	min.par <- h1$par;	min.val <- h1$value; }
			}
			
			parin <- min.par * runif(6, 0.95, 1.05);
			loop <- loop + 1;
			
			cat(min.val, " ", min.par, "\n");
	}
	par_h1 <- min.par;
	val_h1 <- min.val;
	
	pars <- c(par_h1[1], par_h1[3]);
	parin <- c(par_h1[2], par_h1[4], par_h1[5:6]);
	loop <- 1;
	min.par <- c();
	min.val <- Inf;
	while(loop <= 5){
			if(loop %% 2 == 0){ h1 <- try( optim( pars, get_loglike1_a, parin=parin, y=y, markers=markers, xj=xj, method="Nelder-Mead" ), TRUE); }
			if(loop %% 2 == 1){ h1 <- try( optim( pars, get_loglike1_a, parin=parin, y=y, markers=markers, xj=xj, method="BFGS" ), TRUE); }
			if( class(h1)=="try-error" ){
					if(length(min.par)==0) { pars <- c( pars * runif(2)); }
          else { pars <- min.par * runif(2, 0.98, 1.02); }
          next;		
			}
			
			if( length(min.par) == 0 ){	min.par <- h1$par;	min.val <- h1$value;	}
			else{
					if( h1$value <= min.val ){	min.par <- h1$par;	min.val <- h1$value; }
			}
			
			pars <- min.par * runif(2, 0.95, 1.05);
			loop <- loop + 1;
			
			cat("optim a", min.val, " ", min.par, "\n");
	}
	
	par_var <- parin[3:4]
	parin   <- c( min.par[1], parin[1], min.par[2], parin[2] );
	
	#par_var <- par_h1[5:6];
	#parin   <- par_h1[1:4];
	h1 <- try( optim(par_var, get_loglike1_var, parin=parin, y=y, markers=markers, xj=xj), TRUE );
	if( h1$value <= min.val ){	min.par <- c(par_h1[1:4], h1$par);	min.val <- h1$value; }
	par_h1 <- min.par;
	val_h1 <- min.val;
	cat(min.val, " ", min.par, "\n");
	
	##########-------- calculate LR2 ----------------#########
	final <- c( 2*( val_h0 - val_h1 ), val_h1, par_h1 );
	return ( final );
}

#########---- real analysis ----------#################

N_marker <- dim(maps)[2];
n <- dim(renew_all_pheno)[1];
#perm <- sample( seq(1:n), n, replace=FALSE );
#perm_pheno <- renew_all_pheno[perm,];
#dat <- list(phenos = perm_pheno, genos = renew_all_geno);
dat <- list(phenos = renew_all_pheno, genos = renew_all_geno);

id_xy <- seq(1:N_marker);
xj <- env_indice;
parin <- c( 153.5724, 1.007277, 0.8, 550 );
LR <- matrix( 0, length(id_xy), 8 );
for(i in 1:length(id_xy)){
		cat("marker:",i,"\t");
		LR[i, ] <- qtlmodel.get_est_LR2(dat, id_xy[i], parin, xj);
		if( i %% 10 == 0 ){ save.image("single_1marker.Rdata"); }
}

cat("max of LR2", max(LR[,1]), "\n");
save.image("single_1marker.Rdata");

