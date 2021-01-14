# generate table #
rm(list = ls())

if (!require("MASS")) install.packages("MASS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("akima")) install.packages("akima")
if (!require("plotly")) install.packages("plotly")
if (!require("reshape2")) install.packages("reshape2")
if (!require("plm")) install.packages("plm")

set.seed(100)
  
#######################
#####     DGP     #####
#######################

DGP1 <- function(T_, N, beta_true){
  # Set parameters
  K <- 3
  range_x <- c(0,1)
  alpha_i <- 1
  mu_u_i <- rep(0, T_)
  sigma_u_i <- diag(rep(1, T_))
  # Simulate data 
  X <- lapply(as.list(1:N), function(i) array(runif(K*T_,range_x[1],range_x[2]), dim = c(T_, K)))
  U <- lapply(as.list(1:N), function(i) mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i))
  Y <- mapply(function(X_i, u_i) alpha_i + X_i %*% beta_true + u_i, X, U, SIMPLIFY=F)
     # Alternatively, we can write in the following way:
     # X_i <- array(runif(K*T_,range_x[1],range_x[2]), dim = c(T_, K))
     # u_i <- mvrnorm(n = 1, mu = mu_u_i, Sigma = sigma_u_i)
     # Y_i <- alpha_i + X_i %*% beta_true + u_i
  # Save all results to data frame
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = NA) %>% cbind(matrix(NA,nrow=N*T_,ncol=K))
  colnames(df)[4:(3+K)] <- paste0("x_it_",c(0:(K-1)))
  # Loop over N & T_
  for(i in 1:N){
    X_i <- X[[i]]
    u_i <- U[[i]]
    Y_i <- Y[[i]]
    for(t in 1:T_){
      df[(i-1)*T_ + t, 4:(3+K)] <- X_i[t,]
      df$y_it[(i-1)*T_ + t] <- Y_i[t]
    }
  }
  return(list(X_list = X, Y_list = Y, df=df))
}

DGP2 <- function(T_, N, beta_true){
  if(T_ < 2 || N < 2)
    stop("T or N less than 2")
  alpha_i <- runif(N, min = 0, max = 1)
  eps_t <- rnorm(T_, mean = 0, sd = 1)
  z0 <- 0
  z_t <- c(0.5*z0 + eps_t[1])
  for(i in 2:T_)
    z_t[i] <- 0.5*z_t[i-1] + eps_t[i]
  v_i <- rnorm(N, mean = 0, sd = 2)
  e_t <- rnorm(T_ + 1, mean = 0, sd = 1)
  mu <- 5
  w_t <- mu + e_t[2:(T_+1)] + e_t[1:T_] / 3
  x0 <- 1
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = NA,
                   x_it_0 = x0,
                   x_it_1 = NA)
  for(i in 1:N){
    for(t in 1:T_){
      x_it <- sqrt(alpha_i[i]) * z_t[t]
      df$x_it_1[(i-1)*T_ + t] <- x_it
      u_it <- v_i[i] + w_t[t]
      df$y_it[(i-1)*T_ + t] <- alpha_i[i] + beta_true[1] * x0 + beta_true[2] * x_it + u_it
    }
  }
  # Save results to lists
  X_list <- list()
  Y_list <- list()
  K <- 2
  for(i in 1:N){
    X_list[[i]] <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+K)])
    Y_list[[i]] <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
  }
  return(list(df=df, X_list=X_list, Y_list=Y_list))
}

DGP3 <- function(T_, N, beta_true){
  # Set parameters
  r<-2
  K<-3
  mu <- beta_true[1]
  gamma <- 2
  delta <- 4
  iota <- rep(1, r)
  mu1 <- mu2 <- c1 <- c2 <- 1
  # Generate variables
  F_t <- matrix(rnorm(n = r*T_, mean = 0, sd = 1), nrow=r, ncol=T_)
  Lambda_i <- matrix(rnorm(n = r*N, mean = 0, sd = 1), nrow=r, ncol=N)
  Eta_it_1 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eta_it_2 <- matrix(rnorm(n = T_*N, mean = 0, sd = 1), nrow=N, ncol=T_)
  Eps_it <- matrix(rnorm(n = T_*N, mean = 0, sd = 4), nrow=N, ncol=T_)
  e_i <- rnorm(n = N, mean = 0, sd = 1)
  eta_t <- rnorm(n = T_, mean = 0, sd = 1)
  # Calculate intermediate variables
  iota_Lambda_i <- crossprod(iota, Lambda_i)
  iota_F_t <- crossprod(iota, F_t)
  dim(iota_Lambda_i) <- c(N, 1)
  dim(iota_F_t) <- c(1, T_)
  iota_Lambda_it <- iota_Lambda_i %*% rep(1, T_)
  iota_F_it <- rep(1, N) %*% iota_F_t
  x_i <- iota_Lambda_i + e_i
  w_t <- iota_F_t + eta_t
  gamma_x_it <- gamma * x_i %*% rep(1, T_)
  delta_w_it <- delta * rep(1, N) %*% w_t
  Lambda_F_it <- crossprod(Lambda_i, F_t)
  # Simulate data
  X_it_1 <- mu1 + c1 * Lambda_F_it + iota_Lambda_it + iota_F_it + Eta_it_1
  X_it_2 <- mu2 + c2 * Lambda_F_it + iota_Lambda_it + iota_F_it + Eta_it_2
  Y_it <- beta_true[2]*X_it_1 + beta_true[3]*X_it_2 + mu + gamma_x_it + delta_w_it + 
    Lambda_F_it + Eps_it
  # Save all results to data frame
  df <- data.frame(i = rep(c(1:N), each = T_),
                   t = rep(c(1:T_), times = N),
                   y_it = as.vector(t(Y_it)),
                   x_it_0 = 1,
                   x_it_1 = as.vector(t(X_it_1)),
                   x_it_2 = as.vector(t(X_it_2)))
  # Save results to lists
  X_list <- list()
  Y_list <- list()
  for(i in 1:N){
    X_list[[i]] <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+K)])
    Y_list[[i]] <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
  }
  return(list(df=df, X_list=X_list, Y_list=Y_list))
}

#####################
####   Method    ####
#####################

#####  Introduction   #####

### Fixed Effects Model ###
OLS_FE <- function(df){
  # Initialize
  N <- length(unique(df$i))
  T_ <- length(unique(df$t))
  K <- ncol(df) - 3
  A <- array(0, dim = c(K,K))
  B <- rep(0, K)
  # Loop over N & T_
  for(i in 1:N){
    X_i <- as.matrix(df[((i-1)*T_+1):(i*T_), 4:(3+K)])
    Y_i <- as.matrix(df$y_it[((i-1)*T_+1):(i*T_)])
    x_i_mean <- colMeans(X_i)
    y_i_mean <- mean(Y_i)
    for(t in 1:T_){
      A <- A + (X_i[t,] - x_i_mean) %*% t(X_i[t,] - x_i_mean)
      B <- B + (X_i[t,] - x_i_mean) * (Y_i[t] - y_i_mean)
    }
  }
  beta_hat_fe <- solve(A) %*% B
  return(list(beta_hat = beta_hat_fe))
}

# alternatively, we can use plm package for estimation
OLS_FE2 <- function(df){
  K <- ncol(df) - 3
  data <- pdata.frame(df,index=c("i","t"))
  formulate <- reformulate(response = c("y_it"), termlabels =
                             paste0("x_it_",c(0:(K-1))) %>% paste(collapse = " + "))
  result <- plm(formulate, data=data,
                effect = "individual",model="within")
  return(list(beta_hat = as.matrix(result$coefficients)))
}


#####  Interacctive Fixed Effect Methods  #####

### Least Squares Model ###
#Step 1:define funtion to caculate F_hat, dim of F_hat is (T_, r)
caculate_F_hat <- function(X_list, Y_list, beta_hat, r){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  K <- dim(X_list[[1]])[2]
  
  WWT <- matrix(0, nrow=T_, ncol=T_)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    W_i <- (Y_i - X_i %*% beta_hat)
    WWT <- WWT + W_i %*% t(W_i)
  }
  eig <- eigen(WWT)
  F_hat <- sqrt(T_)*eig$vectors[,order(eig$values, decreasing = T)[1:r]]
  return(F_hat)
}

#Step 2:define funtion to caculate Lambda_hat, dim of Lambda_hat is (r, N)
caculate_Lambda_hat <- function(X_list, Y_list, beta_hat, F_hat, r){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  K <- dim(X_list[[1]])[2]
  
  Lambda_hat <- matrix(NA, nrow = r, ncol = N)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    Lambda_hat[,i] <- t(F_hat) %*% (Y_i - X_i %*% beta_hat) / T_
  }
  return(Lambda_hat)
}

#Step 3:define funtion to caculate Beta_hat
caculate_beta_hat <- function(X_list, Y_list, F_, Lambda){
  N <- length(X_list)
  T_ <- dim(X_list[[1]])[1]
  K <- dim(X_list[[1]])[2]
  
  A <- matrix(0, nrow=K, ncol=K)
  B <- matrix(0, nrow=K, ncol=1)
  for(i in 1:N){
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    A <- A + t(X_i) %*% X_i
    B <- B + t(X_i) %*% (Y_i - F_ %*% Lambda[,i])
  }
  beta_hat <- solve(A) %*% B
  return(beta_hat)
}

#Step 5:caculate Beta_hat by iterations
#Input is the data {X_list, Y_list, df}, the tolerence level, and the number of factors r            
            
least_squares <- function(X_list, Y_list, df, tolerance, r){
  # Initialize
  K <- dim(X_list[[1]])[2]
  formulate <- reformulate(response = c("y_it"), termlabels =
                             paste0("x_it_",c(0:(K-1))) %>% paste(collapse = " + "))
  beta_hat_0 <- plm(formulate, data=df, model="pooling")$coefficients %>% as.matrix()
  
  beta_hat_list <- list(beta_hat_0)
  e <- Inf
  while (e > tolerance) {
    F_hat <- caculate_F_hat(X_list, Y_list, beta_hat_0, r)
    Lambda_hat <- caculate_Lambda_hat(X_list, Y_list, beta_hat_0, F_hat, r)
    beta_hat <- caculate_beta_hat(X_list, Y_list, F_hat, Lambda_hat)
    
    beta_hat_list[[length(beta_hat_list)+1]] <- beta_hat
    e <- norm(beta_hat - beta_hat_0, type = "F")
    beta_hat_0 <- beta_hat
  }
  return(list(beta_hat=beta_hat, beta_hat_list=beta_hat_list))
}

##########################
#####     Result     #####
##########################

#####  Simulation  #####
dgps <- list(DGP1=DGP1, DGP3=DGP3)
methods <- list(OLS_FE = function(X_list, Y_list, df) OLS_FE(df),  # make these methods have same params
                OLS_FE2 = function(X_list, Y_list, df) OLS_FE2(df),
                least_squares = 
                  function(X_list, Y_list, df)
                    least_squares(X_list, Y_list, df, tolerance=0.005), 4)

sim <- function(dgp, method, beta_true, all_N, all_T, nsims){
  # Initialize
  K <- length(beta_true)
  # Data frame to save every beta_hat
  df_beta_hat <- data.frame(T_ = rep(all_T, times = length(all_N), each = nsims),
                            N = rep(all_N, each = length(all_T) * nsims),
                            sim = rep(1:nsims, times = length(all_N) * length(all_T)),
                            beta = matrix(NA,ncol=K))
  # Loop over all_N and all_T and c(1:nsims) for simulation
  count_df_beta_hat <- 1
  for(i in 1:length(all_N)){
    N <- all_N[i]
    for(j in 1:length(all_T)){
      T_ <- all_T[j]
      for(h in 1:nsims){
        sim_data <- dgp(T_=T_, N=N, beta_true=beta_true)
        result <- method(sim_data$X_list, sim_data$Y_list, sim_data$df)
        beta_hat <- result$beta_hat
        df_beta_hat[count_df_beta_hat, 4:(3+K)] <- beta_hat
        count_df_beta_hat <- count_df_beta_hat + 1
      }
    }
  }
  return(list(df_beta_hat=df_beta_hat))
}

#####  Compute mean squared error  #####

#input: list of estimations and real parameters
#output: mean squared error
mse <- function(est_list, real_para){
  mse <- 0
  N <- length(est_list) 
  for (i in 1:N){
    mse <- mse + norm(est_list[[i]]-real_para, type="2")^2
  }
  mse <- 1/N * mse
  return(mse)
}

#####  Statistics  #####
statistics <- function(df_beta_hat, beta_true, all_N, all_T, nsims){
  # Initialize
  K <- length(beta_true)
  # Data frame to save statistics variable
  df_statistic <- data.frame(T_ = rep(all_T, times = length(all_N)),
                             N = rep(all_N, each = length(all_T)),
                             mean_x1 = NA,
                             var_x1 = NA,
                             bias_x1 = NA,
                             mse = NA)
  
  count_df_statistic <- 1
  for(i in 1:length(all_N)){
    N <- all_N[i]
    for(j in 1:length(all_T)){
      T_ <- all_T[j]
      row_range <- c((count_df_statistic*nsims-nsims+1) : (count_df_statistic*nsims)) # rows range of beta_hat for N and T_
      
      df_statistic$mean_x1[count_df_statistic] <- mean(df_beta_hat$beta.2[row_range])
      df_statistic$var_x1[count_df_statistic] <- var(df_beta_hat$beta.2[row_range])
      df_statistic$bias_x1[count_df_statistic] <- 
        abs(df_statistic$mean_x1[count_df_statistic] - beta_true[2]) / beta_true[2]
      beta_hat_list <- as.list(as.data.frame(t(df_beta_hat[row_range, 4:(3+K)])))
      df_statistic$mse[count_df_statistic] <- mse(beta_hat_list, beta_true)
      count_df_statistic <- count_df_statistic + 1
    }
  }
  return(list(df_statistic = df_statistic))
}

#####  Test  #####
dgp2 <- DGP2(10, 100, beta_true = c(1,2))
OLS_FE2(dgp2$df)

dgp1 <- DGP1(T_ = 3,
             N = 400,  beta_true = c(1,2,3))
OLS_FE(dgp1$df)$beta_hat
OLS_FE2(dgp1$df)$beta_hat
all.equal(OLS_FE(dgp1$df)$beta_hat, OLS_FE2(dgp1$df)$beta_hat)

                
              
T_ <- 1000
N <- 1000
tol <-0.005
beta_true <- c(5,1,3)
dgp3 <- DGP3(T_,N,beta_true)
df <-  dgp3$df
X_list <- dgp3$X_list
Y_list <- dgp3$Y_list

ls <- least_squares(X_list, Y_list, df, tol)
ls$beta_hat
ls$beta_hat_list
ls_mse <- mse(ls$beta_hat_list, c(5,1,3))

##### Generate table #####
all_N <- c(100)
all_T <- c(10,20,50,100)
nsims <- 100
beta_true <- c(5,1,3)
sim1 <- sim(dgps$DGP3, methods$least_squares, beta_true, all_N, all_T, nsims)
stat1 <- statistics(sim1$df_beta_hat, beta_true, all_N, all_T, nsims)
stat1$df_statistic

###########################
####   Visualization   ####
###########################

##### sim1 #####
all_n <- c(1000, 2000, 3000)
nsims <- 100
beta_true <- c(1:3)
sim1 <- sim(K = 3, beta_true = beta_true, all_n = all_n, all_T = c(3), nsims = nsims)
beta_hat_1 <- sim1$beta_hat_1
beta_true_1 <- beta_true[1]
beta_hat_1_df <- data.frame(t(beta_hat_1[,1,]))
colnames(beta_hat_1_df) <- all_n
beta_hat_1_df <- melt(beta_hat_1_df, measure.vars=colnames(beta_hat_1_df), variable.name = "N")
# generate point and box plot for every N
plot.heigth <- max(beta_true_1 - min(beta_hat_1_df$value), max(beta_hat_1_df$value) - beta_true_1)
point_plot <- ggplot(data = beta_hat_1_df) +
  geom_point(aes(y = value, x = c(1:(nsims*length(all_n))), color = N, shape = N)) +
  geom_hline(yintercept = beta_true_1, color = I("black")) +
  ylim(beta_true_1 - plot.heigth, beta_true_1 + plot.heigth) +
  labs(title = "point plot for every N", x = "iteration", y = "beta_hat_1") +
  theme(axis.text.x = element_blank())
print(point_plot)

box_plot <- ggplot(data = beta_hat_1_df) +
  geom_boxplot(aes(y = value, x = N, color = N)) +
  ylim(min(beta_hat_1_df$value), max(beta_hat_1_df$value)) +
  labs(title = "box plot for every N", y = "beta_hat_1")
print(box_plot)

##### sim2 #####
sim2 <- sim(K = 3, beta_true = c(1:3), all_n = seq(1000, 10000, 1000), all_T = seq(3,30,3), nsims = 1)
beta_hat_1 <- sim2$beta_hat_1
beta_true_1 <- 1
summary_beta_hat_1 <- sim2$summary_beta_hat_1

im <- with(summary_beta_hat_1, interp(T_,N,bias))
with(im,image(x,y,z, xlab = "T_", ylab = "N"))

plot_ly(x=im$x, y=im$y, z=im$z, type="surface") %>% add_surface() %>% layout(scene = list(
  xaxis = list(title = 'T'),
  yaxis = list(title = 'N'),
  zaxis = list(title = 'bias')))
