library(Matrix)
library(rgl)
library(glmnet)


simulation <- function(init_num,k_mat,duration,variation=0){
  result <- matrix(t(init_num))
  for(i in 1:duration){
    if(is.list(k_mat)){
      result <- cbind(result,k_mat[[i]]%*%result[,i])
    }else{
      result <- cbind(result,k_mat%*%result[,i])
    }
  }
  colnames(result) <- 1:ncol(result)
  rownames(result) <- c('G1','S','G2','Other')
  result <- result+matrix(rnorm(nrow(result)*ncol(result),sd = variation),nrow = nrow(result),ncol = ncol(result))
  return(result)
}

solve_mat <- function(result1,model_prior=NULL,w=NULL){
  matrix_all <- matrix(0,nrow = nrow(result1)-1,ncol = nrow(result1)^2)
  ref_ind <- (0:(nrow(result1)-1))*nrow(result1)
  for(i in 1:nrow(matrix_all)){
    matrix_all[i,i+ref_ind] <- 1
  }
  matrix_all[3,3] <- 1/2
  vector_all <- c(rep(1,nrow(matrix_all)),as.numeric(result1[,-1]))
  colnames(matrix_all) <- letters[1:ncol(matrix_all)]
  
  for(i in 2:ncol(result1)-1){
    matrix_x <- matrix(0,nrow = nrow(result1),ncol = nrow(result1)^2)
    for(j in 1:nrow(matrix_x)){
      matrix_x[j,(4*j-3):(4*j)] <- result1[,i]
    }
    matrix_all <- rbind(matrix_all,matrix_x)
  }
  
  #prior based on cell cycling model
  if(is.null(model_prior)){
    model_prior <- rep(0,ncol(matrix_all))
    model_prior[c(2,7,9)] <- 1  #unlikely transition rates
    model_prior[c(1,6,11,16)] <- 0 #stable rates
    model_prior[c(3,5,10,4,8,12,13,14,16)] <- 0 #likely transition rates(cell cycle and other)

  }
  if(is.null(w)){
    w2 <- c(rep(1,3),rep(1,nrow(matrix_all)-3))
  }else{
    w2 <- c(rep(1,3),rep(w,each=4))
  }
  hbound <- rep(1,ncol(matrix_all))
  hbound[3] <- 2
  glm_fit <- glmnet(matrix_all,vector_all,penalty.factor = model_prior,lower.limits = 0,upper.limits = hbound,weights = w2,
                    intercept = F,exclude = c(2,7,9),standardize=F)
  result_mat <- matrix(as.matrix(glm_fit$beta)[,which.max(glm_fit$dev.ratio)],nrow(result1),nrow(result1),byrow = T)
}




get_k_mat <- function(exp_data,u=3){
  ind <- as.numeric(colnames(exp_data)[-ncol(exp_data)])
  result_mat_list <- list()
  previous <- NULL
  for(i in 2:ncol(exp_data)-1){
    w <- exp(-(1/(2*u^2))*(ind-i)^2)
    result_mat_list[[i]] <- solve_mat(exp_data,w = w,previous = previous)
    previous <- as.numeric(t(result_mat_list[[i]]))
  }
  return(result_mat_list)
}

#extract G1-S,S-G2,G2-M rate:
plot_rate <- function(result_mat_list,spar=0.35,ind_list =list(c(1,3),c(2,1),c(3,2)),main='Cell cycle transitions'){
  plot(1:length(result_mat_list),type="n",ylim = c(0,1),xlab = 'Time',ylab = 'Rate',main = main)
  for(i in 1:length(ind_list)){
    rate_data <- unlist(lapply(result_mat_list,function(x)x[ind_list[[i]][1],ind_list[[i]][2]]))
    if(min(ind_list[[i]] == c(1,3))==1) rate_data <- rate_data/2
    lines(smooth.spline(1:length(result_mat_list), rate_data, spar = spar),col=i)
  }
}

plot_count <- function(exp_data){
  plot(1:ncol(exp_data),apply(exp_data,2,max),xlab = 'Time',ylim = c(0,max(apply(exp_data,2,max))),ylab = 'Cell count in each phase',type = 'n')
  lines(1:ncol(exp_data),exp_data[1,],col=3)
  lines(1:ncol(exp_data),exp_data[2,],col=1)
  lines(1:ncol(exp_data),exp_data[3,],col=2)
  lines(1:ncol(exp_data),exp_data[4,],col=4)
  legend('topleft',rownames(exp_data),lty=1,col = c(3,1,2,4))
}

lsmooth <- function(x,span){
  l_model <- loess(count ~ time, span=span, data.frame(time=1:length(x), count=x))
  return(predict(l_model,newdata = 1:length(x)))
}

data_smooth <- function(d,span=0.3){
  dsmooth <- t(apply(d,1,function(x)lsmooth(x,span = span)))
  colnames(dsmooth) <- 1:ncol(dsmooth)
  return(dsmooth)
}
