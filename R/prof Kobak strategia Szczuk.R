Strategia <- function(mi, sigma, t, S0, r, q, K, nsym, pbankructwa){
  start_time <- Sys.time()
  y <- optim(par = 1, fn = function(x) 
    (abs(Symulacja(mi,sigma,t,S0,r,q,K,nsym,x)[1] - pbankructwa)), method = "Brent",
    lower = -0.999, upper = 10^2)
  end_time <- Sys.time()
  print(paste("Up³ynê³o: ", round(end_time - start_time, 2), " sekund"))
  return(y)
}


Symulacja <- function(mi, sigma, t, S0, r, q, K, nsym, stawkaL){
  #start_time <- Sys.time()
  Wt <- replicate(nsym, vector())
  for (j in 1:nsym){
    Wt[[j]] <- rnorm(t-1)
  }
  S <- replicate(nsym, vector())
  Bogactwo <- replicate(nsym, vector())
  for(j in 1:nsym){
    S[[j]][[1]] <- S0
    Bogactwo[[j]][[1]] <- K*S0
    for(i in (2:t)){
      S[[j]][[i]] <- S[[j]][[i-1]]*exp((mi-sigma^2/2)*(1/365) + sigma*Wt[[j]][[i-1]]*sqrt(1/365))
      Bogactwo[[j]][[i]] <- (1+stawkaL)*K*S0*(1+r)^((i-1)/365)-stawkaL*K*S[[j]][[i]]*exp(q*(i-1)/365)
    }
  }
  Bogactwobezbank <- Bogactwo[which(lapply(Bogactwo, function(x) any(x < 0)) == FALSE)]
  bankructwa <- sum(unlist(lapply(S, function(x) 
    sum(any((head(x,1)/x)*(1+r)^((0:(t-1))/365) <= stawkaL/(1+stawkaL))))))
  #end_time <- Sys.time()
  #print(paste("Up³ynê³o: ", round(end_time - start_time, 2), " sekund"))
  #print(paste("Procent bankructw to: ", bankructwa/nsym))
  #print(paste("Œrednie bogactwo to: ", round(mean(unlist(lapply(Bogactwo, function(x) tail(x,1)))),3)))
  #print(paste("Œrednie bogactwo bez bankructwa to: ",
  #            round(mean(unlist(lapply(Bogactwobezbank, function(x) tail(x,1)))),3)))
  return(c(Procentbankructw = bankructwa/nsym, 
           Œredniebogactwo = round(mean(unlist(lapply(Bogactwo, function(x) tail(x,1)))),3), 
           Œredniebogactwobezbankructwa = round(mean(unlist(lapply(Bogactwobezbank, function(x) tail(x,1)))),3)))
}

Symulacja(-0.1,0.3,1800,4,0.05,0.02,1,1000,1.107751)

Strategia(-0.1,0.3,1800,4,0.05,0.02,1,1000,0.15)




