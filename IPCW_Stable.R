IPCW_Stable<- function(Gg){
  G<-Gg
  G<-G[order(G[,5]),]
  n<-length(G[,2])
  A<- G[,1]
  U<- G[,5]#?۲⵽??ʱ??
  Z<- G[,10]#Э??��
  Ga<-G[,11]#??׼Э??��
  delta1<- G[,6]#I(T<=pmin(C1,C2))
  delta2<- G[,7]#I(C2<=pmin(C1,T))
  X_t<- G[,8] 
  W_t<- G[,9] 
  
  
  U_matrix <- matrix(1,n,1)%*%t(as.matrix(U))
  Y_t <- matrix(-1,n,n) #matrix
  UU_matrix1 <- as.matrix(U)%*%matrix(1,1,n)#
  Y_t <- I(UU_matrix1>=U_matrix)*1
  
  Xt0_matrix <- as.matrix(X_t)%*%matrix(1,1,n)
  Wt0_matrix <- as.matrix(W_t)%*%matrix(1,1,n)
  X_tt <- (Xt0_matrix >=U_matrix)*1##????ɾʧЭ??��X??t?? ###X_t1<-Xt0_matrix 
  W_tt <- (Wt0_matrix >=U_matrix)*1###????ɾʧЭ??��W??t??
  #W_t1 <- as.matrix(W_t)%*%matrix(1,1,n0)
  
  dNC_t<-(I(UU_matrix1 == U_matrix)*1)
  dNC_t[which(delta2!=1),]=0
  
  U_0<-c(0,U[-n])
  Time_clculas_M<-matrix(1,n,1)%*%t(as.matrix(U-U_0))  
  
  Indix1<- which(A==1)
  Indix0<- which(A==0)
  Col1<- c(1:which((A*U)==max(A*U)))
  Col0<- c(1:which(((1-A)*U)==max((1-A)*U)))
  
  n1<-sum(A==1)
  n0<-sum(A==0)
  
  dNC_t1<- dNC_t[Indix1,Col1]
  dNC_t0<- dNC_t[Indix0,Col0]
  
  Time_clculas_M1<-Time_clculas_M[Indix1,Col1]
  Time_clculas_M0<-Time_clculas_M[Indix0,Col0]
  
  Y_t1<- Y_t[Indix1,Col1]
  Y_t0<- Y_t[Indix0,Col0]
  
  X_t1 <- X_tt[Indix1,Col1]##????ɾʧЭ??��X??t?? ###X_t1<-Xt0_matrix 
  X_t0 <- X_tt[Indix0,Col0]
  
  W_t1 <- W_tt[Indix1,Col1]
  W_t0 <- W_tt[Indix0,Col0]
  
  
  fun_C1<-function(theta){
    ###---------------------?????��Ʒ???-------------------
    Wn_t<-W_t1/exp(theta[2]*X_t1)                     ###????\tilde(W)^C_ij(t)??һά??
    Xn_t<-X_t1                                        ###????\tilde(X)^C_ij(t)?ڶ?ά??
    
    Wn_bar_t<-colSums(Y_t1*W_t1)/colSums(Y_t1*exp(theta[2]*X_t1))                     ####????\bar(W)_j(t)??һά??
    Xn_bar_t<-colSums(Y_t1*exp(theta[2]*X_t1)*X_t1)/colSums(Y_t1*exp(theta[2]*X_t1))  ####????\bar(X)_j(t)?ڶ?ά??
    
    fC_Wn_t<-Wn_t-matrix(1,n1,1)%*%t(as.matrix(Wn_bar_t))   ###?????��Ʒ??̵ı???????D_ij(t)-bar(D)_j(t)
    fC_Xn_t<-Xn_t-matrix(1,n1,1)%*%t(as.matrix(Xn_bar_t)) 
    
    S_beta_inte1<- sum(fC_Wn_t*dNC_t1)###S_beta(\tau;beta)ǰ??????
    S_beta_inte2<- sum( (Y_t1*fC_Wn_t*(theta[1]*W_t1))%*%((U-U_0)[Col1] )  )###S_beta(\tau;beta)????????
    S_beta_inte<- S_beta_inte1-S_beta_inte2  ###S_beta(\tau;beta)  ###S_beta(\tau;beta)
    
    S_gamma_inte1<- sum(fC_Xn_t*dNC_t1)
    S_gamma_inte2<- sum( (Y_t1*fC_Xn_t*theta[1]*W_t1)%*%((U-U_0)[Col1]  ))
    S_gamma_inte<-S_gamma_inte1-S_gamma_inte2    ###S_gamma(\tau;beta)
    
    ### S_theta<-matrix(c(S_beta_inte,S_gamma_inte),2,1)####?��Ʒ???
    ###-------------------?��Ʒ??̣?��??ά?ȣ?????????---------------   
    c(S_beta_inte^2+S_gamma_inte^2)
  }
  
  result_C1=optim(c(0.3,1.5),fn=fun_C1)
  result_C1
  es_theta1<-result_C1$par
  
  #theta=c(0.3,2.5)
  fun_C0<-function(theta){
    ###---------------------?????��Ʒ???-------------------
    Wn_t<-W_t0/exp(theta[2]*X_t0)                     ###????\tilde(W)^C_ij(t)??һά??
    Xn_t<-X_t0                                        ###????\tilde(X)^C_ij(t)?ڶ?ά??
    
    Wn_bar_t<-colSums(Y_t0*W_t0)/colSums(Y_t0*exp(theta[2]*X_t0))                     ####????\bar(W)_j(t)??һά??
    Xn_bar_t<-colSums(Y_t0*exp(theta[2]*X_t0)*X_t0)/colSums(Y_t0*exp(theta[2]*X_t0))  ####????\bar(X)_j(t)?ڶ?ά??
    
    fC_Wn_t<-Wn_t-matrix(1,n0,1)%*%t(as.matrix(Wn_bar_t))   ###?????��Ʒ??̵ı???????D_ij(t)-bar(D)_j(t)
    fC_Xn_t<-Xn_t-matrix(1,n0,1)%*%t(as.matrix(Xn_bar_t)) 
    
    S_beta_inte1<- sum(fC_Wn_t*dNC_t0)###S_beta(\tau;beta)ǰ??????
    S_beta_inte2<- sum( (Y_t0*fC_Wn_t*(theta[1]*W_t0))%*%((U-U_0)[Col0] )  )###S_beta(\tau;beta)????????
    S_beta_inte<- S_beta_inte1-S_beta_inte2  ###S_beta(\tau;beta)  ###S_beta(\tau;beta)
    
    S_gamma_inte1<- sum(fC_Xn_t*dNC_t0)
    S_gamma_inte2<- sum( (Y_t0*fC_Xn_t*theta[1]*W_t0)%*%((U-U_0)[Col0]  ))
    S_gamma_inte<-S_gamma_inte1-S_gamma_inte2    ###S_gamma(\tau;beta)
    
    ### S_theta<-matrix(c(S_beta_inte,S_gamma_inte),2,1)####?��Ʒ???
    ###-------------------?��Ʒ??̣?��??ά?ȣ?????????---------------   
    c(S_beta_inte^2+S_gamma_inte^2)
  }
  
  result_C0=optim(c(0.3,2.5),fn=fun_C0)
  result_C0
  es_theta0<-result_C0$par
  
  
  ###----------------????ɾʧ???ջع?ϵ???��ƽ???-------------------------
  
  
  ###-----------------????ɾʧ??׼???պ????��?---------------------------
  C_hazardfun<- function(Y_t,X_t,W_t,dNC_t,es_theta,n,Time_clculas_M){
    r_C_j_t<-   1/colSums(Y_t*exp( es_theta[2]*X_t )  ) 
    r_C_j_t_matrix<-   matrix(1,n,1)%*%t(as.matrix(r_C_j_t))
    
    C_baseline_0<- (colSums(dNC_t)-colSums(Y_t*es_theta[1]*W_t*Time_clculas_M))*r_C_j_t
    #C_baseline_0[which(C_baseline_0<0)]=0
    #C_baseline<-cumsum(C_baseline_0)
    return(C_baseline_0)}  
  ###plot(U[1:ncol(Y_t1)],C_baseline1)
  C_baseline_1<- C_hazardfun(Y_t1,X_t1,W_t1,dNC_t1,es_theta1,n1,Time_clculas_M1)#differential calculus
  C_baseline_0<- C_hazardfun(Y_t0,X_t0,W_t0,dNC_t0,es_theta0,n0,Time_clculas_M0)
  
  C_baseline1<-cumsum(C_baseline_1) #cumulative hazard function value
  C_baseline0<-cumsum(C_baseline_0)
  
  C_baseline_1_M<-matrix(1,n,1)%*%t(as.matrix(C_baseline_1)) 
  C_baseline_0_M<-matrix(1,n,1)%*%t(as.matrix(C_baseline_0))  
  
  C_baseline1_M<-matrix(1,n,1)%*%t(as.matrix(C_baseline1)) 
  C_baseline0_M<-matrix(1,n,1)%*%t(as.matrix(C_baseline0))  
  ###-----------------????ɾʧ??׼???պ????��ƽ???-------------------------
  
  ###--------------------?��?ÿ???????ķ??պ???-----------------------------
  A1<-matrix(1,ncol(Y_t1),ncol(Y_t1))
  A1[lower.tri(A1)]=0
  A0<-matrix(1,ncol(Y_t0),ncol(Y_t0))
  A0[lower.tri(A0)]=0
  
  C_hazard1<-  (exp(es_theta1[2]*X_tt[,Col1])* C_baseline_1_M)%*%A1+(es_theta1[1]*W_tt[,Col1]*Time_clculas_M[,Col1])%*%A1
  C_hazard0<-  (exp(es_theta0[2]*X_tt[,Col0])* C_baseline_0_M)%*%A0+(es_theta0[1]*W_tt[,Col0]*Time_clculas_M[,Col0])%*%A0
  
  ###--------------------ÿ???????ķ??պ????��ƽ???---------------------------
  
  ###--------------------------Ȩ?????��?----------------------------------
  ##rho_t<-exp(C_hazard1)[1,]???Ƚ?Ȩ??
  rho_t1<-exp(C_hazard1[Indix1,])/exp( C_baseline1_M[Indix1,])
  rho_t0<-exp(C_hazard0[Indix0,])/exp( C_baseline0_M[Indix0,])
  rho_t1[rho_t1>10]<-10
  rho_t0[rho_t0>10]<-10
  max(rho_t1)
  max(rho_t0)
  ###-------------------------Ȩ?????��ƽ???-----------------------------
  
  ###----------------------?��?????ģ?ͷ??ջع?ϵ??---------------------------
  
  dN_t<-(I(UU_matrix1 == U_matrix)*1)
  dN_t[which(delta1!=1),]=0
  
  Z_M<-  as.matrix(Z)%*% matrix(1,1,n)
  Ga_M<-   as.matrix(Ga)%*% matrix(1,1,n)
  
  dN_t1<- dN_t[Indix1,Col1]
  dN_t0<- dN_t[Indix0,Col0]
  
  Z_M1<-  Z_M[Indix1,Col1]
  Z_M0<-  Z_M[Indix0,Col0]
  Ga_M1<- Ga_M[Indix1,Col1]
  Ga_M0<- Ga_M[Indix0,Col0]
  
  #beta<-c(0.3,1.2)
  funT1<- function(beta){
    ###---------------?????��Ʒ???-------------------------
    Gan_t<- Ga_M1/exp(beta[2]*Z_M1)
    Zn_t<-  Z_M1
    
    Gan_bar_t<- colSums(rho_t1*Y_t1*Ga_M1)/colSums( rho_t1*Y_t1*exp(beta[2]*Z_M1)   )
    Zn_bar_t<- colSums(rho_t1*Y_t1*Z_M1*exp(beta[2]*Z_M1))/colSums(rho_t1*Y_t1*exp(beta[2]*Z_M1) )
    
    fT_Gan_t<-  Gan_t-matrix(1,n1,1)%*%t(as.matrix(Gan_bar_t))
    fT_Zn_t<-  Zn_t-matrix(1,n1,1)%*%t(as.matrix(Zn_bar_t))
    
    S_alpha_inte1<-  sum(rho_t1*fT_Gan_t*dN_t1)  
    S_alpha_inte2<-  sum((rho_t1*fT_Gan_t*Y_t1*beta[1]*Ga_M1)%*%((U-U_0)[Col1]))
    
    S_alpha<-  S_alpha_inte1- S_alpha_inte2
    
    S_eta_inte1<- sum(rho_t1*fT_Zn_t*dN_t1)  
    S_eta_inte2<-  sum((rho_t1*fT_Zn_t*Y_t1*beta[1]*Ga_M1)%*%((U-U_0)[Col1]))
    
    S_eta<-  S_eta_inte1- S_eta_inte2
    
    ###----------------?��Ʒ??̽???----------------------
    c(S_alpha^2+S_eta^2)
  }
  
  funT0<- function(beta){
    ###---------------?????��Ʒ???-------------------------
    Gan_t<- Ga_M0/exp(beta[2]*Z_M0)
    Zn_t<-  Z_M0
    
    Gan_bar_t<- colSums(rho_t0*Y_t0*Ga_M0)/colSums( rho_t0*Y_t0*exp(beta[2]*Z_M0)   )
    Zn_bar_t<- colSums(rho_t0*Y_t0*Z_M0*exp(beta[2]*Z_M0))/colSums(rho_t0*Y_t0*exp(beta[2]*Z_M0) )
    
    fT_Gan_t<-  Gan_t-matrix(1,n0,1)%*%t(as.matrix(Gan_bar_t))
    fT_Zn_t<-  Zn_t-matrix(1,n0,1)%*%t(as.matrix(Zn_bar_t))
    
    S_alpha_inte1<-  sum(rho_t0*fT_Gan_t*dN_t0)  
    S_alpha_inte2<-  sum((rho_t0*fT_Gan_t*Y_t0*beta[1]*Ga_M0)%*%((U-U_0)[Col0]))
    
    S_alpha<-  S_alpha_inte1- S_alpha_inte2
    
    S_eta_inte1<- sum(rho_t0*fT_Zn_t*dN_t0)  
    S_eta_inte2<-  sum((rho_t0*fT_Zn_t*Y_t0*beta[1]*Ga_M0)%*%((U-U_0)[Col0]))
    
    S_eta<-  S_eta_inte1- S_eta_inte2
    
    ###----------------?��Ʒ??̽???----------------------
    c(S_alpha^2+S_eta^2)
  }
  
  
  result_T1=optim(c(0.15,-0.3),funT1)
  result_T0=optim(c(0.15,-0.4),funT0)
  #result_T<- multiroot(funT,c(0.15,-1),1e-8,maxiter=10000)
  result_T1
  result_T0
  es_beta1<-result_T1$par
  es_beta0<-result_T0$par
  
  ###-------------------------beta?????��?------------------
  ###------------------- A_T????????????---------------   
  AA_T<- function(es_beta,Ga_M,Z_M,Y_t,rho_t,n,dN_t,Col){
    Gan_t<- Ga_M/exp(es_beta[2]*Z_M)
    Zn_t<-  Z_M
    
    Gan_bar_t<- colSums(rho_t*Y_t*Ga_M)/colSums( rho_t*Y_t*exp(es_beta[2]*Z_M)   )
    Zn_bar_t<- colSums(rho_t*Y_t*Z_M*exp(es_beta[2]*Z_M))/colSums(rho_t*Y_t*exp(es_beta[2]*Z_M) )
    
    fT_Gan_t<-  Gan_t-matrix(1,n,1)%*%t(as.matrix(Gan_bar_t))
    fT_Zn_t<-  Zn_t-matrix(1,n,1)%*%t(as.matrix(Zn_bar_t))
    ###------???????????ջع?ϵ??Э??????-----
    S_1_alpha_alpha<-  sum(    (rho_t*Y_t*exp(es_beta[2]*Z_M)*fT_Gan_t*fT_Gan_t)%*%((U-U_0)[Col])   )###?????��Ʒ??̵?΢????ʽ
    
    r_j_t<-   colSums( rho_t*Y_t*exp(es_beta[2]*Z_M )  )       ###?????????  r^C_j(t)?ĵ???
    Yn_3<- colSums(   Y_t*exp(es_beta[2]*Z_M)* fT_Gan_t* fT_Zn_t*rho_t) /r_j_t
    Yn_3_matrix<-  matrix(1,n,1)%*%t(as.matrix(Yn_3)) 
    
    Yn_4<- colSums(rho_t*Y_t*exp(es_beta[2]*Z_M)*fT_Zn_t*fT_Zn_t)/r_j_t
    Yn_4_matrix<-  matrix(1,n,1)%*%t(as.matrix(Yn_4)) 
    
    S_1_eta_eta<- sum(Yn_4_matrix*dN_t*rho_t)-sum( (rho_t*Y_t*Yn_4_matrix*es_beta[1]*Ga_M)%*%((U-U_0)[Col] ))        ###?????��Ʒ??̵?΢????ʽ
    
    S_1_alpha_eta<- sum(Yn_3_matrix*dN_t*rho_t)-sum( (rho_t*Y_t*Yn_3_matrix*es_beta[1]*Ga_M)%*%((U-U_0)[Col] ))
    
    S_1_eta_alpha<- sum(    (rho_t*Y_t*exp(es_beta[2]*Z_M)*fT_Zn_t*fT_Gan_t)%*%((U-U_0)[Col] ))###?????��Ʒ??̵?΢????ʽ
    
    ###------------------A_T??????------------------------
    A_T<-matrix(c(S_1_alpha_alpha,S_1_eta_alpha,S_1_alpha_eta,S_1_eta_eta),2,2)
    
    ###------------------------?��Ʒ??̼???-----------------------
    S_alpha<-rowSums(rho_t*fT_Gan_t*dN_t)-(rho_t*fT_Gan_t*Y_t*es_beta[1]*Ga_M)%*%((U-U_0)[Col] )##?????��Ʒ???
    S_eta<-rowSums(rho_t*fT_Zn_t*dN_t)-(rho_t*Y_t*es_beta[1]*Ga_M*fT_Zn_t)%*%((U-U_0)[Col] ) ##?????��Ʒ???
    ###------------------------?��Ʒ??̽???-----------------------
    V_T<- matrix(c(sum(S_alpha^2),sum(S_alpha*S_eta),sum(S_eta*S_alpha),sum(S_eta^2)),2,2)
    Sigma_T<- solve(A_T)%*%V_T%*%t(solve(A_T))
    list(Sigma_T,A_T,S_alpha,S_eta)
  }
  Sigma_T1<- AA_T(es_beta1,Ga_M1,Z_M1,Y_t1,rho_t1,n1,dN_t1,Col1)
  #es_beta=es_beta1;Ga_M=Ga_M1;Z_M=Z_M1;Y_t=Y_t1;rho_t=rho_t1;n=n1;dN_t=dN_t1;Col=Col1;
  Sigma_T0<- AA_T(es_beta0,Ga_M0,Z_M0,Y_t0,rho_t0,n0,dN_t0,Col0)
  
  ###-----------------------beta??Э?????????��ƽ???--------------------
  
  ###-------------------------?��?T??׼???պ???-----------------
  T_hazardfun<- function(es_beta,rho_t,Y_t,Z_M,n,Ga_M,dN_t,Time_clculas_M){
    r_j_t<-   1/colSums(  rho_t*Y_t*exp( es_beta[2]*Z_M )  ) 
    r_j_t_matrix<-   matrix(1,n,1)%*%t(as.matrix(r_j_t))
    
    T_baseline_0<- (colSums(rho_t*dN_t)-colSums(rho_t*Y_t*es_beta[1]*Ga_M*Time_clculas_M))*r_j_t
    return(T_baseline_0)
  }
  
  T_baseline_1<- T_hazardfun(es_beta1,rho_t1,Y_t1,Z_M1,n1,Ga_M1,dN_t1,Time_clculas_M1)
  T_baseline_0<- T_hazardfun(es_beta0,rho_t0,Y_t0,Z_M0,n0,Ga_M0,dN_t0,Time_clculas_M0)
  
  T_baseline1<-cumsum(T_baseline_1)  ###?ۼƻ?׼???պ???
  T_baseline0<-cumsum(T_baseline_0)  ###?ۼƻ?׼???պ???
  
  ###plot(U[1:ncol(Y_t1)],T_baseline1)
  ###par(mfrow=c(1,2))
  
  T_baseline_1_M<-matrix(1,n,1)%*%t(as.matrix(T_baseline_1)) 
  T_baseline_0_M<-matrix(1,n,1)%*%t(as.matrix(T_baseline_0)) 
  
  T_baseline1_M<-matrix(1,n,1)%*%t(as.matrix(T_baseline1))  
  T_baseline0_M<-matrix(1,n,1)%*%t(as.matrix(T_baseline0))  
  
  ###------------------------T??׼???պ????��ƽ???-------------------
  
  ###--------------------?��?ÿ???????????????պ???----------------------------
  Zz_M1<- Z_M[,Col1]
  Gaga_M1<- Ga_M[,Col1]
  
  Zz_M0<- Z_M[,Col0]
  Gaga_M0<- Ga_M[,Col0]
  
  Time_clculas_Mm1<-Time_clculas_M[,Col1]
  Time_clculas_Mm0<-Time_clculas_M[,Col0]
  
  T_hazard1<-  (exp(es_beta1[2]*Zz_M1)* T_baseline_1_M)%*%A1+(es_beta1[1]*Gaga_M1*Time_clculas_Mm1)%*%A1
  T_hazard0<-  (exp(es_beta0[2]*Zz_M0)* T_baseline_0_M)%*%A0+(es_beta0[1]*Gaga_M0*Time_clculas_Mm0)%*%A0
  
  ###--------------------ÿ???????ķ??պ????��ƽ???---------------------------
  S_T1<- exp(-T_hazard1)###????????
  S_T0<- exp(-T_hazard0)###????????
  
  ###-----------------------------?޶???ֵ?????Ĺ��?----------------------------------
  Limit_time<- 10
  L<- Limit_time
  Time_up<- U[1:min(ncol(Y_t1),ncol(Y_t0))]
  if(L<max(U)){
    if(length(which(U==L))==0){
      nl <- min(which(U>L))-1
      DN_1_L_M<- cbind(dN_t[Indix1,1:nl],matrix(0,n1,1))
      DN_0_L_M<- cbind(dN_t[Indix0,1:nl],matrix(0,n0,1))
    }else{nl <- min(which(U>=L))-1
    DN_1_L_M<- cbind(dN_t[Indix1,1:nl],dN_t[Indix1,nl+1])      
    DN_0_L_M<- cbind(dN_t[Indix0,1:nl],dN_t[Indix0,nl+1])}
    if(nl==0){
      S_T_L1<-S_T1[,1]
      S_T_L0<-S_T0[,1]
      Y_t1_L<- Y_t1[,1]
      Y_t0_L<- Y_t0[,1]
      L_time<- c(L)
      L_time_0 <- c(0)
      S_1<-colSums(S_T_L1)/n
      S_0<-colSums(S_T_L0)/n
      M_time1 <- sum(S_1*(L))
      M_time0 <- sum(S_0*(L))
    }else{
      S_T_L1<-cbind(S_T1[,1:nl],S_T1[,nl])
      S_T_L0<-cbind(S_T0[,1:nl],S_T0[,nl])
      Y_t1_L<- cbind(Y_t[Indix1,1:nl],Y_t[Indix1,nl+1]) 
      Y_t0_L<- cbind(Y_t[Indix0,1:nl],Y_t[Indix0,nl+1]) 
      S_1<-colSums(S_T_L1)/n
      S_0<-colSums(S_T_L0)/n
      L_time<- c(U[1:nl],L)
      L_time_0 <- c(0,L_time[1:nl])
      M_time1 <- sum(S_1*(L_time-L_time_0))
      M_time0 <- sum(S_0*(L_time-L_time_0))}
  }else{
    nl<-n
    DN_1_L_M<- cbind(dN_t[Indix1,1:nl],rep(0,n1))      
    DN_0_L_M<- cbind(dN_t[Indix0,1:nl],rep(0,n0))
    S_T_L1<-cbind(S_T1[,1:nl],S_T1[,nl])
    S_T_L0<-cbind(S_T0[,1:nl],S_T0[,nl])
    Y_t1_L<- cbind(Y_t[Indix1,1:nl],rep(0,n1))
    Y_t0_L<- cbind(Y_t[Indix0,1:nl],rep(0,n0)) 
    
    S_1<-colSums(S_T_L1)/n
    S_0<-colSums(S_T_L0)/n
    L_time <- c(U,L)
    L_time_0 <- c(0,U)
    M_time1 <- sum(S_1*(L_time-L_time_0))
    M_time0 <- sum(S_0*(L_time-L_time_0))}
  

  ###-----------------------------??ֵ?????��ƽ???---------------------------------
  
  ###---------------------------delta?????��?-----------------------
  mu_L1<-  S_T_L1%*%as.matrix(L_time-L_time_0) #Ϊһ??ʱ??LΪֹ??Sij??t???????Ļ???
  mu_L_M1<-  as.matrix(mu_L1)%*%matrix(1,1,nl+1) 
  
  Time_clculas_L_M<-matrix(1,n,1)%*%t(as.matrix(L_time-L_time_0))  
  #   DN_0_L_M<- cbind(DN_t0[,1:nl],matrix(0,n0,1)) ####???????????????
  
  AA1<-matrix(1,nl+1,nl+1)
  AA1[lower.tri(AA1)]=0
  Z_M_L<-  as.matrix(Z)%*% matrix(1,1,nl+1)
  Ga_M_L<-   as.matrix(Ga)%*% matrix(1,1,nl+1)
  ####   Y_t1_L<- cbind(Y_t1[,1:nl],Y_t1[,nl])  ####???????????????
  rho_t_L1<- cbind(rho_t1[,1:nl],rho_t1[,nl])####???????????????
  
  mu_t_L<-  (S_T_L1*Time_clculas_L_M)%*%AA1  ###????u_ij(t)
  
  Yn_5<- (colSums( (mu_L_M1- mu_t_L)*exp(es_beta1[2]*Z_M_L)*Z_M_L)/n)/colSums(Y_t1_L*rho_t_L1*exp(es_beta1[2]*Z_M_L[Indix1,]))
  ####   Yn_5_M<- matrix(1,n0,1)%*%t(as.matrix(Yn_5))
  
  delta_eta_1<-  sum(Yn_5*colSums(rho_t_L1*DN_1_L_M))
  delta_eta_2<-  sum(Yn_5*colSums(Y_t1_L*rho_t_L1*es_beta1[1]*Ga_M_L[Indix1,]*Time_clculas_L_M[Indix1,])  )
  delta_eta_inte1<- -delta_eta_1+ delta_eta_2 #-sum((S_T_L0*T_baseline1_M[,1:(nl+1)]*exp(es_beta0[2]*Z_M_L)*Z_M_L)%*%(L_time-L_time_0))/n
  
  delta_alpha_inte1<- -sum( colSums(S_T_L1* Ga_M_L)*L_time*(L_time-L_time_0))/n
  ###  delta_eta_inte1<-     -sum(colSums(S_T_L*exp(es_beta[2]*Z_M_L)*Z_M_L )*c(T_baseline[1:nl],T_baseline[nl])*(L_time-L_time_0))/n
  
  a_g_clculas<- colSums(rho_t_L1*Ga_M_L[Indix1,]*Y_t1_L)/colSums(rho_t_L1*exp(es_beta1[2]*Z_M_L[Indix1,])*Y_t1_L)
  a_h_clculas<- colSums(rho_t_L1*exp(es_beta1[2]*Z_M_L[Indix1,])*Z_M_L[Indix1,]*Y_t1_L)/colSums(rho_t_L1*Y_t1_L*exp(es_beta1[2]*Z_M_L[Indix1,]))
  
  delta_alpha_inte2<-     sum(colSums((mu_L_M1- mu_t_L)*exp(es_beta1[2]*Z_M_L))* a_g_clculas*(L_time-L_time_0))/n
  
  d_Lambda_T_L<- colSums(rho_t_L1*(DN_1_L_M-Y_t1_L*es_beta1[1]*Ga_M_L[Indix1,]*Time_clculas_L_M[Indix1,])  )/colSums(Y_t1_L*rho_t_L1*exp(es_beta1[2]*Z_M_L[Indix1,]))                     
  
  delta_eta_inte2<-   sum( colSums((mu_L_M1- mu_t_L)*exp(es_beta1[2]*Z_M_L))* a_h_clculas*d_Lambda_T_L)/n  #####????T_baseline0
  delta_inte1<-  c(delta_alpha_inte1+delta_alpha_inte2,delta_eta_inte1+delta_eta_inte2)
  
  
  E_ee<-  -colSums((mu_L_M1- mu_t_L)*exp(es_beta1[2]*Z_M_L))/n
  a_0_clculas<-  1/colSums(rho_t_L1*Y_t1_L*exp(es_beta1[2]*Z_M_L[Indix1,]))
  a_0_clculas_M<- matrix(1,n1,1)%*%t(as.matrix(a_0_clculas))
  D_lambda_L<- colSums(rho_t_L1*(DN_1_L_M-Y_t1_L*es_beta1[1]*Ga_M_L[Indix1,]*Time_clculas_L_M[Indix1,]))
  D_lambda_L_M<- matrix(1,n1,1)%*%t(as.matrix( D_lambda_L))
  
  ff<- E_ee*a_0_clculas
  ff_M<- matrix(1,n1,1)%*%t(as.matrix(ff))
  DM_t1_L<- DN_1_L_M-Y_t1_L*es_beta1[1]*Ga_M_L[Indix1,]*Time_clculas_L_M[Indix1,]-Y_t1_L*exp(es_beta1[2]*Z_M_L[Indix1,])*a_0_clculas_M* D_lambda_L_M  
  
  delta_inte2<- rowSums(ff_M*DM_t1_L) ##### 
  
  coff<- t(as.matrix(delta_inte1) )%*%solve(Sigma_T1[[2]])
  Inte1<- coff[1,1]*Sigma_T1[[3]]+coff[1,2]*Sigma_T1[[4]]
  del<- (mu_L1-M_time1)/n
  phi1<-rbind(Inte1+ delta_inte2+ del[Indix1],as.matrix(del[Indix0]))
  
  ####.............Estimator of G0..................
  mu_L0<-  S_T_L0%*%as.matrix(L_time-L_time_0) #Ϊһ??ʱ??LΪֹ??Sij??t???????Ļ???
  mu_L_M0<-  as.matrix(mu_L0)%*%matrix(1,1,nl+1) 
  
  rho_t_L0<- cbind(rho_t0[,1:nl],rho_t0[,nl])####???????????????
  
  mu_t_L<-  (S_T_L0*Time_clculas_L_M)%*%AA1  ###????u_ij(t)
  
  Yn_5<- (colSums( (mu_L_M0- mu_t_L)*exp(es_beta0[2]*Z_M_L)*Z_M_L)/n)/colSums(Y_t0_L*rho_t_L0*exp(es_beta0[2]*Z_M_L[Indix0,]))
  ####   Yn_5_M<- matrix(1,n0,1)%*%t(as.matrix(Yn_5))
  
  delta_eta_1<-  sum(Yn_5*colSums(rho_t_L0*DN_0_L_M))
  delta_eta_2<-  sum(Yn_5*colSums(Y_t0_L*rho_t_L0*es_beta0[1]*Ga_M_L[Indix0,]*Time_clculas_L_M[Indix0,])  )
  delta_eta_inte1<- -delta_eta_1+ delta_eta_2 #-sum((S_T_L0*T_baseline1_M[,1:(nl+1)]*exp(es_beta0[2]*Z_M_L)*Z_M_L)%*%(L_time-L_time_0))/n
  
  delta_alpha_inte1<- -sum( colSums(S_T_L0* Ga_M_L)*L_time*(L_time-L_time_0))/n
  ###  delta_eta_inte1<-     -sum(colSums(S_T_L*exp(es_beta[2]*Z_M_L)*Z_M_L )*c(T_baseline[1:nl],T_baseline[nl])*(L_time-L_time_0))/n
  
  a_g_clculas<- colSums(rho_t_L0*Ga_M_L[Indix0,]*Y_t0_L)/colSums(rho_t_L0*exp(es_beta0[2]*Z_M_L[Indix0,])*Y_t0_L)
  a_h_clculas<- colSums(rho_t_L0*exp(es_beta0[2]*Z_M_L[Indix0,])*Z_M_L[Indix0,]*Y_t0_L)/colSums(rho_t_L0*Y_t0_L*exp(es_beta0[2]*Z_M_L[Indix0,]))
  
  delta_alpha_inte2<-     sum(colSums((mu_L_M0- mu_t_L)*exp(es_beta0[2]*Z_M_L))* a_g_clculas*(L_time-L_time_0))/n
  
  d_Lambda_T_L<- colSums(rho_t_L0*(DN_0_L_M-Y_t0_L*es_beta0[1]*Ga_M_L[Indix0,]*Time_clculas_L_M[Indix0,])  )/colSums(Y_t0_L*rho_t_L0*exp(es_beta0[2]*Z_M_L[Indix0,]))                     
  
  delta_eta_inte2<-   sum( colSums((mu_L_M0- mu_t_L)*exp(es_beta0[2]*Z_M_L))* a_h_clculas*d_Lambda_T_L)/n  #####????T_baseline0
  delta_inte1<-  c(delta_alpha_inte1+delta_alpha_inte2,delta_eta_inte1+delta_eta_inte2)
  
  
  E_ee<-  -colSums((mu_L_M0- mu_t_L)*exp(es_beta0[2]*Z_M_L))/n
  a_0_clculas<-  1/colSums(rho_t_L0*Y_t0_L*exp(es_beta0[2]*Z_M_L[Indix0,]))
  a_0_clculas_M<- matrix(1,n0,1)%*%t(as.matrix(a_0_clculas))
  D_lambda_L<- colSums(rho_t_L0*(DN_0_L_M-Y_t0_L*es_beta0[1]*Ga_M_L[Indix0,]*Time_clculas_L_M[Indix0,]))
  D_lambda_L_M<- matrix(1,n0,1)%*%t(as.matrix( D_lambda_L))
  
  ff<- E_ee*a_0_clculas
  ff_M<- matrix(1,n0,1)%*%t(as.matrix(ff))
  DM_t0_L<- DN_0_L_M-Y_t0_L*es_beta0[1]*Ga_M_L[Indix0,]*Time_clculas_L_M[Indix0,]-Y_t0_L*exp(es_beta0[2]*Z_M_L[Indix0,])*a_0_clculas_M* D_lambda_L_M  
  
  delta_inte2<- rowSums(ff_M*DM_t0_L) ##### 
  
  coff<- t(as.matrix(delta_inte1) )%*%solve(Sigma_T0[[2]])
  Inte1<- coff[1,1]*Sigma_T0[[3]]+coff[1,2]*Sigma_T0[[4]]
  del<- (mu_L0-M_time1)/n
  phi0<-rbind(as.matrix(del[Indix1]),Inte1+ delta_inte2+ del[Indix0])
  ###---------------------------delta?????��?---------------------------
  phi<- sum((phi1-phi0)^2)
  list(M_time=M_time0-M_time1,phi=phi)
  #list(theta=es_theta,beta=es_beta,sigma_theta<-sigma_theta,sigma_beta<-sigma_beta,M_time=M_time,phi=phi)
}