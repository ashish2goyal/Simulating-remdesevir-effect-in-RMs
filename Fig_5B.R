rm(list = ls())

## first set to current directory using setwd
#setwd("C:/....")

## Calling some libraries
library(deSolve)
library(DEoptim)

#### Viral dynamics model
COV_VL_model_exchange <- function(t,x,params){
  with(as.list(x),{   
    
    if(t>0.5+nu_BAL){
      EC50bx=EC50b
    } else{
      EC50bx=10^20
    }
    
    if(t>0.5+nu){
      EC50x=EC50
    } else{
      EC50x=10^20
    }
    
    
    if(t>tau){
      rbx=rb
    }else{
      rbx=0
    }

    ddt_S=  - beta*V*S 
    ddt_I= beta*V*S  - delta*I*I^k   - rho*I
    ddt_V= pi*I *(1 - ((A3T(t)/VD)/((A3T(t)/VD)+EC50x)))     - c*V     #   *(1 - ((A3T(t)/VD)/((A3T(t)/VD)+EC50)))
    
    ddt_Sb=  rbx*(Sb)*(1-((Sb+Ib+Rb)/Sb_0))- betab*Vb*Sb   - rhob*Sb*Ib
    ddt_Ib= betab*Vb*Sb  - deltab*Ib*Ib^kb  
    ddt_Vb= pib*Ib*(1 - ((A3T(t)/VD)/((A3T(t)/VD)+EC50bx)))     - c*Vb  #  *(1 - ((A3T(t)/VD)/((A3T(t)/VD)+EC50b)))
    ddt_Rb = rhob*Sb*Ib + rbx*(Rb)*(1-((Sb+Ib+Rb)/Sb_0))
    ddt_Gb= deltab*Ib*Ib^kb -  rbx*(Sb+Rb)*(1-((Sb+Ib+Rb)/Sb_0))
 
    der <- c( ddt_S, ddt_I, ddt_V,  ddt_Sb,  ddt_Ib,  ddt_Vb, ddt_Rb , ddt_Gb)
    
    list(der)
  })       
}

RDV_PK_model <- function(t,x,params){
  with(as.list(x),{   
    
    
    ddt_A1 = - (24*0.69/1)*A1   - k12*A1     - k1T*A1     + k1e*A1T                 #	 prodrug GS-5734      (half-life _ 1 hr, https://sidp.org/resources/Documents/COVID19/Matt%20Remdesivir%20Handouts%204.7.2020.pdf)   
    
    ddt_A2 =   -   (24*0.69/24)*A2   +     k12*A1   - k23*A2     -   k2T *  A2   + k2e*A2T  #  	k2=(24*0.69/24)	  alanine metabolite (GS-704277)   (half-life -- 24 hr, https://sidp.org/resources/Documents/COVID19/Matt%20Remdesivir%20Handouts%204.7.2020.pdf)
    
    ddt_A3 =   - k3* A3   +  k23 *  A2   -   k3T *  A3    + k3e*A3T   #	    nucleoside (GS-441524)  
    
    
    ##   ;; https://sidp.org/resources/Documents/COVID19/Matt%20Remdesivir%20Handouts%204.7.2020.pdf
    
    ##  ; GS-5734 inhibited MHV more potently than GS-44152 with more efficient metabolism of the prodrug into the active nucleoside triphosphate by bypassing the rate-limiting first phosphorylation step (https://mbio.asm.org/content/mbio/9/2/e00221-18.full.pdf)
    
    ddt_A1T =  k1T*A1          -   k12*A1T          -   k1e*A1T               # 	Tissue	 prodrug GS-5734        
    
    ddt_A2T =  k2T *  A2      +   k12*A1T           -  k23*A2T        -   k2e*A2T      -  k24T*A2T  #	Tissue   alanine metabolite (GS-704277)
    
    ddt_A3T =   k3T *  A3     + k23*A2T            -   k34T *  A3T      - k3e*A3T    # Tissue  nucleoside (GS-441524) 
    
    ddt_A4T =  k24T*A2T        + k34T* A3T             -   (24*0.69/24)*A4T     # Tissue GS-441524 requires intracellular phosphorylation via cellular kinases to pharmacologically active metabolite of remdesivir which is the triphosphate of GS-441524    (nucleoside triphosohate )
    #; the active triphosphate metabolite has a prolonged intracellular half-life of more than 35 hours (https://www.who.int/ebola/drc-2018/summaries-of-evidence-experimental-therapeutics.pdf?ua=1)	
    #; also, In PBMCs, NTP represents the predominant metabolite and was persistent with a t1/2 of 14?h (NATURE EBOLA PAPER)		
    #; also NTP half life of approximately 20 hours in humans (https://www.ncbi.nlm.nih.gov/pubmed/28659436)											
    
    
    der <- c( ddt_A1, ddt_A2, ddt_A3,  ddt_A1T,  ddt_A2T,  ddt_A3T,  ddt_A4T)
    
    list(der)
  })       
}



#### Load Viral load data file (first nasal and then BAL)
BAL_data <- read.csv("BAL_data_treated.csv", header=TRUE,stringsAsFactors=FALSE) # Viral_Loads.csv is the file that records viral loads observed for each patient over the course of their infection
#time  	ID   	Reg_treatment	   Y    	variable    	Treatment_status

NASAL_data <- read.csv("NASAL_data_treated.csv", header=TRUE,stringsAsFactors=FALSE) # Viral_Loads.csv is the file that records viral loads observed for each patient over the course of their infection
#time  	ID   	Reg_treatment	   Y    	variable    	Treatment_status


## Determine all unique RM IDs in the dataset
IDs = unique(BAL_data[,"ID"])
print(IDs)

#### Defining final fig structure
nrows = 2
ncolumns = 3
par(mfrow = c(nrows,ncolumns))
par(mfrow = c(nrows,ncolumns),     
    oma = c(0, 0, 0, 0), # 
    mar = c(3.5, 3.5, 0.0, 0.0), 
    mgp = c(2.5, 1, 0),    
    xpd = FALSE)


#### Loading Parameter set that fits each RM
PAR_VL_all=read.csv("Estimated_params_treated_RDV_PK.csv", header=TRUE,stringsAsFactors=FALSE)

ijk=1

# starting plotting  viral load data and simulated data for all RMs
for(ID in IDs){   # [c(2,3)]){
  i = which(IDs %in% ID)
  PAR_VL = PAR_VL_all[i,] # Finding parameters estimated for a RM with identifier ID
  #print(PAR_VL)
  
  tzero=0
  tend=8 
  #print("tzero")
  
  #Loading all RDV PK parameters for plotting simulated data
  k2=PAR_VL$k2_mode
  k1T=PAR_VL$k1T_mode
  k2T=PAR_VL$k2T_mode
  k3T=PAR_VL$k3T_mode
  k1e=PAR_VL$k1e_mode
  k2e=PAR_VL$k2e_mode
  k3e=PAR_VL$k3e_mode
  k12=PAR_VL$k12_mode
  k23=PAR_VL$k23_mode
  k24T=PAR_VL$k24T_mode
  k34T=PAR_VL$k34T_mode
  #print(k34T)
  k3=PAR_VL$k3_mode
  V1=PAR_VL$V1_mode
  V2=PAR_VL$V2_mode
  V3=PAR_VL$V3_mode
  V3T=PAR_VL$V3T_mode
  V4T=PAR_VL$V4T_mode
  
  VA=V1 # Volume of distribution of GS5734 in blood
  VB=V2 # Volume of distribution of ALAMET in blood
  VC=V3 # Volume of distribution of NUC (GS441524) in blood
  VD=V3T  # Volume of distribution of NUC in Tissue 
  
  
  ### Dosing event
  Dose1 = 10  
  Dose2 = 5  
  
  dtreat_ini=0.5
  dtreat_dur=6
  
  dtreat= c(dtreat_ini,seq(from=dtreat_ini+0.5,to=dtreat_dur,by=1))
  
  eventdat <- data.frame(var = rep("A1",length(dtreat)),
                         time = c(dtreat) ,
                         value = c(Dose1,rep(Dose2,length(dtreat)-1)),
                         method = rep("add",length(dtreat)))
  
  #####  Initial conditions for PK RDV model
  A1_0 = 0
  A2_0 = 0
  A3_0 = 0
  A1T_0 = 0
  A2T_0 = 0
  A3T_0 = 0
  A4T_0 = 0
  
  # Solving ODE
  init.x <- c(A1=A1_0,A2=A2_0,A3=A3_0,A1T=A1T_0,A2T=A2T_0,A3T=A3T_0,A4T=A4T_0)
  t.out <- seq(tzero,tend,by=0.001)  
  outRDV_PK <- as.data.frame(lsoda(init.x,t.out,RDV_PK_model,params,
                                   events = list(data=eventdat)))
  
  #print(outRDV_PK)
  
  A3T <- approxfun(outRDV_PK[,c(1,7)], rule = 2)
  
 
  #Loading all RDV PK parameters for plotting simulated data
  beta=10^PAR_VL$beta
  betab=beta
  
  delta=PAR_VL$delta
  deltab=PAR_VL$deltab
  
  k=PAR_VL$k
  pi=10^PAR_VL$pi
  
  pib=10^PAR_VL$pib
  c=15
  
  kb=k
  
  rhob=PAR_VL$rhob
  rb=PAR_VL$rb
  
  rho=0
  
  tau=PAR_VL$tau
  
  EC50=10^PAR_VL$EC50
  EC50b=10^PAR_VL$EC50b
  
  nu=PAR_VL$tau_RDV
  nu_BAL=PAR_VL$tau_RDV_BAL
  

#### LOading VL data as observed Order (#time  	ID   	Reg_treatment	   Y    	variable    	Treatment_status)
  
  # Starting with animal with identifier ID  
  pBAL_unt = 10^as.numeric(BAL_data[which(ID==BAL_data[,"ID"]),"Y"])
  timeBAL_unt = as.numeric(BAL_data[which(ID==BAL_data[,"ID"]),"time"])
  
  pNASAL_unt = 10^as.numeric(NASAL_data[which(ID==NASAL_data[,"ID"]),"Y"])
  timeNASAL_unt = as.numeric(NASAL_data[which(ID==NASAL_data[,"ID"]),"time"])

  LD=10^1 # limit of detection (assumed)
  
  #####  Initial conditions for PK RDV model
  moi=10^-5
  S_0=3*3.7*10^8 / 430  ## ;;;  the average fixed lung volume was 4,300 ml  (https://pubmed.ncbi.nlm.nih.gov/7103258/) and 30mL in nasal cavity (https://pubmed.ncbi.nlm.nih.gov/12563939/) . taking proportional
  ### ; in a 5 Jlm-thick tissue section from just one level of  the nasal cavity, the number ofrespiratory epithelial  cells exceeds 6,000 (https://journals.sagepub.com/doi/pdf/10.1177/019262339001800104)
  Sb_0=3.7*10^8   ## ;;;the mean number of alveolar type II cells was 37 ± 5 billion     https://pubmed.ncbi.nlm.nih.gov/7103258/
  
  I_0= moi* S_0  ### ; based on PFU and (MOI) m<<1, giving rise to prob of inf cells as ~m . With m=0.001, we have I0=0.001*S0 
  Ib_0=  4*moi* S_0  ### ;; TCID intratracheal is 4 times intranasal
  
  V_0=  pi*I_0/c
  Vb_0=pib*Ib_0/c
  
  Rb_0 = 0
  Gb_0 = 10^-2
  
  # Solving ODE
  init.x <- c(S=S_0,I=I_0,V=V_0,Sb=Sb_0,Ib=Ib_0,Vb=Vb_0,Rb=Rb_0,Gb=Gb_0)
  t.out <- seq(tzero,tend,by=0.05)  
  outCov_VL <- as.data.frame(lsoda(init.x,t.out,COV_VL_model_exchange,params))
  
  t1=outCov_VL$time
  NASAL_out=outCov_VL$V
  BAL_out=outCov_VL$Vb
  
  if(ID=="RM4" || ID=="RM5" || ID=="RM6"){
    xlabel="Time (days)"
    
  } else{
    xlabel= ""
      
  }
  
  if(ID=="RM1" || ID=="RM4"){
    ylabel="log10 SARS CoV-2 RNA"
    
  } else{
    ylabel= ""
    
  }
  
  #dev.off()   cbPalette=c("#E7298A", "#7570B3")
  plot(t1,NASAL_out,log="y", bty="n", type="l",xaxt="n",yaxt="n",lwd=1,
       xlim=c(tzero,tend),
       ylim=10^c(0,9),
       cex.main=1.4,cex.axis=1.4,cex.lab=1.4,
       xlab=xlabel,
       ylab=ylabel,
       col="#E7298A")
  axis(1,at=axTicks(1),labels=as.character(axTicks(1)),cex.axis=1.4,font=1,lwd=1)#,las=1)
  axis(2,at=c(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10),labels=expression("0","1 ","2"," 3","4"," 5","6","7","8","9 ","10"),
       cex.axis=1.4,font=1,lwd=1)#,las=1)
  # axis(2,at=c(10^0,10^2,10^4,10^6,10^8,10^10),labels=expression("0","2","4","6","8","10"),
  #     cex.axis=1.4,font=1,lwd=1)#,las=1)
  
  points(timeNASAL_unt,pNASAL_unt,pch=21,bg="#E7298A",cex=1.2)
  
  lines(t1,BAL_out,lwd=1,col="#7570B3")
  
  points(timeBAL_unt,pBAL_unt,pch=22,bg="#7570B3",cex=1.2)
  
  abline(h=LD,lty=2)
  
  legend(4.5,10^9,c("Nasal Swabs","BAL"),bty="n",
         pch=c(21,22),pt.bg=c("#E7298A","#7570B3"),cex=1.2)
  text(1,10^8.8,ID,cex=1.2)
  
  
}
  



