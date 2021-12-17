

rm(list = ls())

## first set to current directory using setwd
#setwd("C:/....")


## Calling some libraries
library(deSolve)
library(DEoptim)



#### Viral dynamics model

COV_VL_model <- function(t,x,params){
  with(as.list(x),{   
    
    if(t>tau){
      rbx=rb
    }else{
      rbx=0
    }
    
    
    ddt_S=  - beta*V*S 
    ddt_I= beta*V*S  - delta*I*I^k  - rho*I
    ddt_V= pi*I     - c*V     #   *(1 - ((A3T/V3T)/((A3T/V3T)+EC50)*Tx))
    
    ddt_Sb=  rbx*(Sb)*(1-((Sb+Ib+Rb)/Sb_0))- betab*Vb*Sb   - rhob*Sb*Ib
    ddt_Ib= betab*Vb*Sb  - deltab*Ib*Ib^kb  
    ddt_Vb= pib*Ib     - c*Vb  #  *(1 - ((A3T/V3T)/((A3T/V3T)+EC50b)*Tx))
    ddt_Rb = rhob*Sb*Ib + rbx*(Rb)*(1-((Sb+Ib+Rb)/Sb_0))
    ddt_Gb= deltab*Ib*Ib^kb -  rbx*(Sb+Rb)*(1-((Sb+Ib+Rb)/Sb_0)) 
    
 
    der <- c( ddt_S, ddt_I, ddt_V,  ddt_Sb,  ddt_Ib,  ddt_Vb, ddt_Rb , ddt_Gb)
    
    list(der)
  })       
}



#### Load Viral load data file (first nasal and then BAL)
BAL_data <- read.csv("BAL_data_untreated.csv", header=TRUE,stringsAsFactors=FALSE) # Viral_Loads.csv is the file that records viral loads observed for each patient over the course of their infection
#time  	ID   	Reg_treatment	   Y    	variable    	Treatment_status

NASAL_data <- read.csv("NASAL_data_untreated.csv", header=TRUE,stringsAsFactors=FALSE) # Viral_Loads.csv is the file that records viral loads observed for each patient over the course of their infection
#time  	ID   	Reg_treatment	   Y    	variable    	Treatment_status


## Determine all unique RM IDs in the dataset
IDs = unique(NASAL_data[,"ID"])
print(IDs)

#### Defining final fig structure
nrows = 4
ncolumns = 4
par(mfrow = c(nrows,ncolumns))
par(mfrow = c(nrows,ncolumns),     
    oma = c(0, 0, 0, 0), # 
    mar = c(3.5, 3.5, 0.0, 0.0), 
    mgp = c(2.5, 1, 0),    
    xpd = FALSE)


#### Loading Parameter set that fits each RM
PAR_VL_all=read.csv("Estimated_params_untreated_RDV_PK.csv", header=TRUE,stringsAsFactors=FALSE)
head(PAR_VL_all)

ijk=1

count_entry=1

# starting plotting  viral load data and simulated data for all RMs
for(ID in IDs){   # [c(2,3)]){
  i = which(IDs %in% ID)
  print(i)
  print(ID)
  PAR_VL = PAR_VL_all[i,] # Finding parameters estimated for a RM with identifier ID
 
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
  
  tau=PAR_VL$tau
  
  rho=0
  
#### LOading VL data as observed Order (#time  	ID   	Reg_treatment	   Y    	variable    	Treatment_status)
  
  # Starting with animal with identifier ID  
  pBAL_unt = 10^as.numeric(BAL_data[which(ID==BAL_data[,"ID"]),"Y"])
  timeBAL_unt = as.numeric(BAL_data[which(ID==BAL_data[,"ID"]),"time"])
  
  pNASAL_unt = 10^as.numeric(NASAL_data[which(ID==NASAL_data[,"ID"]),"Y"])
  timeNASAL_unt = as.numeric(NASAL_data[which(ID==NASAL_data[,"ID"]),"time"])

  tzero=0
  tend=22
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
  outCov_VL <- as.data.frame(lsoda(init.x,t.out,COV_VL_model,params))
  
  t1=outCov_VL$time
  NASAL_out=outCov_VL$V
  BAL_out=outCov_VL$Vb
  

  if(ID=="RM17" || ID=="RM18" || ID=="RM19" || ID=="RM20"){
    xlabel="Time (days)"
    
  } else{
    xlabel= ""
    
  }
  
  if(ID=="RM7" || ID=="RM11" || ID=="RM15" || ID=="RM19"){
    ylabel="log10 SARS CoV-2 RNA"
    
  } else{
    ylabel= ""
    
  }
  
plot(t1,NASAL_out,log="y", bty="n", type="l",xaxt="n",yaxt="n",lwd=1,
       xlim=c(tzero,tend),
       ylim=c(10^0,10^9),
       cex.main=1.4,cex.axis=1.4,cex.lab=1.4,
       xlab=xlabel,
       ylab=ylabel,
       col="#E7298A")
  axis(1,at=axTicks(1),labels=as.character(axTicks(1)),cex.axis=1.4,font=1,lwd=1)#,las=1)
  axis(2,at=c(10^0,10^1,10^2,10^3,10^4,10^5,10^6,10^7,10^8,10^9,10^10),labels=expression("0"," ","2"," ","4"," ","6"," ","8"," ","10"),
       cex.axis=1.4,font=1,lwd=1)#,las=1)

  points(timeNASAL_unt,pNASAL_unt,pch=21,bg="#E7298A",cex=1.2)
  
  lines(t1,BAL_out,lwd=1,col="#7570B3")
  
  points(timeBAL_unt,pBAL_unt,pch=22,bg="#7570B3",cex=1.2)
  
  abline(h=LD,lty=2)
  if(ID=="RM20"){
  legend(15,10^10,c("Nasal Swabs","BAL"),bty="n",
         pch=c(21,22),pt.bg=c("#E7298A","#7570B3"),cex=1.2)
  }
  text(18,10^8.8,ID,cex=1.2)
  
  
  ijk=ijk+1
  
  

}
  
