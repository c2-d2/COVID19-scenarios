require(EnvStats)
require(ggplot2)
require(deSolve)
require(caTools)
require(incidence)

# function for estimating cases given input parameters 

#### Inputs
# nsim is number of simulations
# N is the total population (simulations stop if this many people get infected)
# T is length of simulation in days

# num_introductions is number of expected introductions with timing calibrated to an exponentially increasing model until travel ban and then flat
# inc_shape and scale are parameters for Weibull distribution of incubation period (average ~5.8 for baseline parameters)

# pre_symp is proportion of transmission occurring pre_symptomatically
# inf_max is max days post symptoms of infectiousness used for triangular distributin of onward infection time relative to symptom onset; 
# R0_reduce is relative R0 for those without symptoms
# R0_symp is basic reproductive number for those with symptoms

# isolate_shape and isolate_rate are parameters for gamma distribution for time from symptom onset to self observation
# isolate_prob is probability an infected person will self-observe and limit contacts once symptoms develop; assume transmission stops once isolated
# symp_prob is proportion of people with symptoms

estimate_cases <- function(nsim,N,T,num_introductions,inc_shape,inc_scale,
                           pre_symp,inf_max,R0_reduce,R0_symp,isolate_shape,isolate_rate,isolate_prob,symp_prob){
  
  # create data frame to store results
  epidemic_curve_master <- as.data.frame(cbind("Sim"=double(),"Day_infected"=double(),"Cases"=double(),"Imports"=double(),"Severe_cases"=double(),"Cum_cases"=double(),"Cum_severe_cases"=double()))

  # keep track of serial interval (time from symptom onset to symptom onset for infector-infectee), doubling times and growth rates
  SI <- rep(NA,nsim)
  Doubling_times <- rep(NA,nsim)
  Growth_rates <- rep(NA,nsim)
  total_cases <- rep(NA,nsim)
  
  # function for infectiousness
  # triangular distribution - choose min based on maximum day of infectiousness possible and % pre_symptomatic
  # make mode slightly off center so not completely symmetric distribution
  rtrizero_t_round <- function(n, pre_symp, max){
    if(n<=0){return(NA)}
    min <- round((max-max/(1-pre_symp)))
    return(round(rtri(n,min,max,mode=min+round((max-min)/2+1))))
  }
  
  # function for making sure people can't be infected before their infector; instead make them infected the next day 
  timeline <- function(n,t){
    ifelse(n<=t,t+1,n)
  }
  
  # function to set infection time for asymptomatic (uniform distribution)
  asymp <- function(n,t,inf_max){
    if(is.na(n[1])){
      return(t+round(runif(length(n),1,inf_max)))
    } else{
      return(n)
    }
  }

  
  # epidemic curve used for external force of infection
  model <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      beta <- betahat * (1 - a2/(1 + exp(-a1 * (t - atau))))
      dS <- -beta * S * I / (S+E+I+R)
      dE <- beta * S * I/  (S+E+I+R) - sigma *E
      dI <- sigma * E - gamma * I
      dR <- gamma * I
      list(c(dS,dE,dI,dR))
    })
  }
  
  # number of people in source population
  NS <- 11000000
  # starting values
  y<- c(S=NS-1,E=0,I=1,R=0)
  times<-seq(0,53,1)
  parms<-c(betahat=0.95,a1=0.1,a2=0.6,atau=10,sigma=0.2,gamma=0.2)
  out<-as.data.frame(lsoda(y,times,model,parms))

  # make flat until time T
  out[54:T+1,] <- out[54,]
  out[54:T+1,1] <- seq(54,T,1)

  # Calibrate extF to the number of introductions, given the progression of the epidemic in the source population
  sumsqrt <- sum(sqrt(N))
  extF <- -log(1-num_introductions/(sqrt(N)*sumsqrt))/trapz(seq(0,T,1),out$I)
  
  for (sim in 1:nsim){
    cat(sim,"\n")
    epidemic_curve <- as.data.frame(cbind("Day_infected"=seq(1,T,1)))
    
    # create a dataframe to track cases
    cases <- as.data.frame(cbind("ID"=rep(NA,N),"Day_infected"=rep(NA,N),"Import"=rep(0,N),"Symp_prob"=rep(NA,N),"Symptom_onset_day"=rep(NA,N),
                                 "Isolation"=rep(NA,N),"Isolation_time"=rep(NA,N),"Infector"=rep(NA,N),"Severity"=rep(NA,N)))
    

    # determine # and timing of outside introductions based on external force of infection
    day_introduced <- c()
    for (i in 1:T){
      # Probability of infection
      prob_inf_fromsource <- 1 - exp(-extF*out$I[i])
      num_inf <- rbinom(1,N,prob_inf_fromsource)
      if (!is.na(num_inf)){
        day_introduced <- c(day_introduced,rep(i,num_inf))
      }
    }
    
    # determine if symptomatic
    symp <- rbinom(length(day_introduced),1,symp_prob)
    
    # incubation period if symptomatic
    inc_period <- ifelse(symp==1,day_introduced + round(rweibull(length(day_introduced),inc_shape,inc_scale)),NA)
    
    # determine if isolated and when
    isolation <- ifelse(symp==1,rbinom(length(day_introduced),1,isolate_prob),NA)
    isolation_time <- inc_period + round(ifelse(isolation==1,rgamma(1,isolate_shape,isolate_rate),NA))
    
    # determine if severe
    severe <- ifelse(symp==1,rbinom(length(day_introduced),1,0.18),0)

    if (length(day_introduced)>0){
      cases[1:length(day_introduced),]<- c(seq(1:length(day_introduced)),day_introduced,rep(1,length(day_introduced)),symp,inc_period,isolation,isolation_time,rep(0,length(day_introduced)),severe)
    }
    
    num_infectees <- length(day_introduced)
    
    # find who gets infected on day t and then determine how many people they will infect and when
    for (t in 1:T){
      
      # find who gets infected on day t
      new_infectors <- cases$ID[which(cases$Day_infected==t)]
      num_new_infectors <- length(new_infectors)
      
      # keep track of if/when the infectors develop symptoms and are isolated
      isolate_times <- as.list(cases$Isolation_time[which(cases$Day_infected==t)])
      symptom_times <- as.list(cases$Symptom_onset_day[which(cases$Day_infected==t)])

      if (num_new_infectors >0){
        
          # use R0 distribution to get # of potential secondary infectees for each new infector
          R0_asymp <- R0_symp*R0_reduce
          # choose number to infect based on R0 given symptoms using negative binomial distribution
          num_pot_infectees <- rnbinom(num_new_infectors, size=0.5,  mu=c(ifelse(is.na(unlist(symptom_times)),R0_asymp,R0_symp)))
          #num_pot_infectees <- rpois(num_new_infectors,c(ifelse(is.na(unlist(symptom_times)),R0_asymp,R0_symp))) # Also could make R0 negative binomial to allow for superspreading
          
          if (sum(num_pot_infectees)>0){
            # identify times will be infecting relative to symptom onset
            time_inf<- lapply(num_pot_infectees, rtrizero_t_round, pre_symp, inf_max) 
            
            # add to symptom times
            time_pot_new_infectees_1 <-lapply(seq_along(time_inf),function(i)unlist(time_inf[i])+unlist(symptom_times[i]))
            
            # those without symptom times will have NA for infection time, so instead draw from a uniform distribution between day 1 and maximum day of infectiousness
            time_pot_new_infectees_1 <- lapply(time_pot_new_infectees_1,asymp,t,inf_max)
                                                                    
            # if infection time is before or on same day as that of infectors, make it the next day
            time_pot_new_infectees <- lapply(time_pot_new_infectees_1, timeline, t) 
              
            # if self-observation / limiting contacts is in place, see how many people will be infected before isolation 
            x <- lapply(seq_along(time_pot_new_infectees),function(i)which(unlist(time_pot_new_infectees[i])<unlist(isolate_times[i]) | is.na(unlist(isolate_times[i]))))
            time_new_infectees <- lapply(seq_along(time_pot_new_infectees),function(i)unlist(time_pot_new_infectees[i])[unlist(x[i])])
            
            # get number of new infectees
            num_new_infectees <- length(unlist(time_new_infectees))
            
            # update cases table if there are new cases and the entire population hasn't been infected
            if (num_new_infectees >0 & ((num_infectees + num_new_infectees)<N)){
              # assign ID
              cases$ID[(num_infectees+1):(num_infectees+num_new_infectees)] <- seq((num_infectees+1),(num_infectees+num_new_infectees),1)
              
              # record day infected
              cases$Day_infected[(num_infectees+1):(num_infectees+num_new_infectees)] <- unlist(time_new_infectees)
              
              # record if will be symptomatic
              cases$Symp_prob[(num_infectees+1):(num_infectees+num_new_infectees)] <-  rbinom(num_new_infectees,1,symp_prob)
              
              # use incubation period to determine day of symptom onset
              cases$Symptom_onset_day[(num_infectees+1):(num_infectees+num_new_infectees)] <- 
                cases$Day_infected[(num_infectees+1):(num_infectees+num_new_infectees)]  + 
                ifelse(cases$Symp_prob[(num_infectees+1):(num_infectees+num_new_infectees)]==1,(round(rweibull(num_new_infectees,inc_shape,inc_scale))),NA)
              
              # if symptoms determine if they will self isolate
              cases$Isolation[(num_infectees+1):(num_infectees+num_new_infectees)] <-  
                ifelse(cases$Symp_prob[(num_infectees+1):(num_infectees+num_new_infectees)]==1,rbinom(num_new_infectees,1,isolate_prob),NA)
              
              # determine time after symptom onset isolation/isolation will occur
              cases$Isolation_time[(num_infectees+1):(num_infectees+num_new_infectees)] <-
                cases$Symptom_onset_day[(num_infectees+1):(num_infectees+num_new_infectees)] +
                round(ifelse(cases$Isolation[(num_infectees+1):(num_infectees+num_new_infectees)]==1,rgamma(num_new_infectees,isolate_shape,isolate_rate),NA))
              
              # record their infector (for checking serial interval)
              cases$Infector[(num_infectees+1):(num_infectees+num_new_infectees)] <- rep(new_infectors,c(unlist(lapply(time_new_infectees,length))))

              # record if severe (18% of symptomatic)
              cases$Severity[(num_infectees+1):(num_infectees+num_new_infectees)] <- ifelse(cases$Symp_prob[(num_infectees+1):(num_infectees+num_new_infectees)]==1,rbinom(num_new_infectees,1,0.18),0)
              
              # add to total of number of people infected
              num_infectees <- num_infectees+num_new_infectees
             } 
          }
        }
    }
    
    # record serial interval
    cases$infector_symptom_day[cases$Import!=1] <- cases$Symptom_onset_day[cases$Infector]
    SI[sim] <- mean(cases$Symptom_onset_day-cases$infector_symptom_day,na.rm=TRUE)
    
  
    # update epidemic curve table
    summary <-  as.data.frame(aggregate(cases$ID,by=list(cases$Day_infected),length))
    
    # if no cases, make a data frame with all 0s
    if (nrow(summary)==0){
      summary<- as.data.frame(cbind(seq(1,T,1),rep(0,T)))
    } 
    names(summary) <- c("Day_infected","Cases")
    
    # estimate doubling times and growth rate
    summary_grow <-  as.data.frame(aggregate(cases$ID,by=c(list(cases$Day_infected),list(cases$Import)),length))
    names(summary_grow) <- c("Day_infected","Import_status","Cases")
    # only calculate growth rate once there are at least 5 nonimported cases
    cutoff <- min(summary_grow$Day_infected[which(cumsum(summary_grow$Cases[summary_grow$Import_status==0])>=5)])
    summary_grow2 <- as.data.frame(aggregate(cases$ID,by=list(cases$Day_infected),length))
    names(summary_grow2) <- c("Day_infected","Cases")
    summary_grow2 <- summary_grow2[summary_grow2$Day_infected >= cutoff & summary_grow2$Day_infected<=T,]
    
    # If N cases reached, only calculate growth rate while # cases per day increasing 
      # (but use N/10 as cutoff because not all N cases recorded because cuts off when cases will be >N and some occur after T)
    total_case <- 0
    if (nrow(summary_grow2)>1){
      # total case = 1 if cases by end are greater than 3 times the number of introductions, 0 otherwise
      total_case <- ifelse(sum(summary_grow2$Cases)>(3*num_introductions),1,0)
      if (sum(summary_grow2$Cases)>(N/10)){
        summary_grow2 <- summary_grow2[summary_grow2$Day_infected<=summary_grow2$Day_infected[which.max(summary_grow2$Cases)],]
      }
    }
    total_cases[sim] <- total_case
    # change to dates for incidence function
    growth <- rep(as.Date("11/30/19","%m/%d/%y") + summary_grow2$Day_infected,summary_grow2$Cases)
    
    # estimate doubling time and growth rate
    if (length(unique(growth))>1){
      inc <- incidence(growth)
      analysis <- fit(inc)
      Doubling_times[sim] <- ifelse(length(analysis$info$doubling)>0,analysis$info$doubling,NA)
      Growth_rates[sim] <-  ifelse(length(analysis$info$r)>0,analysis$info$r,NA)
    } 
    epidemic_curve <- merge(epidemic_curve,summary,all.x = TRUE)
    
    # keep track of importations
    import_summary <- as.data.frame(aggregate(cases$Import,by=list(cases$Day_infected),sum))
    if (nrow(import_summary)==0){
      import_summary<- as.data.frame(cbind(seq(1,T,1),rep(0,T)))
    } 
    names(import_summary) <- c("Day_infected","Imports")
    epidemic_curve <- merge(epidemic_curve,import_summary,all.x = TRUE)
    
    # keep track of severity
    severity_summary <- as.data.frame(aggregate(cases$Severity,by=list(cases$Day_infected),sum,na.rm=TRUE))
    if (nrow(severity_summary)==0){
      severity_summary<- as.data.frame(cbind(seq(1,T,1),rep(0,T)))
    } 
    names(severity_summary) <- c("Day_infected","Severe_cases")
    epidemic_curve <- merge(epidemic_curve,severity_summary,all.x = TRUE)
    
    epidemic_curve <- cbind(rep(sim,T),epidemic_curve)
    names(epidemic_curve)[1] <- "Sim"
    
    epidemic_curve$Cases[is.na(epidemic_curve$Cases)] <- 0
    epidemic_curve$Severe_cases[is.na(epidemic_curve$Severe_cases)] <- 0
    
    # calculate cumulative sums
    epidemic_curve$Cum_cases <- cumsum(epidemic_curve$Cases)
    epidemic_curve$Cum_severe_cases <- cumsum(epidemic_curve$Severe_cases)

    # remove new cases after peaks for ones where N people are infected
    if (sum(epidemic_curve$Cases)>(N/10)){
      epidemic_curve$Cases[epidemic_curve$Day_infected>epidemic_curve$Day_infected[which.max(epidemic_curve$Cases)]] <- NA
    }
    
    epidemic_curve_master <- rbind(epidemic_curve_master,epidemic_curve)
  } 
  
  # summarize
  epidemic_curve_summary <- aggregate(epidemic_curve_master$Cases,by=list(epidemic_curve_master$Day_infected),mean,na.rm=TRUE)
  epidemic_curve_summary2 <- aggregate(epidemic_curve_master$Cum_cases,by=list(epidemic_curve_master$Day_infected),mean,na.rm=TRUE)
  epidemic_curve_summary3 <- aggregate(epidemic_curve_master$Severe_cases,by=list(epidemic_curve_master$Day_infected),mean,na.rm=TRUE)
  epidemic_curve_summary4 <- aggregate(epidemic_curve_master$Cum_severe_cases,by=list(epidemic_curve_master$Day_infected),mean,na.rm=TRUE)
  epidemic_curve_summary <- cbind(epidemic_curve_summary,epidemic_curve_summary2[,2],epidemic_curve_summary3[,2],epidemic_curve_summary4[,2])
  names(epidemic_curve_summary) <- c("Day_infected","Average_cases","Average_cumulative_cases","Average_severe_cases","Average_cumulative_severe_cases")
  epidemic_curve_master2<- merge(epidemic_curve_master,epidemic_curve_summary)
  
  epidemic_curve_master2 <- cbind(epidemic_curve_master2,
                                   "isolate_prob"=rep(isolate_prob,nrow(epidemic_curve_master2)),
                                    "symp_prob"=rep(symp_prob,nrow(epidemic_curve_master2)),
                                    "num_introductions"=rep(num_introductions,nrow(epidemic_curve_master2)),
                                    "pre_symp"=rep(pre_symp,nrow(epidemic_curve_master2)),
                                    "R0_reduce"=rep(R0_reduce,nrow(epidemic_curve_master2)),
                                    "R0"=rep(R0_symp,nrow(epidemic_curve_master2)))
  
  # make start day December 1, 2019
  date <- as.Date("11/30/19","%m/%d/%y")
  # convert days to dates from today
  epidemic_curve_master2$date <- epidemic_curve_master2$Day_infected + date
  
  # save file
  write.table(epidemic_curve_master2,paste0("epidemic_curve_master_",R0_symp,"_",R0_reduce,"_",isolate_prob,"_",symp_prob,"_",isolate_shape,"_",num_introductions,"_",pre_symp,".csv"),sep=",",col.names = FALSE)
  
  summary_stats <- cbind(SI,Doubling_times,Growth_rates,total_cases,
                   "isolate_prob"=rep(isolate_prob,nsim),
                   "symp_prob"=rep(symp_prob,nsim),
                   "num_introductions"=rep(num_introductions,nsim),
                   "pre_symp"=rep(pre_symp,nsim),
                   "R0_reduce"=rep(R0_reduce,nsim),
                   "R0"=rep(R0_symp,nsim))
 
  write.table(summary_stats,paste0("SI_",R0_symp,"_",R0_reduce,"_",isolate_prob,"_",symp_prob,"_",isolate_shape,"_",num_introductions,"_",pre_symp,".csv"),sep=",",col.names = FALSE)
  
  return(epidemic_curve_master2)

}

# the parameters without preset values are the ones that can be varied in the shiny app
simulation <- estimate_cases(nsim=50,N=100000,T=152,num_introductions,inc_shape=2.39,inc_scale=6.54,
                             pre_symp,inf_max=6,R0_reduce,R0_symp,isolate_shape=1.5,isolate_rate=0.9,isolate_prob,symp_prob)


