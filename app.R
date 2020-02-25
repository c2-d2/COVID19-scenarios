require(shiny)
require(dplyr)
require(ggplot2)
require(plotly)
require(htmlwidgets)

nCoV_table <- read.csv('./nCoV_table.csv')


ui <- fluidPage(position="left",
  
  titlePanel('Scenario exploration for SARS-CoV-2 infections'),
  h5(tags$a(
    href="https://zenodo.org/badge/latestdoi/242795015",
    tags$img(src="https://zenodo.org/badge/242795015.svg",
             title="DOI"))),
  h5("Developed by the Center for Communicable Disease Dynamics at the Harvard T.H. Chan School of Public Health (Rebecca Kahn, Stephen Kissler, Nishant Kishore and Marc Lipsitch), in collaboration with Health Department colleagues. This is for scenario exploration only. Some scenarios generated may be unrealistic, and users should exercise judgment in interpreting the outputs."), 
  h5("Important note: This app was developed for the specific purpose of considering how many cases might currently be circulating in a setting where testing has not been widespread. It is not designed as a rigorous way to evaluate impacts of interventions."),
  strong(h5("Overview:")),
  tags$ul(tags$li(h5("We calibrated timing of introductions of infected persons into a city based on an epidemic curve with exponential growth beginning on December 1 (based on genomics data, SARS-CoV-2 was estimated to have been first introduced into the human population sometime in November or early December 2019)
      and then leveling off on January 23 (when travel restrictions were implemented in Wuhan) to account for decreased travel from Wuhan but to allow for continued introductions from elsewhere. 
      Light blue dots along the x-axis show dates of introductions. We then modeled the subsequent transmission chains stemming from those introductions.")),
      tags$li(h5("We varied parameters for number of expected introductions, proportion of infections that are symptomatic, 
      probability a person with a symptomatic infection will self-observe and limit contacts, R0, R0 for asymptomatic infection relative to symptomatic infection,
      and proportion of transmission that occurs presymptomatically. 
      Change the inputs on the bars below to see how these affect the outbreak trajectories.")),
      tags$li(h5("The y axis shows the number of either new or cumulative infections in a city, which are a combination of imported infections and their subsequent infections. Results are shown for 5 months."))),
  h5("The black line is the average of the 50 stochastic simulations shown in grey, and the red line is the average number of severe COVID-2019 illnesses. All 50 simulations shown on a given graph have the same parameters. 
     The variation between outbreak trajectories on the same graph comes from the stochasticity in the following:"),
  tags$ul(tags$li(h5("Timing of introductions (the number of infected people introduced on a given day is drawn from a binomial distribution, with the probability of introduction calibrated to the exponential epidemic curve).")), 
          tags$li(h5("Number of secondary infections (drawn from a negative binomial distribution with mean of R0 and dispersion parameter 0.5 (Riou 2020) for each infection). 
                     This distribution allows for superspreaders, which can impact early epidemic dynamics.")), 
          tags$li(h5("Timing of secondary infections (drawn from a triangular distribution, calibrated to the percent of transmission that occurs presymptomatically).")), 
          tags$li(h5("Whether an infected person will have symptoms (drawn from a binomial distribution for each infected person with probability based on the input value selected). Among infections with symptoms, whether or not the infection is severe is drawn from a binomial distribution with probability 0.18 (WHO).")), 
          tags$li(h5("Incubation period for symptomatic infections (drawn from a Weibull distribution with mean 5.8 days for each illness (Backer 2020)).")), 
          tags$li(h5("Whether a person with symptomatic infection will self-observe and limit contacts (drawn from a binomial distribution for each illness with probability based on the input value selected).")), 
          tags$li(h5("The time at which they will start to limit contacts and continue to do so for the duration of their infection (drawn from a gamma distribution with mean 1.5 days after symptoms)."))),
  h5("Important caveats:"),
  tags$ul(tags$li(h5("Onward infections stop if 100,000 people are infected. This branching model only includes infected individuals so does not incorporate depletion of susceptibles or contact patterns. With the exception of potential self-observation and limiting of contacts by symptomatic individuals, it also does not incorporate any control measures.")),
          tags$li(h5("In principle, R0 depends on population density, but there is little evidence that it actually does. Use the slider to change the value of R0.")),
          tags$li(h5("We are working on incorporating flight data to better calibrate the timing of introductions.")),
          tags$li(h5("We have not incorporated any changes in R0 due to seasonality given the limited data and uncertainty. If warmer weather will decrease R0, the total number of infections by the end of April will be overestimated.")),
          tags$li(h5("The table at the top shows the median serial interval (time between symptom onset for infector-infectee pair) and doubling time. 
        Note, doubling times are only calculated once at least 5 local infections have occurred and only for simulations with a total number of infections by the end of April at least 3 times the number of introductions.
        Doubling times between 5-7 days are most consistent with the Wuhan outbreak (Wu 2020), so doubling times shorter than this likely indicate unrealistic parameter combinations. Longer doubling times may reflect impact of self-observation.")),
          tags$li(h5("The static table at the bottom (does not change with inputs) shows the average cumulative infections by the end of February, March and April across values for number of introductions and proportion of presymptomatic transmission for simulations with 50% infections symptomatic, 10% probability of self-observation and limiting contacts, overall R0 of 2, and 0.2 relative R0 for asymptomatic infections."))),
  h5("Findings:"),
  tags$ul(tags$li(h5("The amount of presymptomatic transmission and the extent to which infected people limit contacts (a combination of proportion symptomatic and probability they will self-observe) have a large impact on the total number of infections."))),

  div(style="display:inline-block",shinyWidgets::sliderTextInput(inputId = "num_introductions",
              label = "Number of introductions",
              choices = c(1,2,5,10,25,50,100),selected=1)),
  div(style="display:inline-block",shinyWidgets::sliderTextInput(inputId = "symp_prob",
              label = "Proportion of infections that are symptomatic",
              choices = c(0.5,0.7,0.9), selected=0.5)),
  div(style="display:inline-block",shinyWidgets::sliderTextInput(inputId = "isolate_prob",
              label = "Probability a person with symptomatic infection will self-observe and limit contacts",
              choices = c(0.1,0.5,0.9),selected=0.1)),
  div(style="display:inline-block",shinyWidgets::sliderTextInput(inputId = "R0",
              label = "R0",
              choices = c(2,2.2,2.6,3),selected=2)),
  div(style="display:inline-block",shinyWidgets::sliderTextInput(inputId = "R0_reduce",
              label = "R0 for asymptomatic infection relative to symptomatic infection",
              choices = c(0.2,0.6,1),selected = 0.2)), 
  div(style="display:inline-block",shinyWidgets::sliderTextInput(inputId = "pre_symp",
              label = "Proportion of transmission that occurs presymptomatically",
              choices = c(0,0.1,0.3,0.5),selected=0.1)),
  
  mainPanel(position="right",
            tableOutput('table2'),plotlyOutput("graph2"),plotlyOutput("graph3"),tableOutput('table1')
  )
)

server <- function(input, output) {
  				
  output$table1 <- renderTable({
    names(nCoV_table) <- c("# Introductions",	"Proportion presymptomatic transmission",
    "Average cumulative infections end of February", "Average cumulative infections end of March", "Average cumulative infections end of April",
    "Median serial interval (days)", "Median doubling time (days)")
     nCoV_table
  }, caption="The above table shows results for simulations with 50% infections symptomatic, 10% probability of self-observation and limiting contacts, overall R0 of 2, and 0.2 relative R0 for each asymptomatic infections.")
  
  output$table2 <- renderTable({
    isolate_prob <- input$isolate_prob
    symp_prob <- input$symp_prob
    num_introductions <- input$num_introductions
    pre_symp <- input$pre_symp
    R0_reduce <- input$R0_reduce
    R0_symp <- input$R0
    
    data <- read.csv(paste0('./simulations/Info/SI_',R0_symp,"_",R0_reduce,"_",isolate_prob,"_",symp_prob,"_1.5_",num_introductions,"_",pre_symp,".csv"),header=FALSE)
    names(data) <- c("X1","SI","Doubling_times","Growth_rates","total_cases",
                   "isolate_prob","symp_prob","num_introductions","pre_symp","R0_reduce","R0")
    data$Doubling_times[data$total_cases==0] <- NA
    median <- apply(data,2,median,na.rm=TRUE)
    #mean <- apply(data,2,mean,na.rm=TRUE)
    #data <- cbind(median,mean)
    median <- as.data.frame(rbind(round(median[2],2),round(median[3],2)))
    names(median) <- "Median"
    median[median==Inf] <- NA
    median[median>30] <- ">30"
    Values <- as.data.frame(rbind("serial interval (days)","doubling time (days)"))
    names(Values) <- c("Measure")
    data <- cbind(Values,median)
  })
  
  output$graph2 <- renderPlotly({
    isolate_prob <- input$isolate_prob
    symp_prob <- input$symp_prob
    num_introductions <- input$num_introductions
    pre_symp <- input$pre_symp
    R0_reduce <- input$R0_reduce
    R0_symp <- input$R0
    Runs <- read.csv(paste0('./simulations/Runs/epidemic_curve_master_',R0_symp,"_",R0_reduce,"_",isolate_prob,"_",symp_prob,"_1.5_",num_introductions,"_",pre_symp,".csv"),header=FALSE)
    names(Runs) <- c("X","Day_infected","Sim","Cases","Imports","Severe_cases","Cum_cases","Cum_severe_cases", "Average_cases", "Average_cumulative_cases",
                      "Average_severe_cases","Average_cumulative_severe_cases", "isolate_prob","symp_prob","num_introductions","pre_symp","R0_reduce","R0_symp","date")
    Runs$date <- as.Date(as.character(Runs$date))
    Runs$Imports[Runs$Imports==0] <- NA
    Runs$Imports[!is.na(Runs$Imports)] <- 0
    
    p <- ggplot(Runs) + 
    geom_point(aes(x=date,y=Imports,text=paste0("Introductions on ",date)),color="light blue") + 
    geom_line(aes(x=date,y=Cases,group=factor(Sim),text=paste0("New cases on ",date,": ",round(Cases))),color="grey") +
    geom_line(aes(x=date,y=Average_cases,group=1,text=paste0("Average # new cases on ",date,": ",round(Average_cases))),color="black") +
    #geom_line(aes(x=date,y=Severe_cases),color="tomato1") + 
    xlab("Date of infection") + ylab("New Infections") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position = "none",
    axis.title=element_text(size=8),plot.title = element_text(size=10)) +
    ggtitle("Estimated number of new SARS-CoV-2 infections \nthrough April 2020, by infection date") 
    ggplotly(p, tooltip = "text")
  })
  
  output$graph3 <- renderPlotly({

    isolate_prob <- input$isolate_prob
    symp_prob <- input$symp_prob
    num_introductions <- input$num_introductions
    pre_symp <- input$pre_symp
    R0_reduce <- input$R0_reduce
    R0_symp <- input$R0
    Runs <- read.csv(paste0('./simulations/Runs/epidemic_curve_master_',R0_symp,"_",R0_reduce,"_",isolate_prob,"_",symp_prob,"_1.5_",num_introductions,"_",pre_symp,".csv"),header=FALSE)
    names(Runs) <- c("X","Day_infected","Sim","Cases","Imports","Severe_cases","Cum_cases","Cum_severe_cases", "Average_cases", "Average_cumulative_cases",
                     "Average_severe_cases","Average_cumulative_severe_cases", "isolate_prob","symp_prob","num_introductions","pre_symp","R0_reduce","R0_symp","date")
    Runs$date <- as.Date(as.character(Runs$date))
    Runs$Imports[Runs$Imports==0] <- NA
    Runs$Imports[!is.na(Runs$Imports)] <- 0
    
    p <- ggplot(Runs) + 
    geom_line(aes(x=date,y=Cum_cases,group=factor(Sim),text=paste0("Cumulative cases by ",date,": ",round(Cum_cases))),color="grey") +
    geom_point(aes(x=date,y=Imports),color="light blue") + 
    geom_line(aes(x=date,y=Average_cumulative_severe_cases,group=2,text=paste0("Average # cumulative severe cases by ",date,": ",round(Average_cumulative_severe_cases))),color="tomato1") + 
    geom_line(aes(x=date,y=Average_cumulative_cases,group=1,text=paste0("Average # cumulative cases by ",date,": ",round(Average_cumulative_cases))),color="black") +
    xlab("Date of infection") + ylab("Total Infections") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position = "none",
          axis.title=element_text(size=8),plot.title = element_text(size=10)) +
    ggtitle("Estimated number of cuumulative SARS-CoV-2 infections \nthrough April 2020, by infection date") 
    ggplotly(p, tooltip = "text")
      
  })
 

}

shinyApp(ui = ui, server = server)



