Code for scenario exploration for SARS-CoV-2 infections

Developed by the Center for Communicable Disease Dynamics at the Harvard T.H. Chan School of Public Health (Rebecca Kahn, Stephen Kissler, Nishant Kishore and Marc Lipsitch), in collaboration with Health Department colleagues. This is for scenario exploration only. Some scenarios generated may be unrealistic, and users should exercise judgment in interpreting the outputs.

Overview:

- We calibrated timing of introductions of infected persons into a city based on an epidemic curve with exponential growth beginning on December 1 (based on genomics data, SARS-CoV-2 was estimated to have been first introduced into the human population sometime in November or early December 2019) and then leveling off on January 23 (when travel restrictions were implemented in Wuhan) to account for decreased travel from Wuhan but to allow for continued introductions from elsewhere. Light blue dots along the x-axis show dates of introductions. We then modeled the subsequent transmission chains stemming from those introductions.

- We varied parameters for number of expected introductions, proportion of infections that are symptomatic, probability a person with a symptomatic infection will self-observe and limit contacts, R0, R0 for asymptomatic infection relative to symptomatic infection, and proportion of transmission that occurs presymptomatically. Change the inputs on the bars below to see how these affect the outbreak trajectories.

- The y axis shows the number of either new or cumulative infections in a city, which are a combination of imported infections and their subsequent infections. Results are shown for 5 months.

The black line is the average of the 50 stochastic simulations shown in grey, and the red line is the average number of severe COVID-2019 illnesses. All 50 simulations shown on a given graph have the same parameters. The variation between outbreak trajectories on the same graph comes from the stochasticity in the following:

- Timing of introductions (the number of infected people introduced on a given day is drawn from a binomial distribution, with the probability of introduction calibrated to the exponential epidemic curve).

- Number of secondary infections (drawn from a negative binomial distribution with mean of R0 and dispersion parameter 0.5 (Riou 2020) for each symptomatic infection and a mean = R0*reduction in R0 for asymptomatic infection). This distribution allows for superspreaders, which can impact early epidemic dynamics.

- Timing of secondary infections (drawn from a triangular distribution, calibrated to the percent of transmission that occurs presymptomatically).

- Whether an infected person will have symptoms (drawn from a binomial distribution for each infected person with probability based on the input value selected). Among infections with symptoms, whether or not the infection is severe is drawn from a binomial distribution with probability 0.18 (WHO).

- Incubation period for symptomatic infections (drawn from a Weibull distribution with mean 5.8 days for each illness (Backer 2020)).

- Whether a person with symptomatic infection will self-observe and limit contacts (drawn from a binomial distribution for each illness with probability based on the input value selected).

- The time at which they will start to limit contacts and continue to do so for the duration of their infection (drawn from a gamma distribution with mean 1.5 days after symptoms).

Important caveats:

- Onward infections stop if 100,000 people are infected. This branching model only includes infected individuals so does not incorporate depletion of susceptibles or contact patterns. With the exception of potential self-observation and limiting of contacts by symptomatic individuals, it also does not incorporate any control measures.

- In principle, R0 depends on population density, but there is little evidence that it actually does. Use the slider to change the value of R0.

- We are working on incorporating flight data to better calibrate the timing of introductions.

- We have not incorporated any changes in R0 due to seasonality given the limited data and uncertainty. If warmer weather will decrease R0, the total number of infections by the end of April will be overestimated.

- The table at the top shows the median serial interval (time between symptom onset for infector-infectee pair) and doubling time. Note, doubling times are only calculated once at least 5 local infections have occurred and only for simulations with a total number of infections by the end of April at least 3 times the number of introductions. Doubling times between 5-7 days are most consistent with the Wuhan outbreak (Wu 2020), so doubling times shorter than this likely indicate unrealistic parameter combinations. Longer doubling times may reflect impact of self-observation.

- The static table at the bottom (does not change with inputs) shows the average cumulative infections by the end of February, March and April across values for number of introductions and proportion of presymptomatic transmission for simulations with 50% infections symptomatic, 10% probability of self-observation and limiting contacts, R0 of 2.2 for symptomatic infections, and 0.2 relative R0 for asymptomatic infections.

Findings:

- The amount of presymptomatic transmission and the extent to which infected people limit contacts (a combination of proportion symptomatic and probability they will self-observe) have a large impact on the total number of infections.

To explore scenarios: https://rebeccakahn.shinyapps.io/COVID19/

References

- Introduction time: http://virological.org/t/phylogenetic-analysis-of-23-ncov-2019-genomes-2020-01-23/335/5
  
- Incubation period: https://www.medrxiv.org/content/10.1101/2020.01.27.20018986v1.full.pdf

- Severity: https://www.who.int/docs/default-source/coronaviruse/transcripts/transcript-coronavirus-press-conference-full-07feb2020-final.pdf?sfvrsn=3beba1c0_2

- Dispersion parameter: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.4.2000058

- Doubling time: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30260-9/fulltext
