# CODESIR
CODESIR contains a basic ODE SIR model and variations.

codeSIR.py contains the most basic SIR model using the Euler method. It outputs a plot with susceptible, infected, recovered, and total populations represented.

Version 1 of codeSIRdiscreteq.py contains the SIR model that allows the user to input varying "q factors" that represent the effectiveness of a population's non-pharmaceutical interventions for an initial period of time, followed by a period of no interventions, which can be followed by a period of resumed interventions at whatever level the user inputs. It outputs a plot with susceptible, infected, recovered, and total populations represented.

Version 2 of the discreteq variant of our SIR model now only plots the number of infected people, allows the user to also vary the level of relaxation during the period of no interventions, and allows the user to input a "comparison spacing" parameter that allows the code to plot multiple time periods of relaxation on one plot. It outputs a plot with no relaxation period, 1*"spacing parameter" length relaxation period, all the way to 4*"spacing parameter" length relaxation period to compare the effects of turning on and off certain levels of npi's.

Version 3 of the discreteq variant of our SIR model has removed the user input aspect of the previous versions, and replaced this by defining a function called "quarantine_comparison" that allows the user to keep some of the parameters fixed while varying others without having to repeatedly input those values. It currently outputs 15 graphs that each have a different start date to the lifting of the quarantine. Each graph contains the same 5 comparison plots from the previous version.

Version 4 of the discreteq variant of the SIR model has an improved title, and now includes a curve where no interventions took place as a point of comparison for the other possibilities.

Version 5 of the discreteq variant of the SIR model no longer allows quarantine to be turned back on. It just looks at when the best day to remove restrictions is. It outputs a plot with 6 different end days for restrictions plus a curve with no restrictions and a curve with no end to the restrictions. Currently it produces 10 graphs the vary the "q factor" to look to how that effects when we can relax restrictions.

Version 1 of codesircombmixdisq.py is essentially the same thing as version 5 of discreteq, but now we have the ability to control how much of the population is participating in the restrictions. In this variation, one group does not participate the entire time, and the other initially participates, then stops after a given amount of time. The "exp" variation has both groups initially following restrictions. Then one group stops while the other keeps following the restrictions. Both only output 1 graph with 8 curves currently.

Version 2 of combmixdisq is essentially the same thing as version 4 of discreteq, but now we have the ability to control how much of the population is participating in the restrictions. In this variation, one group does not participate the entire time, and the other initially participates, then stops after a given amount of time, and then resumes participation in restrictions. The "exp" variation has both groups initially quarantine (possibily at different values), then the restrictions are lifted, and then one group refuses to go back to quarantine while the other does go back.

Any of the combmixdisq models can be easily modified to look at differing population dynamics.

Note that codesirqformon.py and codesirqformonb.py are still very much works in progress and still have some issues that need to be worked out.
