This R project gathers all the data and R code needed to reproduce simulations and results analysis that are presented in the paper untitled "Overlooked implications of dispersal in Atlantic salmon: lessons from a demo-genetic agent-based model".

It is the first application study using MetaIBASAM, a demo-genetic agent-based model to simulate spatially structured salmon populations.

MetaIbasam is an extension of the existing IBASAM model (https://github.com/Ibasam/IBASAM/wiki) by incorporating a dispersal process to describe Atlantic salmon metapopulation and its eco-evolutionary dynamics. MetaIBASAM allows an investigation of the consequences of dispersal on local populations and network dynamics at the demographic, phenotypic, and genotypic levels. More generally, it allows to explore eco-evolutionary dynamics by taking into account complex interactions between ecological and evolutionary processes (plasticity, genetic adaptation and dispersal), feedbacks (e.g. genetic <-> demography) and trade-offs (e.g. growth vs survival). By doing so, one can investigate responses to changing environments and alternative management strategies.

In this study, we ran different scenarios of dispersal rates and look at their demo-genetic consequences on local populations and the whole metapopulation.

----

Contacts: 
amaia.lamarins@inrae.fr 
mathieu.buoro@inrae.fr

----

Below is described the workflow:


0. Install metaIbasam package v.0.0.6 through the .tar.gz file and define the directory in a terminal

> cd /folder/Paper_metaIbasam_scripts


1. Run the model for a defined scenario and number of simulations - launch simulations in a terminal:

> Rscript --vanilla metaIbasam.R 1 0 3 1 30 &

Arguments:
#1: scenarioConnect (e.g. 1 for 0 dispersal)
#2: scenarioFishing (e.g. 0 for fishing all populations)
#3: scenarioFishingrates (e.g. 3 for fishing at 7% 1SW and 15% MSW)
#4: first simulation (e.g. 1)
#5: last simulation (e.g. 30)

Each scenario (combination of arguments 1, 2 and 3) has its own temporary folder and results folder. Several simulations can run in parallel, depending on the computer memory.


2. To launch simulations at a specific hour (using the linux package "at"):
# at: could be today, tomorrow,... friday,...
> at 04:00am tomorrow
> pkill -9 -u mbuoro R
> Rscript --vanilla metaIbasam.R 1 0 3 1 30 &
> Rscript --vanilla metaIbasam.R 2 0 3 1 30 &
> Rscript --vanilla metaIbasam.R 3 0 3 1 30 &

then crtl+D



3. Extract the results (e.g. DEMOGRAPHY.R, PHENOGENOTYPE.R) on defined scenarios:

> cd code/
> nohup bash Rextract.sh &


The Rextract.sh file launches the files "DEMOGRAPHY.R" and "PHENOGENOTYPE.R" on each scenario defined in Rextract.sh and number of simulations defined in the files.


4. Plot the figures found in the paper illustrating some model processes and simulations outcomes.

Run the R code in the file "FIGURES.R" in the "code" folder.

