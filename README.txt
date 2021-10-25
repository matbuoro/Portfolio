MetaIBASAM:  A demo-genetic agent-based model to simulate spatially structured salmon populations

MetaIbasam is an extension of the existing IBASAM model (https://github.com/Ibasam/IBASAM/wiki) by incorporating a dispersal process to describe Atlantic salmon metapopulation and its eco-evolutionary dynamics. MetaIBASAM allows an investigation of the consequences of dispersal on local populations and network dynamics at the demographic, phenotypic, and genotypic levels. More generally, it allows to explore eco-evolutionary dynamics by taking into account complex interactions between ecological and evolutionary processes (plasticity, genetic adaptation and dispersal), feedbacks (e.g. genetic <-> demography) and trade-offs (e.g. growth vs survival). By doing so, one can investigate responses to changing environments and alternative management strategies.
----

Contacts: 
amaia.lamarins@inrae.fr 
mathieu.buoro@inrae.fr
-



0. Define the directory in a terminal

> cd Portfolio/

1. To launch simulations in a terminal:

> Rscript --vanilla metaIbasam.R 4 2 3 1 30 &

Arguments:
#1: scenarioConnect (e.g. 4)
#2: scenarioFishing (e.g. 2)
#3: scenarioFishingrates (e.g. 3)
#4: first simulation (e.g. 1)
#5: last simulation (e.g. 30)

Each scenario (combination of arguments 1, 2 and 3) has its own temporary folder and results folder. Several simulations can run in parallel, depending on the computer memory.


2. To launch simulations at a specific hour (using the linux package "at"):
# at: could be today, tomorrow,... friday,...
> at 04:00am tomorrow
> pkill -9 -u mbuoro R
> Rscript --vanilla metaIbasam.R 4 2 1 1 30 &
> Rscript --vanilla metaIbasam.R 7 1 2 1 30 &
> Rscript --vanilla metaIbasam.R 7 2 2 1 30 &

then crtl+D



3. To extract results (e.g. Demography.R):

> cd code/
> nohup bash Rextract.sh &


The Rextract.sh file launches the file "DEMOGRAPHY.R" on each scenario defined in Rextract.sh.


 

 
 
 