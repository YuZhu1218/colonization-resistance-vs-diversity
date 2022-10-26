# Microbial Consumer Resource Model

## Description 
The microbial consumer resource model (MiCRM) that we adapted are in the paper "Niehaus et al (2019). Microbial coexistence through chemical-mediated interactions. Nature Communications. https://doi.org/10.1038/S41467-019-10062-X." 

## File

    CRM_mainscript.m --contains the main script that simulating the micrbial community and introduce an invader after the community is stable, return the composition of species before and after invasion
    BinomialSampling.m --the function that returns a matrix of certain size where the each cell is a bionomial probability with defined p.
    DirichletSampling.m --the function that returns a matrix of nMediator * nMediator size that is sampled from Dirichlet distribution.
    FamilyEncounter.m --the function that assigns the species and resources to either family A or S with certain fraction of species and resources in the each family
    SupplyResource.m --the function that returns the initial fraction of each resources.
    WellmixedInteraction.m --the mixing and dynamic progress of microbial community.


