# Mediator-explicit model

## Description 
Mediator-explicit model is firstly introduced in the paper "Niehaus et al (2019). Microbial coexistence through chemical-mediated interactions. Nature Communications. https://doi.org/10.1038/S41467-019-10062-X". It takes the chemicals that either produced or comsumed by the microbial species as the mediator and the interactions between species are indirectly represented by the interactions between chemicals and species. For more details, take a look of the paper.

## Files

    CommunityResistanceSimulation.m --the main script that simulate 10,000 samples of microbial community and introduce an invader, returns the composition of species before and after invasion.
    DistInteractionStrengthMT_PA.m  --the function to construct the interaction matrix between species and chemicals with 50:50 ratio of positive and negative interactions.
    DistInteractionStrengthMT_PB.m --the function to construct the interaction matrix between species and chemicals with a defined fraction of positive interactions over all interacions.
    WellmixedInteraction_DpMM_ExMT4C.m --the mixing and dynamic progress of microbial community.
    NetworkConfig_Binomial.m --the function that returns the a matrix of certain size where the each cell is a bionomial probability with defined p.



