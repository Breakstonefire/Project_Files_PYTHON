In this folder you will find differents versions of my model about ploting the concentration of lin4, lin14 and LIN14 protein but also some codes about the activator of this gene network in itself.

For every codes, I will describe what do we plot.
------------------------------------------------------------------------------------------
defining_the_activator:
This code is used to compute the initiation rate of the activator of lin4 gene which is composed of NHR23 + NHR85 + hormones + LIN42.
We define all the constants and parameters in order to plot the resulting shape of the activator.
In this one especially I let unset the parameters that are unknown but it is commented so it still can be used to modify a bit hte activator (at least the parameters).

------------------------------------------------------------------------------------------
Fiting_experimental_results_for_NHR23-NHR85-LIN42_parameters_research we plot:
- the NHR23 fit to DNA-binding results from Chris Hammel's lab (in fact it is the concentration directly but the idea remains)
- the NHR85 fit to DNA-binding results from Chris Hammel's lab (in fact it is the concentration directly but the idea remains)
- the LIN42 fit to DNA-binding results from Chris Hammel's lab (in fact it is the concentration directly but the idea remains)

------------------------------------------------------------------------------------------
Fiting_NHR23-NHR85-LIN42_v2_One_Sine_Wave_of_Frequency_1:
In this code we do the same as the previous one but here the sine wave frequency of the fit is a constant an is 1 because we now assume that the concentration oscillates in every larval stage at a frequency = 1).
- the NHR23 fit to DNA-binding results from Chris Hammel's lab (in fact it is the concentration directly but the idea remains)
- the NHR85 fit to DNA-binding results from Chris Hammel's lab (in fact it is the concentration directly but the idea remains)
- the LIN42 fit to DNA-binding results from Chris Hammel's lab (in fact it is the concentration directly but the idea remains)
- the Supposed {NHR23;NHR85;LIN42} concentrations during one Larval Stage

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm:
- the concentration of the activator in the nucleus in time (4 larval stages) considered as a whole molecule
- the concentration of microRNA lin4 in the nucleus in time (4 larval stages)
- the concentration of mRNA lin14 in the nucleus in time (4 larval stages)
- the concentration of protein LIN14 in the cytoplasm in time (4 larval stages)


------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v2:
- the concentration of the activator in the nucleus in time (1 larval stages) considered as a whole molecule
- the concentration of the activator in the nucleus in time (3 larval stages) considered as a whole molecule
- the concentration of microRNA lin4 in the nucleus in time (3 larval stages)
- the concentration of mRNA lin14 in the nucleus in time (3 larval stages)
- the concentration of protein LIN14 in the cytoplasm in time (3 larval stages)

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v3_threshold:
- the concentration of the activator in the nucleus in time (1 larval stages) considered as a whole molecule
	Here we added the threshold on which is based this different model: when the activator is upper than this threshold, the activator acts.
- the concentration of the activator in the nucleus in time (3 larval stages) considered as a whole molecule
- the concentration of microRNA lin4 in the nucleus in time (3 larval stages)
- the concentration of mRNA lin14 in the nucleus in time (3 larval stages)
- the concentration of protein LIN14 in the cytoplasm in time (3 larval stages)
- the parameters that are chosen in the begining of the code are showed at the end also

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v4_get_rid_of_the_complex_and_time_scaled:
In this code we do not consider the presence of an intermediate complex at all formed by lin4 binded on lin14 compared to the previous models.
- the concentration of mRNA lin14 in the nucleus in time (4 larval stages)
- the concentration of protein LIN14 in the cytoplasm in time (4 larval stages)

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v5.0_Normalized_Scaled_Good_colors:
- the concentration of the activator in the nucleus in time (1 larval stages) considered as a function of 2 molecules (NHR23 and NHR85) that is active when the two of them are present in the nucleus only
- the concentration of the activator in the nucleus in time (4 larval stages) considered as a whole molecule
- the concentration of microRNA lin4 in the nucleus in time (4 larval stages)
	Here the numerical solution and the analytical one are shown.
- the concentration of mRNA lin14 in the nucleus in time (4 larval stages)
- the concentration of protein LIN14 in the cytoplasm in time (4 larval stages)

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v5.1_Timings_of_Transcription:
This one is the same one as the next one but it contains overflows due to exponential in the functions that describe the activator.

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v5.2_without_taking_LIN42_effect:
Here the overflows are corrected but the model is still wrong, the dependency on time in the Initiation rate function is not the right one.
The new model with the right way to take time into account is detailed in the report and in the last codes that I modified after I left the lab.

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v5.3_taking_LIN42_effect:
- the NHR23 Pattern and the NHR85 Pattern and the LIN42 pattern during one larval stage
- the occupancies of NHR23 and NHR85 and the hormones in the nucleus
	Note that the occupancy of the hormone is always the same because we consider one cell here, so according to the diffusion models, after one minute, the concentration of hormones is the constant at a given space of the worm.
- the product of the three previous occupancies : Occup.NHR23*Occup.NHR85*Occup.Hm

- the resulting initiation rate "Irate of lin4 gene" in time during the 4 larval stages and the resulting lin4 concentration in time (all are normalized)
- the resulting initiation rate "Irate of lin14 gene" in time during the 4 larval stages and the resulting lin14 concentration 
- the concentration of protein LIN14 in the cytoplasm in time (4 larval stages)

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v6_Gillespie_Normalized_Scaled_Good_colors
This model is not working because I hadn't the time to add gillespie randomness in the model.

------------------------------------------------------------------------------------------
Modeling_mRNA_lin-14_algorithm_v7_Parameters_Range_Searching
This model is not working for time restrictions of my internship but it would have been nice to add some loops in the code to search the parameters that could give the expected results (the expected shapes for lin4 and lin14 and LIN14 curves during all the development).
