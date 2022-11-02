# HL-LHC-extrapolation
Automatic python scripts to execute the extrapolation from an imput datacard to high-luminosity LHC. 

The whole program consists in two main parts:

FIRST PART --> mkTGraph.py
python mkTGraphs.py <datacard-name.txt> --lumi-step <integer>

Datacards creation and submission to condor: from the basic datacard, the program builds other datacards exploiting the Condor versatility and using the information from:

	– the luminosity values chosen in the range [150/fb,4500/fb];
	– the selected scenarios, i.e. reading and interpreting fitted plot on ECAL energy resolution and electron efficiency. The systematic uncertainties involved are CMS_e_scale_2018, CMS_met_scale_2018, CMS_e_eff_2018, and the ones related to the JES. 

The program creates for each datacard a workspace with the specified features and submit it to Condor (it takes „ one hour and an half). Condor will create root files which will contain, for each luminosity of a particular scenario, a likelihood scan with respect to the signal strength parameter. For each likelihood scan 1500 points are set in the multidimensional fit.


SECOND PART --> mkPlots.py
python mkPlots.py --compare 1

Plot maker: this script takes as input the root files containing the likelihood scans and - using the basic rule expressed by the Rao-Cramer theorem - get the signal strength uncertainty. It stores the values obtained in a matrix and plots them in a ROOT canvas.


mkImpacts.py --> the script allows the creation of impact plots, starting from .root files.
Main commands to execute it are reported in the script itself.
