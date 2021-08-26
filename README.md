# Overcoming-dispersion-trains
This repository provides COMSOL and MATLAB model files discussed in the paper: "Method for overcoming temporal dispersion in unmyelinated nerves to image them with Electrical Impedance Tomography (EIT)", I Tarotin et al. 

List of attached files:

 - **Cfibre_model_fin_RepdZ.mph**: FEM model of a mammalian C fibre, COMSOL
 - **Cfibre_model_fin_RepdZ_cuff.mph**: FEM model of a mammalian C fibre with the cuff (see paper), COMSOL
 - **dZ_trains_FEM.mat**: MATLAB database with dZ trains simulated with the FEM model of a single C fibre with the cuff
 - **ap_dZ_shape_C_cuff.mat**: MATLAB database with single AP and dZ simulated with the FEM model with the cuff
 - **Model_matching_with_experiment.m**: MATLAB code used for adjustment of the created statistical model to the experimental recordings
 - **FullStatModel_clean.m**: MATLAB code of the developed statistical model
