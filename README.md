# InsightPRSConstruction
This code provides functions used in the analysis for the paper "Polygenic risk score based on children’s growth curves is a strong predictor of rapid infant weight gain linked to obesity later in life".

This work is maintained by Ana Maria Kenney at Penn State, questions about the code should be sent to her (akenney375@gmail.com).  

### Files
The main codes are:
 * ComputePRS.R - functions to construct our proposed PRS (polygenic risk score) and evaluate its effect on growth curves
 * ProcessGrowthCurves.R - functions to convert discretely observed growth curves to functional objects and switch to "shifted" curves as described in the main paper
 * ScreeningApproach.R - an implementation of the screening method in the paper "Feature Screening for Time-Varying Coefficient Models with Ultrahigh Dimensional Longitudinal Data"
 * RunFlame.R - a function that runs FLAME
 * flm_0.1.tar - zip file to install FLAME, see [link](http://www.personal.psu.edu/mlr36/codes.html) for more information and examples
 * FS_penreg_v2.R - function to perform function-on-scalar regression used to evaluate PRS
