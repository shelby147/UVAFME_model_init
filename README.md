# UVA Forest Model Enhanced: Initialization version
This repository holds the source code for UVAFME and pre- and post-processing scripts for figures and data in Sundquist _et al._ (2024). [https://dx.doi.org/10.1088/2752-664X/ad7d94](https://dx.doi.org/10.1088/2752-664X/ad7d94)

UVAFME is an individual-based forest gap model which simulates the establishment, growth, and mortality of individual trees on independent gaps or "plots" of a forested landscape.

For more information about UVAFME visit the [UVAFME website](https://uvafme.github.io/) and GitHub repository for the [previous version](https://github.com/UVAFME/UVAFME_model/) of the model.

Please contact [Shelby Sundquist](https://orcid.org/0000-0001-5379-0008) at ss3988@nau.edu with any questions. 

Model input and output data and raw and intermediate data products are available upon request. 

## Updates

_March 2025_

Src code updates to implement dynamic management, where UVAFME allocates forest treatments across a landscape according to a vector of likelihoods for different management actions.

Scripts for inputs and analysis pertaining to a new manuscript, "Individual-tree based modeling climate-adaptive silviculture strategies to inform policy in Alaska’s boreal forest".

## Repository folder structure

_Summaries are available in select subdirectories in Description.txt_

scripts - creating input files for UVAFME and analyzing output

scripts/*project* - inputs and analyses related to specific publications 

scripts/*/inputs - creating input files for TVSF and KLC study areas

scripts/*/analyses - creating figures from model output in TVSF and KLC

UVAFME_src/src - source code for version of UVAFME model used in this study
