# Changelog
Changes to the project are listed in this file and are listed chronologically, with the most recent entries at the top.

## Latest updates
### 2024-12-23
**Author:** @Jan-ErikDuhme 

**New:**
- OHaraRudyIso model: adds PKA phosphorylation to multiple different currents as specified in the paper by Heijman et al. 2011 (https://doi.org/10.1016/j.yjmcc.2011.02.007) to the OHaraRudy model. 
The fractions of PKA phosphorylation can be specified in the corresponding ev-File (...P_frac).

- Tomek model  (https://doi.org/10.7554/eLife.48890) with the addition of a stretch activated current and troponin C coupling to the Land17 tension model (can be removed from the model via SAC and TRPN preprocessor macros in the parameter.h file).

**Improved:**
- Added possibility to specify initial values of the Land17 tension model in the CardioMechanics config file.

**Fixed:**
- Changed PCa Multiplier in OHaraRudy mid ev-File from 1.8 to original value 2.5.

**Known issues:**

## Previous Updates

### 2024-12-19
**Author:** @MarieHouillon

**New:**

**Improved:**
- Hardcoded Docker image name has been replaced with the GitHub variable GITHUB_REPOSITORY in lowercase so that the action can be run in forks.
- Replace specific commit with semantic versioning for Docker actions.

**Fixed:**

**Known issues:** 

## How to Use This Changelog
- Entries are grouped by update date.  
- Categories help classify the type of change for quick reference.  
- Review entries to understand the project's development and evolution.