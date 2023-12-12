/*
 * File: TenTusscher2IschemiaSetup.h
 *
 * Institute of Biomedical Engineering, 
 * Karlsruhe Institute of Technology (KIT)
 * https://www.ibt.kit.edu
 * 
 * Repository: https://github.com/KIT-IBT/CardioMechanics
 *
 * License: GPL-3.0 (See accompanying file LICENSE or visit https://www.gnu.org/licenses/gpl-3.0.html)
 *
 */



#ifndef TENTUSSCHER2ISCHEMIASETUP
#define TENTUSSCHER2ISCHEMIASETUP

// ************************************************************************************************************
// ************************************************************************************************************
// ******  USER SPECIFIED SETUP FOR THE SETTINGS OF ISCHEMIA !!!
// ************************************************************************************************************
// ************************************************************************************************************

//// if defined, the newly added ATP channel considering the ADP and ATP concentration is activated
// #define ACTIVATE_IKATP_CHANNEL

//// if defined, the time dependend changing and the different effects during ischemia are activated
//// this only works and makes sense if the ACTIVATE_IKATP_CHANNEL option above is defined
// #define ISCHEMIA

// here you can set the type of ischemia that should be considered
// this only works and makes sense if the ISCHEMIA option above is defined
#ifdef ISCHEMIA
# define HYPERKALEMIA
# define ACIDOSIS
# define HYPOXIA

// #define ISCHEMIA_VERBOSE_OUPUT
#endif  // ifdef ISCHEMIA

// ************************************************************************************************************
// ************************************************************************************************************
// ******  END OF: USER SPECIFIED SETUP FOR THE SETTINGS OF ISCHEMIA !!!
// ************************************************************************************************************
// ************************************************************************************************************


// ************************************************************************************************************
// ************************************************************************************************************
// ******  NO USER SPECIFIED SETUP - LEAVE THE FOLLOWING LINES UNCHANGED !!!
// ************************************************************************************************************
// ************************************************************************************************************
// do not change the following lines! You can set up the ischemia effects in in the lines above!
// the following options are used to (de)activate the time course of different variables used for ischemia
#ifdef ISCHEMIA
# ifdef HYPERKALEMIA
        #  define KO
# endif  // ifdef HYPERKALEMIA
# ifdef ACIDOSIS
        #  define GCAL
        #  define GNA
        #  define dVmNa
# endif  // ifdef ACIDOSIS
# ifdef HYPOXIA
        #  define ATP
# endif  // ifdef HYPOXIA

// #define MGI
#endif  // ifdef ISCHEMIA

// do not change the lines above! You can set up the ischemia effects in the lines above!
// ************************************************************************************************************
// ************************************************************************************************************
// ************************************************************************************************************


#endif  // ifndef TENTUSSCHER2ISCHEMIASETUP


/*-------------------------------------------------------

   Changelog:

 # 23.07.2008 - 1.0.0 (Daniel Weiss)
     Initial release
 */
