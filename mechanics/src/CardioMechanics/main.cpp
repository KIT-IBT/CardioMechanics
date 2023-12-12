/*
 * File: main.cpp
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
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include "DCCtrl.h"
#include "CardioMechanics.h"
#include "filesystem.h"
#include <kaExceptions.h>

using namespace math_pack;

void PrintCardioMechanicsBanner() {
    DCCtrl::print << "\n";
    DCCtrl::print << "--------------------------------------------------------------------------------\n";
    DCCtrl::print << " C a r d i o M e c h a n i c s ";
    DCCtrl::print << " [ build: " << __DATE__ ;
    DCCtrl::print << " ]\n\n";
    DCCtrl::print << " Institute of Biomedical Engineering\n";
    DCCtrl::print << " Karlsruhe Institute of Technology (KIT)\n";
    DCCtrl::print << " Kaiserstr. 12, 76131, Karlsruhe, Germany\n";
    DCCtrl::print << "--------------------------------------------------------------------------------\n\n";
}

void PrintHelpText() {
    std::stringstream options;
    
    options << left << setw(25) << "-help" << "Show this help dialog\n"
    << left << setw(25) << "-settings <xml>" << "Setting File\n"
    << left << setw(25) << "-check_settings <xml>" << "Check settings file\n"
    << left << setw(25) << "-no_try_catch" << "Runs simulation without a try/catch loop\n"
    << left << setw(25) << "-handbrake" << "Wait for keypress before simulation starts, e.g. to attach to debugger\n"
    << left << setw(25) << "-verbose" << "Tell me all your secrets\n"
    << left << setw(25) << "-debug" << "Tell me all your secret secrets\n";
    
    DCCtrl::print << options.str() << std::endl;
}


void ShowImageYouAreStupid()
{
    // Get user name to insult him
    char *p=getenv("LOGNAME");
    std::stringstream cmd;
    cmd << "finger " << p;
    FILE* pipe = popen(cmd.str().c_str(),"r");
    char line[128];
    fgets( line, sizeof(line), pipe);
    std::stringstream ss(line);
    std::string name;
    for(int i=0; i < 4; i++)
        ss >> name;
    
    DCCtrl::print << "\n\n\n\t\t\t\t         \\|||/\n";
    DCCtrl::print << "\t\t\t\t         (o o) \n";
    DCCtrl::print << "\t\t\t\t ,~~~ooO~~(_)~~~~~~~~~,\n";
    DCCtrl::print << "\t\t\t\t | Shame on you,      |\n";
    DCCtrl::print << "\t\t\t\t | " << setw(19) << std::left << name << "|\n";
    DCCtrl::print << "\t\t\t\t | That's too stupid! |\n";
    DCCtrl::print << "\t\t\t\t '~~~~~~~~~~~~~~ooO~~~'\n";
    DCCtrl::print << "\t\t\t\t        |__|__|\n";
    DCCtrl::print << "\t\t\t\t         || ||\n";
    DCCtrl::print << "\t\t\t\t        ooO Ooo\n";
    
}

CBStatus InitAndRunSimulation(CardioMechanics& cardio, bool shouldCheckSettingsOnly) {
    
    cardio.Init1();
    cardio.Init2();
    
    if (shouldCheckSettingsOnly) {
        cardio.PrintParameters();
        DCCtrl::print << "Settings seems to be okay, but check with parameters were requested and ignored during initialization !!!!\n";
        DCCtrl::print << "It is possible that some of the ignored parameters would have been requested during the simulation (e.g. Export.Options.XXX) !!!!\n\n\n\n";
        return CBStatus::SUCCESS;
    }
    
    cardio.Run();
    
    return CBStatus::SUCCESS;
}

void Finalize(bool exitCode) {
    DCCtrl::Finalize();
    exit(exitCode);
}

/// returns false if CardioMechanics should stop after this function
bool DigestCommandline(int argc, std::string parameterFile, bool &shouldCheckSettingsOnly,
                       char* argv[], CardioMechanics& cardio, bool &shouldTryAndCatch) {
    std::vector<std::string> extraParameters;
    int i = 1;
    while (i < argc) {
        if (std::string(argv[i]) == "-handbrake") {
            if (DCCtrl::IsProcessZero()) {
                DCCtrl::print << "Press enter to continue ... " << std::endl;
                std::string dummy;
                getline(std::cin, dummy);
            }
            MPI_Barrier(PETSC_COMM_WORLD);
        }
        if (std::string(argv[i]) == "-help") {
            PrintHelpText();
            return false;
        }
        if (std::string(argv[i]) == "-settings" && i+1 != argc)
            parameterFile = std::string(argv[++i]);
        
        if (std::string(argv[i]) == "-check_settings" && i + 1 != argc) {
            parameterFile = std::string(argv[++i]);
            shouldCheckSettingsOnly = true;
            DCCtrl::print << "\n";
            DCCtrl::print << "         Checking simulation setting ...\n";
            DCCtrl::VerboseOn();
        }
        
        if (std::string(argv[i]) == "-no_try_catch")
            shouldTryAndCatch = false;
        
        if (std::string(argv[i]) == "-verbose")
            DCCtrl::VerboseOn();
        
        if (std::string(argv[i]) == "-debug")
            DCCtrl::DebugOn();
        
        if (std::string(argv[i]) == "-parameter" && i+1 != argc) {
            std::string s = argv[++i];
            extraParameters.push_back(s);
        }
        i++; // iterator over variables
    } // end while
    
    if(parameterFile == "")
    {
        PrintHelpText();
        return false; // back to main, stop after this
    }
    
    cardio.ReadParameterFile(parameterFile);
    for (auto ep:extraParameters) {
        std::string key   = "";
        std::string value = "";
        int pos = ep.find('=');
        if (pos!=-1) {
            key   = ep.substr(0, pos);
            value = ep.substr(pos+1, ep.length());
            cardio.SetParameter(key, value);
        } else {
            DCCtrl::print << "Unrecognized option from commandline found: " << ep << "\n";
        }
    }
    
    return true;
}

int main(int argc, char* argv[])
{
    std::string    parameterFile("");
    
    bool shouldTryAndCatch       = true;
    bool shouldCheckSettingsOnly = false;
    
    DCCtrl::Init(argc, argv);
    
    PrintCardioMechanicsBanner();
    
    CardioMechanics cardio;
    
    
    //===== digest commandline parameters =====//
    try {
        if(sizeof(double) != 8 || sizeof(int) != 4)
            throw std::runtime_error("I expected sizeof(double) to be 8 and sizeof(int) to be 4 but this doesn't seem to be the case, this might get us in trouble, I don't say it will, but it might...");
        
        bool continueToSimulation = DigestCommandline(argc, parameterFile, shouldCheckSettingsOnly, argv, cardio, shouldTryAndCatch);
        
        if (!continueToSimulation) {
            DCCtrl::print << "Stop requested.\n";
            Finalize(0);
        }
    }
    catch(const char* e) {
        if(DCCtrl::IsProcessZero())
            cerr << "\nError during initialization: " << e << "\n";
        Finalize(1);
    }
    catch(std::exception& e) {
        if(DCCtrl::IsProcessZero())
            cerr << "\nError during initialization: " << e.what() << "\n";
        Finalize(1);
    }
    catch(...) {
        if(DCCtrl::IsProcessZero())
            cerr << "\nError during initialization: Exception of unknown type!\n";
        Finalize(1);
    }
    
    
    //===== run simulation =====//
    
    if(shouldTryAndCatch || shouldCheckSettingsOnly)
    {
        try {
            InitAndRunSimulation(cardio, shouldCheckSettingsOnly);
        }
        catch(std::exception& e) {
            if(DCCtrl::IsProcessZero())
            {
                ShowImageYouAreStupid();
                cerr << "\n\tRuntime error: " << e.what() << "\n";
            }
            Finalize(1);
        }
        
        catch(kaBaseException& e) {
            if(DCCtrl::IsProcessZero())
            {
                ShowImageYouAreStupid();
                cerr << "\n\tRuntime error: " << e << "\n";
            }
            Finalize(1);
        }
        
        catch(...)
        {
            if(DCCtrl::IsProcessZero())
            {
                ShowImageYouAreStupid();
                cerr << "\n\tRuntime error with unkown exception type: " << std::endl;
            }
            Finalize(1);
        }
    }
    else
    {
        InitAndRunSimulation(cardio, shouldCheckSettingsOnly);
    }
    cardio.PrintParameters();
    DCCtrl::print << "\n\n\nLooking forward to the next job, Semper fi\n\n";
    cardio.DeInit(); // this is a workaround because cardio's destructor does not get called
    Finalize(0);
}
