////////////////////////////////////////////////////////////////////////////////
///
///   main.cc for blpgen
///
///
/// - Authors: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///            Antonino Cannavo (acannavo@bnl.gov      )
////////////////////////////////////////////////////////////////////////////////



///-----------------------------------------------------------------------------
/// Headers
///-----------------------------------------------------------------------------
#include <unistd.h>
#include <string>
#include <iostream>

//#include "ConMan.hh"
#include "BeamGen.hh"
#include "IPGen.hh"
#include "EvtGen.hh"



///-----------------------------------------------------------------------------
/// Global variables
///-----------------------------------------------------------------------------
//ConMan* gCM = 0;



///-----------------------------------------------------------------------------
///  Declaration of PrintHelp()
///-----------------------------------------------------------------------------
void PrintHelp();



///-----------------------------------------------------------------------------
/// Main function
///-----------------------------------------------------------------------------
int main (int argc, char** argv)
{
	//----------------------------------------------------------
	// Read options
	//----------------------------------------------------------
	int flag_h = 0, flag_c = 0, flag_n = 0;
	const char* optDic = "hc:n:v:"; // Option dictionary
	int option;
	char* conf;
	long long int nevt = 1;
	unsigned short int verbose = 0;
	while ( (option = getopt(argc, argv, optDic)) != -1 ) // -1 means getopt() parses all options.
	{
		switch ( option )
		{
			case 'h' :
				flag_h = 1;
				break;
			case 'c' :
				flag_c = 1;
				conf = optarg;
				break;
			case 'n' :
				flag_n = 1;
				nevt = std::stoll(optarg);
				break;
			case 'v' :
				verbose = std::stoi(optarg);
				break;
			case '?' :
				flag_h = 1;
				break;
		}
	}


	//----------------------------------------------------------
	// Configuration manager 
	//----------------------------------------------------------
//	gCM = new ConMan();
//	if ( flag_c ) gCM -> Load(conf);
//	gCM -> Print();


	//----------------------------------------------------------
	// In case of '-h' option activated
	//----------------------------------------------------------
	if ( flag_h )
	{
		PrintHelp();
		return 0;
	};


	//----------------------------------------------------------
	// Generators
	//----------------------------------------------------------
	// Beam generator
	BeamGen* BG = new BeamGen();

	// Interaction point generation
	IPGen* IG = new IPGen();

	// Event generator
	EvtGen* EG = new EvtGen();


	//----------------------------------------------------------
	// Beam generation configuration
	//----------------------------------------------------------
	//--------------------------------------
	// Set config
	//--------------------------------------
	BG -> SetEnergy(200.);          // In MeV
	BG -> SetEnergyDispersion(.1);  // In MeV
	BG -> SetEmittance(1., 1.);     // ex, ey
	BG -> SetTwiss(1., 1., 1., 1.); // betax, alphax, betay, alphay
	BG -> SetMax(.25, 10.);         // x , y in mm

	//--------------------------------------
	// Generated beam will be stored in this variables
	//--------------------------------------
	double pxi = 0.;
	double pyi = 0.;
	double pzi = 0.;
	double ei  = 0.;
	double xi  = 0.;
	double yi  = 0.;


	//----------------------------------------------------------
	// Target configuration
	//----------------------------------------------------------
	//--------------------------------------
	// Set config
	//--------------------------------------
	IG -> SetTargetDiameter(.5); // In mm

	//--------------------------------------
	// Decided IP will be stored in this variables
	//--------------------------------------
	bool isScattered = false;
	double xip = 0.;
	double yip = 0.;
	double zip = 0.;


	//----------------------------------------------------------
	// Pluto configuration
	//----------------------------------------------------------


	//----------------------------------------------------------
	// Looping over target number of events
	//----------------------------------------------------------
	for ( long long int i = 0; i < nevt; i++ )
	{
		//--------------------------------------
		// Step 1: generate beam particle
		//--------------------------------------
		BG -> Generate(pxi, pyi, pzi, ei, xi, yi);

		if ( verbose )
		{
			std::cout << "[Info] Beam condition   : (pxi, pyi, pzi) = (" << pxi << ", " << pyi << ", " << pzi << ")" << std::endl;
			std::cout << "[Info]                    ei              = "  << ei  << " MeV"                            << std::endl;
			std::cout << "[Info]                    (xi, yi)        = (" << xi  << ", " << yi  << ") mm"             << std::endl;
		}

		//--------------------------------------
		// Step 2: decide location of interaction point inside the target
		//--------------------------------------
		IG -> SetBeamCondition(pxi, pyi, pzi, xi, yi);
		IG -> Generate(isScattered, xip, yip, zip);

		if ( !isScattered )
		{
			if ( verbose )
			{
				std::cout << "[Info] Event skipped without interaction" << std::endl;
				continue;
			}
		}

		if ( verbose )
		{
			std::cout << "[Info] Interaction point: (xip, yip, zip) = (" << xip << ", " << yip << ", " << zip << ") mm" << std::endl;
		}


		//--------------------------------------
		// Step 3: make proton scattered
		//--------------------------------------
		EG -> SetIP(xip, yip, zip);

	}



	//----------------------------------------------------------
	// Finalize
	//----------------------------------------------------------
	delete EG;
	delete IG;
	delete BG;

	return 0;
}



///-----------------------------------------------------------------------------
/// Print help message
///-----------------------------------------------------------------------------
void PrintHelp()
{
	std::cout << "usage: blpgen [-c configfile]" << std::endl;
	std::cout << std::endl;
	std::cout << "Examples:" << std::endl;
	std::cout << "  blpgen -c config.conf # Run with config.conf file." << std::endl;
	std::cout << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "  -c  Import config file" << std::endl;
	std::cout << "  -h  Show help message"  << std::endl;
	std::cout << std::endl;
	std::cout << "bye bye :)" << std::endl;
	std::cout << std::endl;
}
