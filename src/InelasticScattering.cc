////////////////////////////////////////////////////////////////////////////////
///
///   ElastricScattering.cc
///
///   Implementation of inelastic scattering data loading and access
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



///----------------------------------------------------------------------------
/// Headers
///----------------------------------------------------------------------------
#include "PhysicalConstants.hh"
#include "InelasticScattering.hh"

#include "TSystem.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>



///----------------------------------------------------------------------------
/// Constructor
///----------------------------------------------------------------------------
InelasticScattering::InelasticScattering() : mg_inelastic_xs(nullptr), mg_inelastic_AN(nullptr)
{
}



///----------------------------------------------------------------------------
/// Destructor
///----------------------------------------------------------------------------
InelasticScattering::~InelasticScattering()
{
	if ( mg_inelastic_xs ) delete mg_inelastic_xs;
	if ( mg_inelastic_AN ) delete mg_inelastic_AN;
}



///----------------------------------------------------------------------------
/// LoadInelasticData: Load digitized inelastic scattering data
///----------------------------------------------------------------------------
void InelasticScattering::LoadInelasticData()
{
	std::cout << "\n========================================" << std::endl;
	std::cout << "Loading Inelastic Data"                     << std::endl;
	std::cout << "========================================"   << std::endl;


	//----------------------------------------------------------
	// Load Cross Section
	//----------------------------------------------------------
	const char* xs_file = "inelastic_crosssection_443_200MeV.csv";
	std::ifstream file_xs(xs_file);

	if ( ! file_xs . is_open() )
	{
		std::cerr << "ERROR: Cannot open " << xs_file                   << std::endl;
		std::cerr << "Current directory: " << gSystem -> pwd()          << std::endl;
		std::cerr << "Make sure the CSV file is in the same directory!" << std::endl;
		exit(1);
	}

	
	std::vector<double> theta_xs, sigma_xs;
	std::string line;

	
	while ( getline(file_xs, line) )
	{
		if ( line . empty() || line[0] == '#' ) continue;

		double theta, sigma;
		if ( sscanf(line . c_str(), "%lf,%lf" , &theta, &sigma) == 2 ||
		     sscanf(line . c_str(), "%lf, %lf", &theta, &sigma) == 2 )
		{
			theta_xs . push_back(theta);
			sigma_xs . push_back(sigma);
		}
	}
	file_xs . close();

	
	mg_inelastic_xs = new TGraph(theta_xs . size(), &theta_xs[0], &sigma_xs[0]);
	mg_inelastic_xs -> SetName("g_inelastic_xs");
	std::cout << "✓ Cross section: " << theta_xs . size() << " points" << std::endl;


	//----------------------------------------------------------
	// Load Analyzing Power
	//----------------------------------------------------------
	const char* an_file = "inelastic_analyzingpower_443_200MeV.csv";
	std::ifstream file_an(an_file);

	if ( ! file_an . is_open() )
	{
		std::cerr << "ERROR: Cannot open " << an_file << std::endl;
		std::cerr << "Current directory: " << gSystem -> pwd() << std::endl;
		std::cerr << "Make sure the CSV file is in the same directory!" << std::endl;
		exit(1);
	}

	std::vector<double> theta_an, an_values;

	while ( getline(file_an, line) )
	{
		if ( line . empty() || line[0] == '#' ) continue;
		double theta, an;

		if ( sscanf(line . c_str(), "%lf,%lf" , &theta, &an) == 2 ||
		     sscanf(line . c_str(), "%lf, %lf", &theta, &an) == 2 )
		{
			theta_an  . push_back(theta);
			an_values . push_back(an);
		}
	}
	file_an . close();

	mg_inelastic_AN = new TGraph(theta_an . size(), &theta_an[0], &an_values[0]);
	mg_inelastic_AN -> SetName("g_inelastic_AN");

	std::cout << "✓ Analyzing power: " << theta_an . size() << " points" << std::endl;
	std::cout << "========================================\n" << std::endl;
}



///----------------------------------------------------------------------------
/// GetInelasticCrossSection: Interpolate cross section at given angle
///----------------------------------------------------------------------------
double InelasticScattering::GetInelasticCrossSection(double theta_cm_deg)
{
	if ( ! mg_inelastic_xs )
	{
		std::cerr << "ERROR: Call LoadInelasticData() first!" << std::endl;
		exit(1);
	}


	return mg_inelastic_xs -> Eval(theta_cm_deg);
}



///----------------------------------------------------------------------------
/// GetInelasticAnalyzingPower: Interpolate A_N at given angle
///----------------------------------------------------------------------------
double InelasticScattering::GetInelasticAnalyzingPower(double theta_cm_deg)
{
	if ( ! mg_inelastic_AN )
	{
		std::cerr << "ERROR: Call LoadInelasticData() first!" << std::endl;
		exit(1);
	}


	return mg_inelastic_AN -> Eval(theta_cm_deg);
}



///----------------------------------------------------------------------------
/// CalculateCMMomentumInelastic: Calculate CM momentum for inelastic
///----------------------------------------------------------------------------
double InelasticScattering::CalculateCMMomentumInelastic(double ekin, double E_ex)
{
	//----------------------------------------------------------
	// Convert to GeV
	//----------------------------------------------------------
	double ekinGeV = ekin / 1000.;


	//----------------------------------------------------------
	// Initial proton momentum
	//----------------------------------------------------------
	double mom = TMath::Sqrt(2.*ekinGeV*mp + ekinGeV*ekinGeV);


	//----------------------------------------------------------
	// Initial 4-momentum: beam proton + target carbon at rest
	//----------------------------------------------------------
	TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);


	//----------------------------------------------------------
	// Mandelstam s (total energy squared in CM)
	//----------------------------------------------------------
    double s = iState . Mag2();


	//----------------------------------------------------------
	// For inelastic: final state has carbon with excitation energy
	// Effective mass of excited carbon: m_C* = m_C + E_ex
	//----------------------------------------------------------
	double mC_star = mC + E_ex;


	//----------------------------------------------------------
	// CM momentum formula with excited final state
	//----------------------------------------------------------
	double pcm_inel = TMath::Sqrt((s - (mp + mC_star)*(mp + mC_star)) *
                                     (s - (mp - mC_star)*(mp - mC_star)) / (4.*s));


	return pcm_inel;
}



///----------------------------------------------------------------------------
/// CreateInelasticSamplingHistogram: Create histogram for importance sampling
///
/// Warning! This method has potenrial memory leak risk.
///----------------------------------------------------------------------------
TH1D* InelasticScattering::CreateInelasticSamplingHistogram(double theta_cm_min, double theta_cm_max, int nbins)
{
	if ( ! mg_inelastic_xs )
	{
		std::cerr << "ERROR: Call LoadInelasticData() first!" << std::endl;
		return nullptr;
	}

	
	TH1D* h = new TH1D("h_inelastic_sampling",
                       "Inelastic cross section sampling",
                       nbins, theta_cm_min, theta_cm_max);
	
	// Fill histogram with interpolated cross section values
	for ( int ibin = 1; ibin <= nbins; ibin++ )
	{
		double theta_cm = h -> GetBinCenter(ibin);
		double xs = GetInelasticCrossSection(theta_cm);
		h -> SetBinContent(ibin, xs);
	}


	return h;
}
