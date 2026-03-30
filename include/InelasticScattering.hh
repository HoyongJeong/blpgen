////////////////////////////////////////////////////////////////////////////////
///
///   InelasticScattering.hh
///
///   Inelastic scattering data management for p+C -> p+C* (4.43 MeV)
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



#ifndef INELASTICSCATTERING_H
#define INELASTICSCATTERING_H



///-----------------------------------------------------------------------------
/// Headers
///-----------------------------------------------------------------------------
#include "TGraph.h"
#include "TH1D.h"



///-----------------------------------------------------------------------------
/// Class declaration
///-----------------------------------------------------------------------------
class InelasticScattering
{
	public:
	//----------------------------------------------------------
	// Constructor and destructor
	//----------------------------------------------------------
	InelasticScattering();
	~InelasticScattering();


	//----------------------------------------------------------
	// LoadInelasticData: Load digitized inelastic scattering data
	//----------------------------------------------------------
	// Reads CSV files and creates interpolation graphs
	// Call this ONCE at the start of your simulation
	//
	// Required files (same directory as where ROOT was started):
	//   - inelastic_crosssection_443_200MeV.csv
	//   - inelastic_analyzingpower_443_200MeV.csv
	//----------------------------------------------------------
	void LoadInelasticData();


	//----------------------------------------------------------
	// GetInelasticCrossSection: Interpolate cross section at given angle
	//----------------------------------------------------------
	// Input: theta_cm_deg - CM angle in degrees
	// Output: cross section in mb/sr (NOTE: Not log scale!)
	//----------------------------------------------------------
	double GetInelasticCrossSection(double theta_cm_deg);


	//----------------------------------------------------------
	// GetInelasticAnalyzingPower: Interpolate A_N at given angle
	//----------------------------------------------------------
	// Input: theta_cm_deg - CM angle in degrees
	// Output: A_N (dimensionless)
	//----------------------------------------------------------
	double GetInelasticAnalyzingPower(double theta_cm_deg);


	//----------------------------------------------------------
	// Check if data is loaded
	//----------------------------------------------------------
	inline bool IsInelasticDataLoaded()
	{
		return ( mg_inelastic_xs != nullptr && mg_inelastic_AN != nullptr );
	}


	//----------------------------------------------------------
	// CalculateCMMomentumInelastic: Calculate CM momentum for inelastic
	//----------------------------------------------------------
	// For inelastic scattering, less energy is available in the final state
	// because E_excitation goes into exciting the nucleus
	//
	// Input:
	//   - ekin: beam kinetic energy [MeV]
	//   - E_ex: excitation energy [GeV] (use 0.00443 for 4.43 MeV state)
	//
	// Output:
	//   - pcm_inelastic: CM momentum [MeV/c]
	//----------------------------------------------------------
	double CalculateCMMomentumInelastic(double ekin, double E_ex = 4.43);


	//----------------------------------------------------------
	// CreateInelasticSamplingHistogram: Create histogram for importance sampling
	//----------------------------------------------------------
	// Creates histogram with bin contents = interpolated cross sections
	// Used by TH1::GetRandom() for importance sampling (like elastic case)
	//
	// Inputs:
	//   - theta_cm_min, theta_cm_max: angle range in CM frame [degrees]
	//   - nbins: number of bins (default 10000)
	//
	// Output:
	//   - TH1D* ready for GetRandom() sampling
	//----------------------------------------------------------
	TH1D* CreateInelasticSamplingHistogram(double theta_cm_min, double theta_cm_max, int nbins = 10000);


	private:
	TGraph* mg_inelastic_xs; // Cross section   vs theta_CM (degrees)
	TGraph* mg_inelastic_AN; // Analyzing power vs theta_CM (degrees)
};



#endif
