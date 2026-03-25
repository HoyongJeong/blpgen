////////////////////////////////////////////////////////////////////////////////
///
///   EvtGen.hh
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



#ifndef EVTGEN_h
#define EVTGEN_h 1



///-----------------------------------------------------------------------------
/// Headers
///-----------------------------------------------------------------------------
#include "TMath.h"


///-----------------------------------------------------------------------------
/// Class declaration
///-----------------------------------------------------------------------------
class EvtGen
{
	public:
	//----------------------------------------------------------
	// Constructor and destructor
	//----------------------------------------------------------
	EvtGen();
	~EvtGen();


	//----------------------------------------------------------
	// Public methods
	//----------------------------------------------------------
	void SetIP(double xip, double yip, double zip);

	void Generate();


	private:
	//----------------------------------------------------------
	// Interation point
	//----------------------------------------------------------
	double m_xip;
	double m_yip;
	double m_zip;


	//----------------------------------------------------------
	// Physical constants
	//----------------------------------------------------------
	const double mp =   938.2720813; // Proton mass in MeV/c2
	const double mC = 11174.862;     // Carbon-12 mass in MeV/c2


	//----------------------------------------------------------
	// Analyzing Power Parameterization
	// Based on: https://pos.sissa.it/324/024/pdf
	// Wissink et al., Phys. Rev. C 45, R504 (1992)
	//----------------------------------------------------------
	const double T0_AN     = 189.0;   // MeV - energy at maximum A_N
	const double theta0_AN = 17.3;    // degrees - angle at maximum A_N (lab frame)
	const double alpha_AN  = 1.21e-4; // MeV^-2
	const double beta_AN   = 1.61e-3; // MeV^-1 deg^-1
	const double gamma_AN  = 1.00e-2; // deg^-2


	//----------------------------------------------------------
	// Detector acceptance configuration
	//----------------------------------------------------------
	const double DETECTOR_THETA_CENTER = 16.0; // degrees
	const double DETECTOR_THETA_WINDOW = 0.1;  // radians (±5 mrad)
	const double DETECTOR_PHI_WINDOW   = 3.;   // radians (±5 mrad)
	const double POLARIMETER_DISTANCE  = 2.18; // meters

	// Derived constants
	const double DETECTOR_THETA_CENTER_RAD = DETECTOR_THETA_CENTER * TMath::DegToRad();
	const double DETECTOR_THETA_MIN        = DETECTOR_THETA_CENTER_RAD - DETECTOR_THETA_WINDOW;
	const double DETECTOR_THETA_MAX        = DETECTOR_THETA_CENTER_RAD + DETECTOR_THETA_WINDOW;
	const double DETECTOR_PHI_MIN_0        = -DETECTOR_PHI_WINDOW;
	const double DETECTOR_PHI_MAX_0        = DETECTOR_PHI_WINDOW;
	const double DETECTOR_PHI_MIN_180      = TMath::Pi() - DETECTOR_PHI_WINDOW;
	const double DETECTOR_PHI_MAX_180      = -TMath::Pi() + DETECTOR_PHI_WINDOW;


	//----------------------------------------------------------
	// Beam Polarization
	//----------------------------------------------------------
	const double BEAM_POLARIZATION = 0.80;  // 80% polarization


	//----------------------------------------------------------
	// GetAnalyzingPower: Calculate A_N for given energy and lab angle
	//----------------------------------------------------------
	double GetAnalyzingPower(double ekin, double theta_lab_deg);
};



#endif
