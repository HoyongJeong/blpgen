////////////////////////////////////////////////////////////////////////////////
///
///   DetectorConfig.hh
///
/// Detector geometry and acceptance parameters
/// 
/// This header declares global variables for detector configuration.
/// The actual values are defined in DetectorConfig.C and can be 
/// modified at runtime using SetDetectorConfig().
///
/// Usage:
///   SetDetectorConfig(16.2);  // Set detector to 16.4 degrees
///   PrintDetectorConfig();     // Display current configuration
///
/// Dependencies: TMath.h
/// Used by: All modules requiring detector acceptance cuts
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



#ifndef DETECTORCONFIG_H
#define DETECTORCONFIG_H



///-----------------------------------------------------------------------------
/// Headers
///-----------------------------------------------------------------------------
#include "TMath.h"



///-----------------------------------------------------------------------------
/// Class declaration
///-----------------------------------------------------------------------------
class DetectorConfig
{
	public:
	//----------------------------------------------------------
	// Constructor and destructor
	//----------------------------------------------------------
	DetectorConfig();
	~DetectorConfig();


	//----------------------------------------------------------
	// Public methods
	//----------------------------------------------------------
	// Initialize
	void CalculateDerivedConstants();

	// Set detector configuration and recalculate derived constants
	void SetDetectorConfig(Double_t theta_center,         // degrees
	                       Double_t theta_window = 0.005, // radians (default ±4°)
	                       Double_t phi_window   = 0.005, // radians (default ±1°)
	                       Double_t distance     = 2.18); // meters
	
	// Display current detector configuration
	void PrintDetectorConfig();

	// Reset to default configuration (theta = 16.2°)
	void ResetDetectorConfig();


	private:
	//----------------------------------------------------------
	// PRIMARY DETECTOR PARAMETERS
	// These can be modified using SetDetectorConfig()
	//----------------------------------------------------------
	double m_theta_center;     // Center angle [degrees]
	double m_theta_window;     // Angular acceptance [radians]
	double m_phi_window;       // Azimuthal acceptance [radians]
	double m_distance;         // Distance from target [meters]


	//----------------------------------------------------------
	// DERIVED CONSTANTS
	// These are automatically recalculated when SetDetectorConfig() is called
	//----------------------------------------------------------
	double m_theta_center_rad; // Center angle [radians]
	double m_theta_min;        // Min acceptance angle [radians]
	double m_theta_max;        // Max acceptance angle [radians]
	double m_phi_min_0;        // Min phi at 0° [radians]
	double m_phi_max_0;        // Max phi at 0° [radians]
	double m_phi_min_180;      // Min phi at 180° [radians]
	double m_phi_max_180;      // Max phi at 180° [radians]
};



#endif
