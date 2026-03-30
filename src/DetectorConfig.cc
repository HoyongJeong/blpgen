////////////////////////////////////////////////////////////////////////////////
///
///   DetectorConfig.cc
///
/// Implementation of detector configuration management
///
/// This file defines the global detector configuration variables and
/// provides functions to modify them at runtime.
///
/// Key Features:
/// - Dynamic detector positioning (for systematic studies)
/// - Automatic recalculation of derived constants
/// - Configuration validation and reporting
///
/// Dependencies: DetectorConfig.h, TMath.h, iostream
///
/// ====================================================================
/// IMPORTANT NOTE FOR KINEMATICS MODULE
/// ====================================================================
/// The function ComputeCMAngleRange() in Kinematics.C uses a lookup
/// table that samples θ_CM from 0° to 50° (line 126). This range is
/// sufficient for the default detector position (θ = 16.2°) but must
/// be extended if detector is repositioned to larger angles.
///
/// Required modifications for different detector angles:
///   DETECTOR_THETA_CENTER < 35° → No change needed (default 50° OK)
///   DETECTOR_THETA_CENTER > 40° → Change to 90° range in Kinematics.C
///   DETECTOR_THETA_CENTER > 70° → Change to 180° range in Kinematics.C
///
/// To modify: Edit line 126 in Kinematics.C:
///   double theta_cm_deg = i * 50.0 / nPoints;  // Change 50.0
/// ====================================================================
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



///----------------------------------------------------------------------------
/// Headers
///----------------------------------------------------------------------------
#include "DetectorConfig.hh"

#include <iostream>



///----------------------------------------------------------------------------
/// Constructor
///----------------------------------------------------------------------------
DetectorConfig::DetectorConfig() :
	m_theta_center(16.2), // degrees
	m_theta_window(.005), // radians (±1.0°)
	m_phi_window(.005),   // radians (±1.0°)
	m_distance(2.18)      // meters
{
	CalculateDerivedConstants();
}



///----------------------------------------------------------------------------
/// Destructor
///----------------------------------------------------------------------------
DetectorConfig::~DetectorConfig()
{
}



///----------------------------------------------------------------------------
/// Initialize
///----------------------------------------------------------------------------
void DetectorConfig::CalculateDerivedConstants()
{
	m_theta_center_rad = m_theta_center * TMath::DegToRad();
	m_theta_min        = m_theta_center_rad - m_theta_window;
	m_theta_max        = m_theta_center_rad + m_theta_window;
	m_phi_min_0        = - m_phi_window;
	m_phi_max_0        =   m_phi_window;
	m_phi_min_180      = TMath::Pi() - m_phi_window;
	m_phi_max_180      = TMath::Pi() + m_phi_window;
}


///----------------------------------------------------------------------------
/// FUNCTION: SetDetectorConfig
///
/// Purpose: Update detector configuration and recalculate derived constants
///
/// Parameters:
///   theta_center  - Detector center angle [degrees]
///   theta_window  - Angular acceptance window [radians] (default: 0.070 rad ≈ 4°)
///   phi_window    - Azimuthal acceptance window [radians] (default: 0.017 rad ≈ 1°)
///   distance      - Polarimeter distance from target [meters] (default: 2.18 m)
///
/// Usage:
///   SetDetectorConfig(14.0);              // Set to 14.0°, use defaults for others
///   SetDetectorConfig(16.2, 0.05);        // Set to 16.2°, custom theta window
///   SetDetectorConfig(18.0, 0.070, 0.02); // Custom theta and phi windows
///
/// Notes:
///   - All derived constants are automatically recalculated
///   - Configuration is printed to screen for verification
///   - No validation checks (assumes user provides sensible values)
///----------------------------------------------------------------------------
void DetectorConfig::SetDetectorConfig(double theta_center,
                                       double theta_window,
                                       double phi_window,
                                       double distance)
{
	//----------------------------------------------------------
	// Update primary parameters
	//----------------------------------------------------------
	m_theta_center = theta_center;
	m_theta_window = theta_window;
	m_phi_window   = phi_window;
	m_distance     = distance;

	//----------------------------------------------------------
	// Recalculate derived constants
	//----------------------------------------------------------
	CalculateDerivedConstants();
    

	//----------------------------------------------------------
	// Print confirmation
	//----------------------------------------------------------
    std::cout << "\n=========================================="   << std::endl;
    std::cout << "   Detector Configuration Updated"              << std::endl;
    std::cout << "=========================================="     << std::endl;
    std::cout << "  Theta center:  "  << m_theta_center << " deg" << std::endl;
    std::cout << "  Theta window:  ±" << m_theta_window << " rad";
    std::cout << " (±" << m_theta_window * TMath::RadToDeg() << " deg)" << std::endl;
    std::cout << "  Phi window:    ±" << m_phi_window   << " rad";
    std::cout << " (±" << m_phi_window   * TMath::RadToDeg() << " deg)" << std::endl;
    std::cout << "  Distance:      "  << m_distance     << " m"   << std::endl;
    std::cout << "=========================================="     << std::endl;
    std::cout << "  Theta range:   [" << m_theta_min * TMath::RadToDeg();
    std::cout << ", " << m_theta_max * TMath::RadToDeg() << "] deg" << std::endl;
    std::cout << "=========================================="     << std::endl << std::endl;
}



///----------------------------------------------------------------------------
/// FUNCTION: PrintDetectorConfig
///
/// Purpose: Display current detector configuration (detailed output)
///
/// Usage:
///   PrintDetectorConfig();
///
/// Output: Shows all primary and derived configuration parameters
///----------------------------------------------------------------------------
void DetectorConfig::PrintDetectorConfig()
{
	
	std::cout << "\n=========================================="  << std::endl;
	std::cout << "   Current Detector Configuration"             << std::endl;
	std::cout << "=========================================="    << std::endl;
	std::cout << "\nPrimary Parameters:" << std::endl;
	std::cout << "  theta_center = " << m_theta_center << " deg" << std::endl;
	std::cout << "  theta_window = " << m_theta_window << " rad";
	std::cout << " (±" << m_theta_window * TMath::RadToDeg() << " deg)" << std::endl;
	std::cout << "  phi_window   = " << m_phi_window   << " rad";
	std::cout << " (±" << m_phi_window   * TMath::RadToDeg() << " deg)" << std::endl;
	std::cout << "  distance     = " << m_distance     << " m"   << std::endl;
	
	std::cout << "\nDerived Constants:" << std::endl;
	std::cout << "  theta_center_rad = " << m_theta_center_rad << " rad" << std::endl;
	std::cout << "  theta_min        = " << m_theta_min        << " rad";
	std::cout << " (" << m_theta_min * TMath::RadToDeg() << " deg)" << std::endl;
	std::cout << "  theta_max        = " << m_theta_max        << " rad";
	std::cout << " (" << m_theta_max * TMath::RadToDeg() << " deg)" << std::endl;
	std::cout << "  phi_min_0        = " << m_phi_min_0        << " rad" << std::endl;
	std::cout << "  phi_max_0        = " << m_phi_max_0        << " rad" << std::endl;
	std::cout << "  phi_min_18 0     = " << m_phi_min_180      << " rad" << std::endl;
	std::cout << "  phi_max_180      = " << m_phi_max_180      << " rad" << std::endl;

	std::cout << "\nAcceptance Summary:" << std::endl;
	std::cout << "  Theta range: [" << m_theta_min   * TMath::RadToDeg();
	std::cout << ", " << m_theta_max   * TMath::RadToDeg() << "] deg" << std::endl;
	std::cout << "  Phi @ 0°:    [" << m_phi_min_0   * TMath::RadToDeg();
	std::cout << ", " << m_phi_max_0   * TMath::RadToDeg() << "] deg" << std::endl;
	std::cout << "  Phi @ 180°:  [" << m_phi_min_180 * TMath::RadToDeg();
	std::cout << ", " << m_phi_max_180 * TMath::RadToDeg() << "] deg" << std::endl;
	std::cout << "==========================================" << std::endl << std::endl;
}



///----------------------------------------------------------------------------
/// FUNCTION: ResetDetectorConfig
///
/// Purpose: Reset detector configuration to default values
///
/// Usage:
///   ResetDetectorConfig();
///
/// Notes: Default configuration is 16.2° center with standard windows
///----------------------------------------------------------------------------
void DetectorConfig::ResetDetectorConfig()
{
	std::cout << "\nResetting detector configuration to defaults..." << std::endl;
	SetDetectorConfig(16.2, 0.005, 0.005, 2.18);
}
