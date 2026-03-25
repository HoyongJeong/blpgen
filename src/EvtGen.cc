////////////////////////////////////////////////////////////////////////////////
///
///   EvtGen.cc
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



///----------------------------------------------------------------------------
/// Headers
///----------------------------------------------------------------------------
#include "EvtGen.hh"



///----------------------------------------------------------------------------
/// Constructor
///----------------------------------------------------------------------------
EvtGen::EvtGen() : m_xip(0.), m_yip(0.), m_zip(0.)
{
}



///----------------------------------------------------------------------------
/// Destructor
///----------------------------------------------------------------------------
EvtGen::~EvtGen()
{
}



///----------------------------------------------------------------------------
/// Set target diameter
///----------------------------------------------------------------------------
void EvtGen::SetIP(double xip, double yip, double zip)
{
	m_xip = xip;
	m_yip = yip;
	m_zip = zip;

	return;
}



///----------------------------------------------------------------------------
/// Generate
///----------------------------------------------------------------------------
void EvtGen::Generate()
{
	return;
}



///----------------------------------------------------------------------------
/// GetAnalyzingPower: Calculate A_N for given energy and lab angle
///
/// Formula: A_N(T,θ) = 1 - α(T-T₀)² - β(T-T₀)(θ-θ₀) - γ(θ-θ₀)²
///
/// Input:
///   - ekin: beam kinetic energy [MeV]
///   - theta_lab_deg: scattering angle in lab frame [degrees]
///
/// Output:
///   - analyzing power A_N (dimensionless, typically 0 to 1)
///----------------------------------------------------------------------------
double EvtGen::GetAnalyzingPower(double ekin, double theta_lab_deg)
{
	double dT = ekin - T0_AN;
	double dTheta = theta_lab_deg - theta0_AN;
	double AN = 1.0 - alpha_AN*dT*dT - beta_AN*dT*dTheta - gamma_AN*dTheta*dTheta;


	//----------------------------------------------------------
	// Clamp to physical range [-1, 1]
	//----------------------------------------------------------
	if ( AN >  1.0 ) AN =  1.0;
	if ( AN < -1.0 ) AN = -1.0;


	return AN;
}
