////////////////////////////////////////////////////////////////////////////////
///
///   ElastricScattering.cc
///
/// Implementation of elastic scattering cross section conversions/
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



///----------------------------------------------------------------------------
/// Headers
///----------------------------------------------------------------------------
#include "ElasticScattering.hh"

#include <iostream>



///----------------------------------------------------------------------------
/// Constructor
///----------------------------------------------------------------------------
ElasticScattering::ElasticScattering()
{
}



///----------------------------------------------------------------------------
/// Destructor
///----------------------------------------------------------------------------
ElasticScattering::~ElasticScattering()
{
}



///----------------------------------------------------------------------------
/// Convert log cross section to linear scale
/// The input files provide log10(σ), so we convert: σ = 10^(log10(σ))
///----------------------------------------------------------------------------
double ElasticScattering::fXS150Op(double x)
{
	return TMath::Power(10, XSLog150MeV(x));
}

double ElasticScattering::fXS160Op(double x)
{
	return TMath::Power(10, XSLog160MeV(x));
}

double ElasticScattering::fXS170Op(double x)
{
	return TMath::Power(10, XSLog170MeV(x));
}

double ElasticScattering::fXS180Op(double x)
{
	return TMath::Power(10, XSLog180MeV(x));
}

double ElasticScattering::fXS190Op(double x)
{
	return TMath::Power(10, XSLog190MeV(x));
}

double ElasticScattering::fXS200Op(double x)
{
	return TMath::Power(10, XSLog200MeV(x));
}

double ElasticScattering::fXS210Op(double x)
{
	return TMath::Power(10, XSLog210MeV(x));
}

double ElasticScattering::fXS220Op(double x)
{
	return TMath::Power(10, XSLog220MeV(x));
}

double ElasticScattering::fXS230Op(double x)
{
	return TMath::Power(10, XSLog230MeV(x));
}

double ElasticScattering::fXS240Op(double x)
{
	return TMath::Power(10, XSLog240MeV(x));
}



///----------------------------------------------------------------------------
/// Helper function: Get cross section formula string for given energy
/// Returns the appropriate function name based on beam kinetic energy
///----------------------------------------------------------------------------
const char* ElasticScattering::GetXSFormula(double ekin)
{
	if      ( ekin == 150. ) return "ElasticScattering::fXS150Op(x)";
	else if ( ekin == 160. ) return "ElasticScattering::fXS160Op(x)";
	else if ( ekin == 170. ) return "ElasticScattering::fXS170Op(x)";
	else if ( ekin == 180. ) return "ElasticScattering::fXS180Op(x)";
	else if ( ekin == 190. ) return "ElasticScattering::fXS190Op(x)";
	else if ( ekin == 200. ) return "ElasticScattering::fXS200Op(x)";
	else if ( ekin == 210. ) return "ElasticScattering::fXS210Op(x)";
	else if ( ekin == 220. ) return "ElasticScattering::fXS220Op(x)";
	else if ( ekin == 230. ) return "ElasticScattering::fXS230Op(x)";
	else if ( ekin == 240. ) return "ElasticScattering::fXS240Op(x)";
	else
	{
		std::cout << "Warning: Energy " << ekin << " not in optical model data." << std::endl;
		std::cout << "Using cos(theta) distribution instead." << std::endl;
		return "cos(x)";
	}
}



///----------------------------------------------------------------------------
/// GetElasticAnalyzingPower: Calculate A_N for elastic scattering
///----------------------------------------------------------------------------
double ElasticScattering::GetElasticAnalyzingPower(double ekin, double theta_lab_deg)
{
	double dT = ekin - T0_AN;
	double dTheta = theta_lab_deg - theta0_AN;
	double AN = 1. - alpha_AN*dT*dT - beta_AN*dT*dTheta - gamma_AN*dTheta*dTheta;


	//----------------------------------------------------------
	// Clamp to physical range [-1, 1]
	//----------------------------------------------------------
	if ( AN >   1. ) AN =   1.;
	if ( AN < - 1. ) AN = - 1.;


	return AN;
}
