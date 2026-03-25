////////////////////////////////////////////////////////////////////////////////
///
///   BeamGen.hh
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



#ifndef BEAMGEN_h
#define BEAMGEN_h 1



///-----------------------------------------------------------------------------
/// Class declaration
///-----------------------------------------------------------------------------
class BeamGen
{
	public:
	//----------------------------------------------------------
	// Constructor and destructor
	//----------------------------------------------------------
	BeamGen();
	~BeamGen();


	//----------------------------------------------------------
	// Public methods
	//----------------------------------------------------------
	void SetEnergy(double E);
	void SetEnergyDispersion(double dE);
	void SetEmittance(double ex, double ey);
	void SetTwiss(double bx, double ax, double by, double ay);
	void SetMax(double xmax, double ymax);

	void Generate(double& pxi, double& pyi, double& pzi, double& ei, double& xi, double& yi);


	private:
	//----------------------------------------------------------
	// Proton energy
	//----------------------------------------------------------
	// Central
	double m_E;

	// Dispersion
	double m_dE;

	
	//----------------------------------------------------------
	// Beam emittance
	//----------------------------------------------------------
	// x-axis
	double m_ex;

	// y-axis
	double m_ey;


	//----------------------------------------------------------
	// Twiss parameters
	//----------------------------------------------------------
	// x-axis beta
	double m_bx;

	// x-axis alpha
	double m_ax;

	// y-axis beta
	double m_by;

	//y-axis alpha
	double m_ay;


	//----------------------------------------------------------
	// Position limit
	//----------------------------------------------------------
	// x-axis
	double m_xmax;

	// y-axis
	double m_ymax;
};



#endif
