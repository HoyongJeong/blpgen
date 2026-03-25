////////////////////////////////////////////////////////////////////////////////
///
///   IPGen.hh
///
/// - Author: Hoyong Jeong     (hoyong5419@korea.ac.kr)
///           Antonino Cannavo (acannavo@bnl.gov      )
///
////////////////////////////////////////////////////////////////////////////////



#ifndef IPGEN_h
#define IPGEN_h 1



///-----------------------------------------------------------------------------
/// Class declaration
///-----------------------------------------------------------------------------
class IPGen
{
	public:
	//----------------------------------------------------------
	// Constructor and destructor
	//----------------------------------------------------------
	IPGen();
	~IPGen();


	//----------------------------------------------------------
	// Public methods
	//----------------------------------------------------------
	void SetTargetDiameter(double d);
	void SetBeamCondition(double pxi, double pyi, double pzi, double xi, double yi);

	void Generate(bool& isScattered, double& xip, double& yip, double& zip);


	private:
	//----------------------------------------------------------
	// Target
	//----------------------------------------------------------
	// Diameter
	double m_d;

	
	//----------------------------------------------------------
	// Beam
	//----------------------------------------------------------
	// Momentum
	double m_pxi, m_pyi, m_pzi;

	// Cross point at z plane
	double m_xi, m_yi;
};



#endif
