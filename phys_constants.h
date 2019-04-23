#ifndef PHYS_CONSTANTS_H_INCLUDED
#define PHYS_CONSTANTS_H_INCLUDED
const double length_const = 1e-10;
const double mass_const = 1e-23;
const double time_const = 1e-15;
const double volume_const = length_const*length_const*length_const;
const double force_const = mass_const*length_const/(time_const*time_const);//1e-3
const double pressure_const = force_const/(length_const*length_const);//1e+17
const double stress_const = force_const/(length_const*length_const);//1e+17
const double velocity_const = length_const/time_const;//1e-5
const double energy_const = force_const*length_const;//1e-13

const double N_a_PhysConst = 6.022045e23;
const double e_q_PhysConst = 1.6e-19;   // ��
const double k_PhysConst = 1.380662e-23/energy_const;// ��/�
const double _1d_k_PhysConst = 1.0/k_PhysConst;

const double CoordSphereDistance2DTriangle_Const[]= {0, 1, 1.7320508075688772, 2, 2.6457513110645907, 3, 3.4641016151377544, 3.6055512754639891, 4, 4.358898943540674, 4.5825756949558398, 5};

struct VoigtNotationType
{
    double k;
    uint_fast8_t i,j;

};

//const uint_fast8_t VoigtNotation3D[9]={0, 5, 4, 5, 1, 3, 4, 3, 2};
const VoigtNotationType VoigtNotation2D[4][4]={ {{1.0, 0,0}, {0.5 ,0,2}, {0.5, 0,2}, {1.0, 0,1}},
                                                {{0.5, 2,0}, {0.25,2,2}, {0.25,2,2}, {0.5, 2,1}},
                                                {{0.5, 2,0}, {0.25,2,2}, {0.25,2,2}, {0.5, 2,1}},
                                                {{1.0, 1,0}, {0.5, 1,2}, {0.5, 1,2}, {1.0, 1,1}}};

const VoigtNotationType VoigtNotation3D[9][9]={ {{1.0, 0,0}, {0.5 ,0,5}, {0.5, 0,4}, {0.5, 0,5},  {1.0, 0,1}, {0.5 ,0,3}, {0.5 ,0,4}, {0.5 ,0,3},  {1.0, 0,2}},
                                                {{0.5, 5,0}, {0.25,5,5}, {0.25,5,4}, {0.25, 5,5}, {0.5, 5,1}, {0.25,5,3}, {0.25,5,4}, {0.25, 5,3}, {0.5, 5,2}},
                                                {{0.5, 4,0}, {0.25,4,5}, {0.25,4,4}, {0.25, 4,5}, {0.5, 4,1}, {0.25,4,3}, {0.25,4,4}, {0.25, 4,3}, {0.5, 4,2}},
                                                {{0.5, 5,0}, {0.25,5,5}, {0.25,5,4}, {0.25, 5,5}, {0.5, 5,1}, {0.25,5,3}, {0.25,5,4}, {0.25, 5,3}, {0.5, 5,2}},
                                                {{1.0, 1,0}, {0.5 ,1,5}, {0.5, 1,4}, {0.5,  1,5}, {1.0, 1,1}, {0.5 ,1,3}, {0.5 ,1,4}, {0.5,  1,3}, {1.0, 1,2}},
                                                {{0.5, 3,0}, {0.25,3,5}, {0.25,3,4}, {0.25, 3,5}, {0.5, 3,1}, {0.25,3,3}, {0.25,3,4}, {0.25, 3,3}, {0.5, 3,2}},
                                                {{0.5, 4,0}, {0.25,4,5}, {0.25,4,4}, {0.25, 4,5}, {0.5, 4,1}, {0.25,4,3}, {0.25,4,4}, {0.25, 4,3}, {0.5, 4,2}},
                                                {{0.5, 3,0}, {0.25,3,5}, {0.25,3,4}, {0.25, 3,5}, {0.5, 3,1}, {0.25,3,3}, {0.25,3,4}, {0.25, 3,3}, {0.5, 3,2}},
                                                {{1.0, 2,0}, {0.5 ,2,5}, {0.5, 2,4}, {0.5,  2,5}, {1.0, 2,1}, {0.5 ,2,3}, {0.5 ,2,4}, {0.5,  2,3}, {1.0, 2,2}}
                                                };


/*const double N_a = 6.022045e23;
const double m_0 = 472.0/N_a*mass_const;

const double ao_1 = 0.39e-9*length_const;
const double C_1 = sqrt(3.0)*0.5*3e9*ao_1*force_const/(length_const*length_const);				// particle radius
const double D_1 = C_1*ao_1*ao_1/72.0;
						// interparticle force magnitude
const double ao_2 = 0.142e-9*length_const;							// particle radius
const double C_2 = 2.0*730*force_const/length_const;				// particle radius
const double D_2 = C_2*ao_2*ao_2/72.0;

const double To_1 = 2 * PI * sqrt(m_0 / C_1);	// Oscillation period
const double To_2 = 2 * PI * sqrt(m_0 / C_2);	// Oscillation period
const double vo_1 = sqrt(ao_1 * ao_1 * C_1/m_0);			// Long-wave velocity
const double vo_2 = sqrt(ao_2 * ao_2 * C_2/m_0);			// Long-wave velocity
const double Bo_1 = 2 * sqrt(m_0*C_1);		// Critical friction
const double Bo_2 = 2 * sqrt(m_0*C_2);		// Critical friction
const double a0_1	= 1		* ao_1;		// initial distance between particles
const double a0_2	= 1		* ao_2;		// initial distance between particles
const LONG n1	= 15;			// number of particles in a single x-row
const LONG n2	= 15;				// number of particles in a single y-row 1

const double B_1		= 0.05	* Bo_1;	// viscouse friction
const double B_2		= 0.05	* Bo_2;	// viscouse friction
const double dt	= 0.2	* min(To_1,To_2) / 10;	// time step for computation
const double v0	= 0.01	* min(vo_1,vo_2);		// maximum value for initial random velocity
const double a_cut_1 = 1.3	* ao_1;		// cut distance for interatomic forces
const double a_cut_2 = 4.2;//2.958 * ao_2;		// cut distance for interatomic forces
const double P = 0.426e12/time_const; //����������� ����������� ����� � ������� �������*/


#endif // PHYS_CONSTANTS_H_INCLUDED
