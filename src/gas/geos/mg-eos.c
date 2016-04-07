/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*
*				mg-eos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <geos/mg.h>



/*
*	Use a 110 column window or larger to see the table layout into columns 
*
*	When using materials from this table the units are specified as 
*
*	Density         - g/cc
*	Pressure        - Megabars
*	Energy          - Megabars cc = 0.1 Megajoules
*	Specific Energy - Megabars cc / gram = 0.1 Megajoules/gram
*	Velocity        - cm/microsecond = 10 km/sec
*	Temperature     - Kelvin
*
*
*	Reference: Equation of State and Strength Properties of Selected
*	           Materials
*		   Daniel J. Steinberg
*                  Lawrence Livermore National Laboratory Report
*                  UCRL-MA-106439 Change 1
*
*	Note in the material database Steinberg did not list specific heats
*	at constant volume.  The values included in the table are the value
*	quoted for C_p on this melt curve.
*
*	NOTES -
*	MATERIALS: lithium, mercury, graphite, Kel-F, polyethylene (CH2), and
*	           teflon
*	The value of CV for these materials were not in Steinberg's table the
*	value used were obtained using the NIST webbook:
*	http://webbook.nist.gov.  The nominal value of Cp was used to fit
*	CV = Cp/gamma0. A rough HACK beware.
*
*	MATERIALS: thulium, lexan, micarta, nylon, polypentene, and silastic
*	The value of CV for these materials were made up.  No data was
*	given in Steinberg's table or the NIST webbook.
*
*	MATERIALS: Tungsten (98.52%) in Plastic Binder
*	The value of CV for this material was not included in Steinberg's
*	tables.  The value for the pure metal was used instead.
*/

LOCAL MG_params MGMaterial[] = {
/*    name          ,
*    Rho0,  P0,  E0,  T0,     C0,    S1,      S2,     S3, gamma0,      b,       CV, Rho_max, V_min, P_inf, E_inf, TrefSpline
*/
"Aluminum (1100-O)",
    2.707, 0.0, 0.0, 0.0, 0.5386, 1.339,  0.0000, 0.0000, 1.9700, 0.4800, 0.884e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Aluminum (2024-T4)",
    2.785, 0.0, 0.0, 0.0, 0.5328, 1.338,  0.0000, 0.0000, 2.0000, 0.4800, 0.863e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Aluminum (6061-T6)",
    2.703, 0.0, 0.0, 0.0, 0.5240, 1.400,  0.0000, 0.0000, 1.9700, 0.4800, 0.885e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Aluminum (7075-T6)",
    2.804, 0.0, 0.0, 0.0, 0.5200, 1.360,  0.0000, 0.0000, 2.2000, 0.4800, 0.848e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Beryllium (S200)",
    1.850, 0.0, 0.0, 0.0, 0.8000, 1.124,  0.0000, 0.0000, 1.1100, 0.1600, 1.820e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Cadmium",
    8.639, 0.0, 0.0, 0.0, 0.2480, 1.640,  0.0000, 0.0000, 2.5000, 0.5500, 0.231e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Copper (OFHC 1/2 Hard)",
    8.930, 0.0, 0.0, 0.0, 0.3940, 1.489,  0.0000, 0.0000, 2.0200, 0.4700, 0.383e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Gold",
   19.300, 0.0, 0.0, 0.0, 0.3080, 1.560,  0.0000, 0.0000, 2.9900, 0.5900, 0.128e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Gold - 5% Copper",
   18.100, 0.0, 0.0, 0.0, 0.3050, 1.560,  0.0000, 0.0000, 2.9900, 0.5900, 0.141e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Hafnium",
   13.300, 0.0, 0.0, 0.0, 0.2980, 1.448, -3.8520, 7.0000, 1.0900, 0.1500, 0.144e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Lead",
   11.340, 0.0, 0.0, 0.0, 0.2006, 1.429,  0.8506,-1.6400, 2.7400, 0.5400, 0.124e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Lithium",/*see Note*/
    0.534, 0.0, 0.0, 0.0, 0.4770, 1.066,  0.0000, 0.0000, 0.9200, 0.3800, 4.385e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Magnesium (AZ31B-H24)",
    1.780, 0.0, 0.0, 0.0, 0.4520, 1.242,  0.0000, 0.0000, 1.5400, 0.3300, 0.999e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Mercury",/*see Note*/
   13.541, 0.0, 0.0, 0.0, 0.1451, 3.687, -9.7820,13.9500, 2.7400, 0.5700, 0.036e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Molybdenum",
   10.200, 0.0, 0.0, 0.0, 0.5143, 1.255,  0.0000, 0.0000, 1.5900, 0.3000, 0.243e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Nickle",
    8.900, 0.0, 0.0, 0.0, 0.4650, 1.445,  0.0000, 0.0000, 1.9300, 0.5000, 0.401e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Nickel Alloy Monel",
    8.810, 0.0, 0.0, 0.0, 0.4190, 1.540,  0.0000, 0.0000, 1.9500, 0.4900, 0.411e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Niobium",
    8.590, 0.0, 0.0, 0.0, 0.4440, 1.207,  0.0000, 0.0000, 1.6600, 0.3400, 0.257e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Niobium Alloy Fanstell 85",
   10.770, 0.0, 0.0, 0.0, 0.4160, 1.195,  0.0000, 0.0000, 1.6300, 0.3600, 0.210e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Platinum",
   21.440, 0.0, 0.0, 0.0, 0.3640, 1.540,  0.0000, 0.0000, 2.7400, 0.5800, 0.128e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Platinum - 20% Iridium",
   21.670, 0.0, 0.0, 0.0, 0.3660, 1.530,  0.0000, 0.0000, 2.7400, 0.5800, 0.128e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Silver",
   10.490, 0.0, 0.0, 0.0, 0.3270, 1.550,  0.0000, 0.0000, 2.4000, 0.5600, 0.233e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Stainless Steel (21-6-9)",
    7.795, 0.0, 0.0, 0.0, 0.4440, 2.200, -2.9050, 2.5440, 1.9300, 0.5000, 0.426e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Stainless Steel (304)",
    7.900, 0.0, 0.0, 0.0, 0.4570, 1.490,  0.0000, 0.0000, 1.9300, 0.5000, 0.423e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Steel (Vascomax 250)",
    8.129, 0.0, 0.0, 0.0, 0.3980, 1.580,  0.0000, 0.0000, 1.6000, 0.5000, 0.408e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Steel (4340 RC 38)",
    7.810, 0.0, 0.0, 0.0, 0.4578, 1.330,  0.0000, 0.0000, 1.6700, 0.4300, 0.448e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Tantalum",
   16.690, 0.0, 0.0, 0.0, 0.3410, 1.200,  0.0000, 0.0000, 1.6700, 0.4200, 0.135e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Tantalum - 10% Tungsten",
   16.960, 0.0, 0.0, 0.0, 0.3460, 1.200,  0.0000, 0.0000, 1.6700, 0.4200, 0.135e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Thorium",
   11.700, 0.0, 0.0, 0.0, 0.2130, 1.278,  0.0000, 0.0000, 1.4500, 0.2100, 0.112e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Thulium",/*see Note*/
    9.289, 0.0, 0.0, 0.0, 0.2260, 0.737,  0.4739, 0.2361, 1.0800, 0.0000, 0.100e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Tin",
    7.287, 0.0, 0.0, 0.0, 0.2590, 1.490,  0.0000, 0.0000, 2.2700, 0.5400, 0.220e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Titanium",
    4.510, 0.0, 0.0, 0.0, 0.5020, 1.536, -5.1380,10.8200, 1.2300, 0.1700, 0.500e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Titanium - 6%Al - 4%V",
    4.419, 0.0, 0.0, 0.0, 0.5130, 1.028,  0.0000, 0.0000, 1.2300, 0.1700, 0.525e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Tungsten",
   19.300, 0.0, 0.0, 0.0, 0.4030, 1.237,  0.0000, 0.0000, 1.6700, 0.3800, 0.129e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Tungsten - 3.5%Ni - 1.5%Fe",
   18.167, 0.0, 0.0, 0.0, 0.4030, 1.237,  0.0000, 0.0000, 1.6700, 0.3800, 0.143e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Uranium",
   19.050, 0.0, 0.0, 0.0, 0.2480, 1.530,  0.0000, 0.0000, 2.3200, 0.5700, 0.108e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Uranium - 0.75% Titanium",
   18.620, 0.0, 0.0, 0.0, 0.2567, 1.619,  0.0000, 0.0000, 2.3200, 0.5700, 0.111e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Uranium - 5% Molybdenum",
   18.170, 0.0, 0.0, 0.0, 0.2590, 1.560,  0.0000, 0.0000, 2.3200, 0.5700, 0.115e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Uranium - 7.5%Nb - 2.5%Zr",
   16.450, 0.0, 0.0, 0.0, 0.2570, 1.500,  0.0000, 0.0000, 1.9000, 0.4500, 0.113e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Vanadium",
    6.100, 0.0, 0.0, 0.0, 0.5077, 1.201,  0.0000, 0.0000, 1.4000, 0.2200, 0.464e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Zinc",
    7.139, 0.0, 0.0, 0.0, 0.3030, 1.550,  0.0000, 0.0000, 2.2400, 0.5200, 0.389e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Zirconium",
    6.506, 0.0, 0.0, 0.0, 0.3890, 0.292,  2.7520,-1.9750, 0.9300, 0.1200, 0.269e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Graphite",/*see Note*/
    2.200, 0.0, 0.0, 0.0, 0.3900, 2.160,  1.5400,-9.4300, 0.2400, 0.0000, 2.775e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Kel-F (C2F3CL)",/*see Note*/
    2.133, 0.0, 0.0, 0.0, 0.2050, 1.660,  0.4064,-1.0370, 0.6600, 0.0000, 1.264e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Lexan (Polycarbonate)",/*see Note*/
    1.196, 0.0, 0.0, 0.0, 0.1933, 3.490, -8.2000, 9.6000, 0.6100, 0.0000, 0.100e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Lithium Fluoride (100 Single Crystal)",
    2.638, 0.0, 0.0, 0.0, 0.5150, 1.350,  0.0000, 0.0000, 1.6900, 0.3400, 1.560e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Lucite (PMMA)",
    1.182, 0.0, 0.0, 0.0, 0.2180, 2.088, -1.1240, 0.0000, 0.8500, 0.0000, 1.460e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Micarta",/*see Note*/
    1.395, 0.0, 0.0, 0.0, 0.2030, 3.784, -7.6900, 6.7200, 0.3000, 0.0000, 0.100e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Nylon (Type 6/6)",/*see Note*/
    1.140, 0.0, 0.0, 0.0, 0.2208, 2.830, -3.1800, 1.1770, 0.7000, 0.0000, 0.100e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Polyethylene (CH2)",/*see Note*/
    0.952, 0.0, 0.0, 0.0, 0.3050, 1.328,  1.1390,-1.8710, 0.6700, 0.0000, 2.660e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Polypentene (CH2)",/*see Note*/
    0.083, 0.0, 0.0, 0.0, 0.1800, 2.410, -2.2170, 1.2000, 0.5800, 0.0000, 0.100e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Polystyrene (CH)",/*see Note*/
    1.046, 0.0, 0.0, 0.0, 0.1890, 2.965, -4.0690, 2.3280, 0.6700, 0.0000, 1.820e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Silastic (APC 2.5)",/*see Note*/
    1.030, 0.0, 0.0, 0.0, 0.0610, 4.512, -7.7200, 4.6580, 1.0000, 0.0000, 0.100e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Teflon",/*see Note*/
    2.150, 0.0, 0.0, 0.0, 0.1680, 1.123,  3.9830,-5.7970, 0.5900, 0.0000, 1.718e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Tungsten (98.52%) in Plastic Binder",/*see Note*/
   15.180, 0.0, 0.0, 0.0, 0.1030, 7.637,-25.9160,30.8360, 1.2600, 0.0000, 0.129e-5,     0.0,   0.0,   0.0,   0.0, NULL,
"Tungsten Carbide",
   14.900, 0.0, 0.0, 0.0, 0.5190, 1.160,  0.0000, 0.0000, 1.5000, 0.0000, 0.160e-5,     0.0,   0.0,   0.0,   0.0, NULL,
NULL,
    0.000, 0.0, 0.0, 0.0, 0.0000, 0.000,  0.0000, 0.0000, 0.0000, 0.0000, 0.000000,     0.0,   0.0,   0.0,   0.0, NULL,
/*    name          ,
*    Rho0,  P0,  E0,  T0,     C0,    S1,      S2,     S3, gamma0,      b,       CV, Rho_max, V_min, P_inf, E_inf, TrefSpline
*/
};


LOCAL	boolean MG_OOR;

	/* LOCAL Function Prototypes */
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
LOCAL	double	MG_pressure(Locstate);
LOCAL	double	MG_sound_speed_squared(Locstate);
LOCAL	double	MG_specific_internal_energy(Locstate);

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
LOCAL	double	MG_temperature(Locstate);
LOCAL	double	MG_entropy(Locstate);
LOCAL	double	MG_gruneisen_gamma(Locstate);
LOCAL	double	MG_fundamental_derivative(Locstate);
LOCAL	double	MG_C_V(Locstate);
LOCAL	double	MG_C_P(Locstate);
LOCAL	double	MG_K_T(Locstate);

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
LOCAL	void    MG_single_eos_load_pressure_and_sound_speed2(Vec_Gas*,int,int);
LOCAL	void	MG_single_eos_load_pressure_and_gammas(Vec_Gas*,int,int);
LOCAL	void	MG_single_eos_load_pressure(Vec_Gas*,int,int);
LOCAL	void	MG_single_eos_load_sound_speed2(Vec_Gas*,int,int);

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
LOCAL	double	MG_dens_Hugoniot(double,Locstate);

	/* General Wave Curve Functions */
LOCAL	double	MG_mass_flux(double,Locstate);
LOCAL	double	MG_mass_flux_squared(double,Locstate);

	/* Functions for the Evaluation of Riemann Solutions */
LOCAL	double	MG_oned_fan_state(double,Locstate,Locstate,Locstate,int,boolean*);

	/* Functions to Compute Riemann Solutions */
LOCAL	double	MG_riemann_wave_curve(Locstate,double);

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
LOCAL	double	MG_dens_rarefaction(double,Locstate);
LOCAL	double	MG_pressure_rarefaction(double,Locstate);

	/* INITIALIZATION UTILITY FUNCTIONS */
LOCAL	void	MG_fprint_EOS_params(FILE*,Gas_param*);
LOCAL	void	MG_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,Gas_param*);
LOCAL	void	MG_prompt_for_EOS_params(INIT_DATA*,Gas_param*,
					 const char*,const char*);

	/* Utility Functions */
LOCAL	boolean	test_out_of_range(double,Locstate,const char*,boolean);
LOCAL	double	MG_pressure_on_adiabat(double,double,double,MG_params*);
LOCAL	double	MG_real_sound_speed_dens(double,double,double,MG_params*);
LOCAL	double	MG_sound_speed_squared_on_adiabat(double,double,double,
	                                          MG_params*);
LOCAL	double	Gamma(double,MG_params*);
LOCAL	double	dGamma(double,MG_params*);
LOCAL	double	d2Gamma(double,MG_params*);
LOCAL	double	P_ref(double,MG_params*);
LOCAL	double	dP_ref(double,MG_params*);
LOCAL	double	d2P_ref(double,MG_params*);
LOCAL	double	E_ref(double,MG_params*);
LOCAL	double	dE_ref(double,MG_params*);
LOCAL	double	d2E_ref(double,MG_params*);
LOCAL	double	T_ref(double,MG_params*);
LOCAL	double   IrhoG(double,double,MG_params*);
LOCAL	double	pr_Riemann_function(double,Locstate);
LOCAL	double	int_c_over_rho_drho(double,Locstate);
LOCAL	double	max_compression(MG_params*);
LOCAL	void	extend_spline_Tref_curve(double,MG_params*);
LOCAL	void	print_MG_params(MG_params*);
LOCAL	void	initialize_Tref(MG_params*);
LOCAL	void	set_eos_function_hooks(EOS*);

	/* Equation of state domain functions */
LOCAL	double	MG_Min_energy(Locstate);
LOCAL	double	MG_Min_pressure(Locstate);
LOCAL	double	MG_Min_temperature(Locstate);


EXPORT	EOS	*set_MG_eos(
	EOS	*eos)
{
	if (eos == NULL)
	    scalar(&eos,sizeof(MG_EOS));
	(void) set_GENERIC_eos(eos);
	set_eos_function_hooks(eos);
	MG_OOR = (debugging("no_MG_OOR")) ? NO : YES;
	return eos;
}


LOCAL	void	set_eos_function_hooks(
	EOS *eos)
{
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	eos->_pressure = MG_pressure;
	eos->_sound_speed_squared = MG_sound_speed_squared;
	eos->_specific_internal_energy = MG_specific_internal_energy;

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	eos->_temperature = MG_temperature;
	eos->_entropy = MG_entropy;
	eos->_gruneisen_gamma = MG_gruneisen_gamma;
	eos->_fundamental_derivative = MG_fundamental_derivative;
	eos->_C_V = MG_C_V;
	eos->_C_P = MG_C_P;
	eos->_K_T = MG_K_T;

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
	eos->_single_eos_load_pressure_and_sound_speed2 =
	    MG_single_eos_load_pressure_and_sound_speed2;
	eos->_single_eos_load_pressure_and_gammas =
	    MG_single_eos_load_pressure_and_gammas;
	eos->_single_eos_load_pressure = MG_single_eos_load_pressure;
	eos->_single_eos_load_sound_speed2 = MG_single_eos_load_sound_speed2;

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
	eos->_dens_Hugoniot = MG_dens_Hugoniot;

	/* General Wave Curve Functions */
	eos->_mass_flux = MG_mass_flux;
	eos->_mass_flux_squared = MG_mass_flux_squared;

	/* Functions for the Evaluation of Riemann Solutions */
	eos->_oned_fan_state = MG_oned_fan_state;

	/* Functions to Compute Riemann Solutions */
	eos->_riemann_wave_curve = MG_riemann_wave_curve;

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
	eos->_dens_rarefaction = MG_dens_rarefaction;
	eos->_pressure_rarefaction = MG_pressure_rarefaction;

	/* INITIALIZATION UTILITY FUNCTIONS */
	eos->_fprint_EOS_params = MG_fprint_EOS_params;
	eos->_read_print_EOS_params = MG_read_print_EOS_params;
	eos->_prompt_for_EOS_params = MG_prompt_for_EOS_params;

	/* Equation of state domain functions */
	eos->_Min_energy = MG_Min_energy;
	eos->_Min_pressure = MG_Min_pressure;
}


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/


/*
*			MG_pressure():
*
*	Returns the thermodynamic pressure of a state.
*
*				     dE  |
*			     P = -  ---- |
*		                     dV  |S
*
*	Where E = specific internal energy,  V = specific volume,  and
*	S = specific entropy.
*/

LOCAL	double	MG_pressure(
	Locstate state)
{
    	MG_params *mg;
	double pr, rho, e, T;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double pmin;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	if (is_obstacle_state(state))
	    return 0.0;

	rho = Dens(state);
	mg = &MG_Eos(state)->MGparams;
	test_out_of_range(rho,state,"MG_pressure",YES);

	switch (state_type(state)) 
	{
	case GAS_STATE:
	    e = (Energy(state) - kinetic_energy(state))/rho;
	    pr = P_ref(rho,mg) + rho*Gamma(rho,mg)*(e - E_ref(rho,mg));
	    break;

	case EGAS_STATE:
	    e = Energy(state);
	    pr = P_ref(rho,mg) + rho*Gamma(rho,mg)*(e - E_ref(rho,mg));
	    break;

	case TGAS_STATE:
	case VGAS_STATE:
	    pr = Press(state);
	    break;

	case FGAS_STATE:
	    T = Temperature(state);
	    pr = P_ref(rho,mg) + mg->CV*rho*Gamma(rho,mg)*(T - T_ref(rho,mg));
	    break;

	default:
	    screen("ERROR in MG_pressure(), no such or unimplemented "
		   "state type\n");
	    clean_up(ERROR);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pmin = Min_pressure(state);
	if (pr < pmin)
	    pr = pmin;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pr;
}		/*end MG_pressure*/


/*
*			MG_sound_speed_squared():
*
*	Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*			c = ---- |
*                           drho |S
*/

LOCAL	double	MG_sound_speed_squared(
	Locstate state)
{
    	MG_params *mg;
	double c, c2ref, c2, rho, e, p, T, Tr, Gm, dGm;
	double pr, dpr, er, der;

	if (is_obstacle_state(state))
	    return 0.0;
	if (state_type(state) == VGAS_STATE) 
	{
	    c = Sound_speed(state);
	    return c*c;
	}

	mg = &MG_Eos(state)->MGparams;
	rho = Dens(state);
	test_out_of_range(rho,state,"MG_sound_speed_squared",YES);
	Gm = Gamma(rho,mg);
	dGm = dGamma(rho,mg);
	pr = P_ref(rho,mg);
	dpr = dP_ref(rho,mg);
	der = dE_ref(rho,mg);
	c2ref = (Gm/rho)*(pr+der) - dpr/(rho*rho);

	switch (state_type(state)) 
	{
	case GAS_STATE:
	    e = (Energy(state) - kinetic_energy(state))/rho;
	    er = E_ref(rho,mg);
	    c2 = (e - er)*(Gm*(Gm+1.0) - dGm/rho) + c2ref;
	    break;

	case EGAS_STATE:
	    e = Energy(state);
	    er = E_ref(rho,mg);
	    c2 = (e - er)*(Gm*(Gm+1.0) - dGm/rho) + c2ref;
	    break;

	case TGAS_STATE:
	    p = Press(state);
	    c2 = (p - pr)*((Gm+1.0) - dGm/(Gm*rho))/rho + c2ref;
	    break;

	case FGAS_STATE:
	    T = Temperature(state);
	    Tr = T_ref(rho,mg);
	    c2 = mg->CV*(T - Tr)*(Gm*(Gm+1.0) - dGm/rho) + c2ref;
	    break;

	default:
	    screen("ERROR in MG_sound_speed_squared(), "
		    "no such, or unimplemented, state type %d\n",
		    state_type(state));
	    clean_up(ERROR);
	}
	if (c2 < 0.0)
	{
	    double emin, pmin, Tmin;
	    switch (state_type(state)) 
	    {
	    case GAS_STATE:
	        e = (Energy(state) - kinetic_energy(state));
	        er = E_ref(rho,mg);
		emin = MG_Min_energy(state);
		if (e < emin)
		    e = emin;
	        c2 = (e/rho - er)*(Gm*(Gm+1.0) - dGm/rho) + c2ref;
	        break;
	    case EGAS_STATE:
	        e = Energy(state);
	        er = E_ref(rho,mg);
		emin = MG_Min_energy(state);
		if (rho*e < emin)
		    e = emin/rho;
	        c2 = (e - er)*(Gm*(Gm+1.0) - dGm/rho) + c2ref;
	        break;
	    case TGAS_STATE:
	        p = Press(state);
		pmin = MG_Min_pressure(state);
		if (p < pmin)
		    p = pmin;
	        c2 = (p - pr)*((Gm+1.0) - dGm/(Gm*rho))/rho + c2ref;
	        break;
	    case FGAS_STATE:
	        T = Temperature(state);
		Tmin = MG_Min_temperature(state);
		if (T < Tmin)
		    T = Tmin;
	        Tr = T_ref(rho,mg);
	        c2 = mg->CV*(T - Tr)*(Gm*(Gm+1.0) - dGm/rho) + c2ref;
	        break;
	    default:
		break;
	    }

	    if ((c2 < 0.0) && MG_OOR)
	    {
	        boolean bin_out;
	        screen("ERROR in MG_sound_speed_squared(), "
		       "imaginary sound speed\n");
	        (void) printf("c2 = %"FFMT"\n",c2);
	        switch (state_type(state)) 
	        {
	        case GAS_STATE:
	            e = (Energy(state) - kinetic_energy(state))/rho;
	            (void) printf("GAS_STATE - rho = %"FFMT", e = %"FFMT"\n",
				  rho,e); 
	            break;
	        case EGAS_STATE:
	            e = Energy(state);
	            (void) printf("EGAS_STATE - rho = %"FFMT", e = %"FFMT"\n",
				  rho,e);
	            break;
	        case TGAS_STATE:
	            p = Press(state);
	            (void) printf("TGAS_STATE - rho = %"FFMT", p = %"FFMT"\n",
				  rho,e);
	            break;
	        case FGAS_STATE:
	            T = Temperature(state);
	            (void) printf("FGAS_STATE - rho = %"FFMT", T = %"FFMT"\n",
				  rho,e);
	            break;
	        default:
		    break;
	        }
	        bin_out = is_binary_output();
	        set_binary_output(NO);
	        fprint_Gas_param(stdout,Params(state));
	        set_binary_output(bin_out);
	        clean_up(ERROR);
	    }
	    c2 = fabs(c2);
	}
	return c2;
}		/*end MG_sound_speed_squared*/


/*
*			MG_specific_internal_energy():
*
*	Returns the specific internal energy = internal energy per unit
*	mass of the state.
*/

LOCAL	double	MG_specific_internal_energy(
	Locstate state)
{
    	MG_params *mg;
	double rho, e;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double Emin;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	if (is_obstacle_state(state))
	    return 0.0;

	rho = Dens(state);
	mg = &MG_Eos(state)->MGparams;
	test_out_of_range(rho,state,"MG_specific_internal_energy",YES);
	switch (state_type(state)) 
	{
	case GAS_STATE:
	    e = (Energy(state) - kinetic_energy(state))/rho;
	    break;

	case EGAS_STATE:
	    e = Energy(state);
	    break;

	case TGAS_STATE:
	    e = E_ref(rho,mg) + (Press(state)-P_ref(rho,mg))/(rho*Gamma(rho,mg));
	    break;

	case VGAS_STATE:
	    return Int_en(state);

	case FGAS_STATE:
	    e = E_ref(rho,mg) + mg->CV*(Temperature(state)-T_ref(rho,mg));
	    break;

	default:
	    screen("ERROR in MG_specific_internal_energy(), "
		   "no such or unimplemented state type\n");
	    clean_up(ERROR);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Emin = Min_energy(state);
	if (rho*e < Emin)
	    e = Emin/rho;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return e;
}		/*end MG_specific_internal_energy*/


/***************END PRIMARY THERMODYNAMIC FUNCTIONS ************************/
/***************SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS ***********/


/*
*			MG_temperature():
*
*	Returns the thermodynamic temperature of a state.
*
*                            dE |
*			T = --- |
*                            dS |V
*/

LOCAL	double	MG_temperature(
	Locstate state)
{
    	MG_params *mg;
    	double T, p, e;
	double rho;

	if (is_obstacle_state(state))
	    return HUGE_VAL;

	mg = &MG_Eos(state)->MGparams;
	rho = Dens(state);
	test_out_of_range(rho,state,"MG_temperature",YES);
	switch (state_type(state)) 
	{
	case GAS_STATE:
	    e = (Energy(state) - kinetic_energy(state))/rho;
	    T = T_ref(rho,mg) + (e - E_ref(rho,mg))/mg->CV;
	    break;

	case EGAS_STATE:
	    e = Energy(state);
	    T = T_ref(rho,mg) + (e - E_ref(rho,mg))/mg->CV;
	    break;

	case VGAS_STATE:
#if defined(VERBOSE_GAS_PLUS)
	    T = Temp(state);
#else /*defined(VERBOSE_GAS_PLUS)*/
	    e = Int_en(state);
	    T = T_ref(rho,mg) + (e - E_ref(rho,mg))/mg->CV;
#endif /*defined(VERBOSE_GAS_PLUS)*/
	    break;

	case TGAS_STATE:
	    p = Press(state);
	    T = T_ref(rho,mg) + (p - P_ref(rho,mg))/(mg->CV*rho*Gamma(rho,mg));
	    break;

	case FGAS_STATE:
	    T = Temperature(state);
	    break;

	default:
	    screen("ERROR in MG_temperature(), no such or unimplemented "
		   "state type\n");
	    clean_up(ERROR);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (T < MACH_EPS) /*TOLERANCE*/
	    T = MACH_EPS;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return T;
}		/*end MG_temperature*/

/*
*			MG_entropy():
*
*	Returns the specific entropy of a state.
*/

LOCAL	double	MG_entropy(
	Locstate state)
{
    	MG_params *mg;
    	double S, T;
	double rho;

	if (is_obstacle_state(state))
	    return HUGE_VAL;

	rho = Dens(state);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Dens(state) < Vacuum_dens(state))
	    rho = Vacuum_dens(state);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	mg = &MG_Eos(state)->MGparams;
	test_out_of_range(rho,state,"MG_entropy",YES);
	switch (state_type(state)) 
	{
	case GAS_STATE:
	case EGAS_STATE:
	case TGAS_STATE:
	case FGAS_STATE:
	    T = temperature(state);
	    S = mg->CV*log(IrhoG(mg->Rho0,rho,mg)*T);
	    break;

	case VGAS_STATE:
	    S = Entropy(state);
	    break;

	default:
	    screen("ERROR in MG_entropy(), no such or unimplemented "
		   "state type\n");
	    clean_up(ERROR);
	}

	return S;
}		/*end MG_entropy*/


/*
*			MG_gruneisen_gamma():
*
*	Returns the dimensionless Gruneisen exponent
*
*
*                                                 dP/dE |
*		GAMMA = - d(log T)/d(log V) |  =  -----  V
*                                            S     rho
*
*	As usual P = thermodynamic pressure,  V = specific volume
*	rho = density, E = specific internal energy,
*	and  S = specific entropy.
*
*
*/

LOCAL	double	MG_gruneisen_gamma(
	Locstate state)
{
    	double rho;
    	MG_params *mg;

	if (is_obstacle_state(state))
	    return HUGE_VAL;

	rho = Dens(state);
	mg = &MG_Eos(state)->MGparams;
	test_out_of_range(rho,state,"MG_gruneisen_gamma",YES);
	return Gamma(rho,mg);
}		/*end MG_gruneisen_gamma*/

/*
*			MG_fundamental_derivative():
*
*	Returns the fundamental derivative of gas dynamics for the state.
*	This quantity is defined by the formula
*
*			    2      2
*		           d P / dV  |
*                                    |S
*             G = -0.5 V -----------------
*                          dP / dV |
*                                  |S
*
*	Where P is the thermodynamic pressure,  V is the specific volume
*	and S is the specific entropy.  Both derivatives are taken at
*	constant S.
*/

LOCAL	double	MG_fundamental_derivative(
	Locstate state)
{
    	MG_params *mg;
	double e, G, rho, c2, c2r, Gr;
	double pr, dpr, d2pr, er, der, d2er, Gm, dGm, d2Gm;
	double k, dk;

	if (is_obstacle_state(state))
	    return 0.0;

	rho = Dens(state);
	mg = &MG_Eos(state)->MGparams;
	test_out_of_range(rho,state,"MG_fundamental_derivative",YES);

	e = MG_specific_internal_energy(state);
	Gm = Gamma(rho,mg);
	dGm = dGamma(rho,mg);
	d2Gm = d2Gamma(rho,mg);
	c2 = MG_sound_speed_squared(state);
	pr = P_ref(rho,mg);
	dpr = dP_ref(rho,mg);
	d2pr = d2P_ref(rho,mg);
	er = E_ref(rho,mg);
	der = dE_ref(rho,mg);
	d2er = d2E_ref(rho,mg);
	k = pr + der;
	dk = dpr + d2er;
	c2r = (Gm/rho)*k - dpr/(rho*rho);
	if (c2r > 0.0)
	{
	    Gr = 1.0 + 0.5*Gm +
	        (0.5/(rho*rho*c2r))*(-2.0*dGm*k-Gm*dk+(Gm+2.0)*dpr + d2pr/rho);

	    G = Gr + ((e - er)/(rho*rho*2.0*c2*c2r))*(
		    k*(2.0*Gm*dGm - 2.0*dGm*dGm/rho + Gm*d2Gm/rho) +
		    dk*(Gm*Gm*(Gm + 1.0) - Gm*dGm/rho) -
		    dpr*(Gm*(2.0 + Gm*(3.0 + Gm) - 3.0*dGm/rho) -
		         2.0*dGm/rho + d2Gm/(rho*rho)) -
		    d2pr*((Gm*(Gm+1.0) - dGm/rho)/rho));
	}
	else
	    G = 0.5*(Gm + 2.0);
	return G;
}		/*end MG_fundamental_derivative*/

/*
*			MG_C_V():
*
*	Specific heat at constant volume.
*
*                        dS  |
*		C_V = T ---- |
*                        dT  | V
*/

LOCAL	double	MG_C_V(
	Locstate state)
{
    	MG_params *mg;

	if (is_obstacle_state(state))
	    return 0.0;

	mg = &MG_Eos(state)->MGparams;
	test_out_of_range(Dens(state),state,"MG_C_V",YES);
	return mg->CV;
}	/* end MG_C_V */

/*
*			MG_C_P():
*
*	Specific heat at constant pressure.
*
*
*                        dS  |
*		C_P = T ---- |
*                        dT  | P
*/

LOCAL	double	MG_C_P(
	Locstate state)
{
    	double Cp, Cv;
	double T, rho;
	double pr, dpr, Tr, der, Gm, dGm;
    	MG_params *mg;

	if (is_obstacle_state(state))
	    return 0.0;

	mg = &MG_Eos(state)->MGparams;
	rho = Dens(state);
	test_out_of_range(rho,state,"MG_C_P",YES);
	Cv = mg->CV;
	T = temperature(state);
	Gm = Gamma(rho,mg);
	dGm = dGamma(rho,mg);
	pr = P_ref(rho,mg);
	dpr = dP_ref(rho,mg);
	Tr = T_ref(rho,mg);
	der = dE_ref(rho,mg);
	Cp = Cv*(1.0 +
		Cv*Gm*Gm*rho*T/
	        (Gm*(pr+der-Cv*Gm*rho*Tr)-dpr/rho+(Gm-dGm/rho)*Cv*rho*(T-Tr)));
	return Cp;
}	/* end MG_C_P */

/*
*			MG_K_T():
*
*	Isothermal compressibility.
*
*                        1   dV  |
*		K_T = - --- ---- |
*                        V   dP  | T
*/

LOCAL	double	MG_K_T(
	Locstate state)
{
    	double K_T, Cv;
	double rho, e, er, der, pr, dpr, Tr;
	double Gm, dGm;
    	MG_params *mg;

	if (is_obstacle_state(state))
	    return 0.0;

	rho = Dens(state);
	mg = &MG_Eos(state)->MGparams;
	test_out_of_range(rho,state,"MG_K_T",YES);
	e = MG_specific_internal_energy(state);
	Cv = mg->CV;
	pr = P_ref(rho,mg);
	dpr = dP_ref(rho,mg);
	er = E_ref(rho,mg);
	der = dE_ref(rho,mg);
	Tr = T_ref(rho,mg);
	Gm = Gamma(rho,mg);
	dGm = dGamma(rho,mg);

	K_T = rho/
	    (Gm*rho*(pr+der-Cv*Gm*rho*Tr)-dpr-(rho*dGm-rho*rho*Gm)*(e-er));

    	return K_T;
}	/* end MG_K_T */


/***************END SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS *******/
/***************VECTORIZED THERMODYNAMIC FUNCTIONS *************************/

/*
*		MG_single_eos_load_pressure_and_sound_speed2():
*
*	Loads a vector of pressures and sound speeds into the
*	appropriate fields of the Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	MG_single_eos_load_pressure_and_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	double     *rho = vst->rho + offset;
	double     *p = vst->p + offset;
	double     *c2 = vst->c2 + offset;
	double     *e = vst->e + offset;
	double     *FD;
	double     pr, dpr, d2pr, er, der, d2er, Gm, dGm, d2Gm;
	double     c2r, Gr;
	double     k, dk;
	double     Rho_max;
	int       i, istart, iend;
    	MG_params *mg;

	FD = (vst->FD != NULL) ? vst->FD + offset : NULL;

	mg = &MG_Eos(vst->state[offset])->MGparams;
	Rho_max = mg->Rho_max;
	for (istart = 0; istart < vsize; ++istart)
	    if (rho[istart] <= Rho_max)
		break;
	if (istart == vsize)
	{
	    if (MG_OOR)
	    {
	        screen("ERROR in "
		       "MG_single_eos_load_pressure_and_sound_speed2(), "
		       "all densities out of EOS range\n");
	        (void) printf("vsize = %d\n",vsize);
	        (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	        for (istart = 0; istart < vsize; ++istart)
	            (void) printf("    rho[%d] = %"FFMT"\n",istart,rho[istart]);
	        (void) printf("\n");
	        clean_up(ERROR);
	    }
	    else
		istart = 0;
	}
	for (iend = vsize-1; iend >= 0; --iend)
	    if (rho[iend] <= Rho_max)
		break;
	++iend;
	if (iend == 0)
	{
	    if (MG_OOR)
	    {
	        screen("ERROR in "
		       "MG_single_eos_load_pressure_and_sound_speed2(), "
		       "all densities out of EOS range\n");
	        (void) printf("vsize = %d\n",vsize);
	        (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	        for (iend = vsize-1; iend >= 0; --iend)
	            (void) printf("    rho[%d] = %"FFMT"\n",iend,rho[iend]);
	        (void) printf("\n");
	        clean_up(ERROR);
	    }
	    else
		iend = vsize;
	}
	for (i = istart; i < iend; ++i)
	{
    	    if (rho[i] > Rho_max)
	    {                       
		if ((istart < i) && ((i+1) < iend))
		{
		    if ((rho[i-1] < Rho_max) && (rho[i+1] < Rho_max))
		    {
			rho[i] = 0.5*(rho[i-1] + rho[i+1]);
		    }
		    else
			rho[i] = min(rho[i-1],rho[i+1]);
		}
		else if (MG_OOR)
		{
		    int j;
	            screen("ERROR in "
		           "MG_single_eos_load_pressure_and_sound_speed2(), "
		           "interior region state out of range\n");
	            (void) printf("vsize = %d\n",vsize);
	            (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	            (void) printf("rho[%d] = %"FFMT"\n",i,rho[i]);
		    (void) printf("istart = %d, iend = %d\n",istart,iend);
	            for (j = 0; j < vsize; ++j)
	                (void) printf("    rho[%d] = %"FFMT"\n",j,rho[j]);
	            (void) printf("\n");
		    clean_up(ERROR);
		}
	    }
	    pr = P_ref(rho[i],mg);
	    dpr = dP_ref(rho[i],mg);
	    er = E_ref(rho[i],mg);
	    der = dE_ref(rho[i],mg);
	    Gm = Gamma(rho[i],mg);
	    dGm = dGamma(rho[i],mg);
	    p[i] = pr + rho[i]*Gm*(e[i] - er);
	    c2r = (Gm/rho[i])*(pr+der) - dpr/(rho[i]*rho[i]);
	    c2[i] = (e[i] - er)*(Gm*(Gm+1.0) - dGm/rho[i]) + c2r;
	    if (c2[i] < 0.0)
		c2[i] = -c2[i];
	    if (FD != NULL)
	    {
	        if (c2r > 0.0)
	        {
	            d2pr = d2P_ref(rho[i],mg);
	            d2er = d2E_ref(rho[i],mg);
	            d2Gm = d2Gamma(rho[i],mg);
	            k = pr + der;
	            dk = dpr + d2er;
	            Gr = 1.0 + 0.5*Gm +
			(0.5/(rho[i]*rho[i]*c2r))*
			(-2.0*dGm*k-Gm*dk+(Gm+2.0)*dpr + d2pr/rho[i]);

	            FD[i] = Gr + ((e[i]-er)/(rho[i]*rho[i]*2.0*c2[i]*c2r))*(
		        k*(2.0*Gm*dGm-2.0*dGm*dGm/rho[i] + Gm*d2Gm/rho[i]) +
		        dk*(Gm*Gm*(Gm + 1.0) - Gm*dGm/rho[i]) -
		        dpr*(Gm*(2.0 + Gm*(3.0 + Gm) - 3.0*dGm/rho[i]) -
		        2.0*dGm/rho[i] + d2Gm/(rho[i]*rho[i])) -
		        d2pr*((Gm*(Gm+1.0) - dGm/rho[i])/rho[i]));
	        }
	        else
	            FD[i] = 0.5*(Gm + 2.0);
	    }
	}
	/*
	*   When the states at the end of the uni_arrays are out of range
	*   it is usually due to irregular stencil points that will be
	*   overwritten with their correct values in the irregular interior
	*   state sweep.  Copy out the adjacent interior values to these
	*   temporary values.
	*/
	for (i = 0; i < istart; ++i)
	{
	    p[i] = p[istart];
	    c2[i] = c2[istart];
	    if (FD != NULL)
		FD[i] = FD[istart];
	}
	for (i = iend; i < vsize; ++i)
	{
	    p[i] = p[iend-1];
	    c2[i] = c2[iend-1];
	    if (FD != NULL)
		FD[i] = FD[iend-1];
	}
}		/*end MG_single_eos_load_pressure_and_sound_speed2*/


/*
*		MG_single_eos_load_pressure_and_gammas():
*
*	Loads the pressure, adiabatic exponent, and Gruneisen
*	coefficient uni_arrays of the Vec_Gas state vst.
*	This function assumes that the specific internal energy
*	uni_array vst->e is already loaded.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_gammas.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	MG_single_eos_load_pressure_and_gammas(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	double     *rho = vst->rho + offset;
	double     *p = vst->p + offset;
	double     *c2 = vst->c2 + offset;
	double     *GAM = vst->GAM + offset;
	double     *e = vst->e + offset;
	double     *FD;
	double     pr, dpr, d2pr, er, der, d2er, Gm, dGm, d2Gm;
	double     c2r, Gr;
	double     k, dk;
	double     Rho_max;
	int       i, istart, iend;
    	MG_params *mg;

	FD = (vst->FD != NULL) ? vst->FD + offset : NULL;

	mg = &MG_Eos(vst->state[offset])->MGparams;
	Rho_max = mg->Rho_max;
	for (istart = 0; istart < vsize; ++istart)
	    if (rho[istart] <= Rho_max)
		break;
	if (istart == vsize)
	{
	    if (MG_OOR)
	    {
	        screen("ERROR in MG_single_eos_load_pressure_and_gammas(), "
		       "all densities out of EOS range\n");
	        (void) printf("vsize = %d\n",vsize);
	        (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	        for (istart = 0; istart < vsize; ++istart)
	            (void) printf("    rho[%d] = %"FFMT"\n",istart,rho[istart]);
	        (void) printf("\n");
	        clean_up(ERROR);
	    }
	    else
		istart = 0;
	}
	for (iend = vsize-1; iend >= 0; --iend)
	    if (rho[iend] <= Rho_max)
		break;
	++iend;
	if (iend == 0)
	{
	    if (MG_OOR)
	    {
	        screen("ERROR in MG_single_eos_load_pressure_and_gammas(), "
		       "all densities out of EOS range\n");
	        (void) printf("vsize = %d\n",vsize);
	        (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	        for (iend = vsize-1; iend >= 0; --iend)
	            (void) printf("    rho[%d] = %"FFMT"\n",iend,rho[iend]);
	        (void) printf("\n");
	        clean_up(ERROR);
	    }
	    else
		iend = vsize;
	}
	for (i = istart; i < iend; ++i)
	{
    	    if (rho[i] > Rho_max)
	    {                       
		if ((istart < i) && ((i+1) < iend))
		{
		    if ((rho[i-1] < Rho_max) && (rho[i+1] < Rho_max))
		    {
			rho[i] = 0.5*(rho[i-1] + rho[i+1]);
		    }
		    else
			rho[i] = min(rho[i-1],rho[i+1]);
		}
		else if (MG_OOR)
		{
		    int j;
	            screen("ERROR in "
		           "MG_single_eos_load_pressure_and_gammas(), "
		           "interior region state out of range\n");
	            (void) printf("vsize = %d\n",vsize);
	            (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	            (void) printf("rho[%d] = %"FFMT"\n",i,rho[i]);
		    (void) printf("istart = %d, iend = %d\n",istart,iend);
	            for (j = 0; j < vsize; ++j)
	                (void) printf("    rho[%d] = %"FFMT"\n",j,rho[j]);
	            (void) printf("\n");
		    clean_up(ERROR);
		}
	    }
	    pr = P_ref(rho[i],mg);
	    dpr = dP_ref(rho[i],mg);
	    er = E_ref(rho[i],mg);
	    der = dE_ref(rho[i],mg);
	    Gm = Gamma(rho[i],mg);
	    dGm = dGamma(rho[i],mg);
	    p[i] = pr + rho[i]*Gm*(e[i] - er);
	    c2r = (Gm/rho[i])*(pr+der) - dpr/(rho[i]*rho[i]);
	    c2[i] = (e[i] - er)*(Gm*(Gm+1.0) - dGm/rho[i]) + c2r;
	    if (c2[i] < 0.0)
		c2[i] = -c2[i];
	    GAM[i] = Gm;
	    if (FD != NULL)
	    {
	        if (c2r > 0.0)
	        {
	            d2pr = d2P_ref(rho[i],mg);
	            d2er = d2E_ref(rho[i],mg);
	            d2Gm = d2Gamma(rho[i],mg);
	            k = pr + der;
	            dk = dpr + d2er;
	            Gr = 1.0 + 0.5*Gm +
			(0.5/(rho[i]*rho[i]*c2r))*
			(-2.0*dGm*k-Gm*dk+(Gm+2.0)*dpr + d2pr/rho[i]);

	            FD[i] = Gr + ((e[i]-er)/(rho[i]*rho[i]*2.0*c2[i]*c2r))*(
		        k*(2.0*Gm*dGm-2.0*dGm*dGm/rho[i] + Gm*d2Gm/rho[i]) +
		        dk*(Gm*Gm*(Gm + 1.0) - Gm*dGm/rho[i]) -
		        dpr*(Gm*(2.0 + Gm*(3.0 + Gm) - 3.0*dGm/rho[i]) -
		        2.0*dGm/rho[i] + d2Gm/(rho[i]*rho[i])) -
		        d2pr*((Gm*(Gm+1.0) - dGm/rho[i])/rho[i]));
	        }
	        else
	            FD[i] = 0.5*(Gm + 2.0);
	    }
	}
	/*
	*   When the states at the end of the uni_arrays are out of range
	*   it is usually due to irregular stencil points that will be
	*   overwritten with their correct values in the irregular interior
	*   state sweep.  Copy out the adjacent interior values to these
	*   temporary values.
	*/
	for (i = 0; i < istart; ++i)
	{
	    p[i] = p[istart];
	    c2[i] = c2[istart];
	    GAM[i] = GAM[istart];
	    if (FD != NULL)
		FD[i] = FD[istart];
	}
	for (i = iend; i < vsize; ++i)
	{
	    p[i] = p[iend-1];
	    c2[i] = c2[iend-1];
	    GAM[i] = GAM[iend-1];
	    if (FD != NULL)
		FD[i] = FD[iend-1];
	}
}		/*end MG_single_eos_load_pressure_and_gammas*/

/*
*			MG_single_eos_load_pressure():
*
*	Loads a vector of pressures into the appropriate field of the 
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

/*ARGSUSED*/
LOCAL	void	MG_single_eos_load_pressure(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	double     *rho = vst->rho + offset;
	double     *p = vst->p + offset;
	double     *e = vst->e + offset;
	double     pr, er, Gm;
	double     Rho_max;
	int       i, istart, iend;
    	MG_params *mg;

	mg = &MG_Eos(vst->state[offset])->MGparams;
	Rho_max = mg->Rho_max;
	for (istart = 0; istart < vsize; ++istart)
	    if (rho[istart] <= Rho_max)
		break;
	if (istart == vsize)
	{
	    if (MG_OOR)
	    {
	        screen("ERROR in MG_single_eos_load_pressure(), "
		       "all densities out of EOS range\n");
	        (void) printf("vsize = %d\n",vsize);
	        (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	        for (istart = 0; istart < vsize; ++istart)
	            (void) printf("    rho[%d] = %"FFMT"\n",istart,rho[istart]);
	        (void) printf("\n");
	        clean_up(ERROR);
	    }
	    else
		istart = 0;
	}
	for (iend = vsize-1; iend >= 0; --iend)
	    if (rho[iend] <= Rho_max)
		break;
	++iend;
	if (iend == 0)
	{
	    if (MG_OOR)
	    {
	        screen("ERROR in MG_single_eos_load_pressure(), "
		       "all densities out of EOS range\n");
	        (void) printf("vsize = %d\n",vsize);
	        (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	        for (iend = vsize-1; iend >= 0; --iend)
	            (void) printf("    rho[%d] = %"FFMT"\n",iend,rho[iend]);
	        (void) printf("\n");
	        clean_up(ERROR);
	    }
	    else
		iend = vsize;
	}
	for (i = istart; i < iend; ++i)
	{
    	    if (rho[i] > Rho_max)
	    {                       
		if ((istart < i) && ((i+1) < iend))
		{
		    if ((rho[i-1] < Rho_max) && (rho[i+1] < Rho_max))
		    {
			rho[i] = 0.5*(rho[i-1] + rho[i+1]);
		    }
		    else
			rho[i] = min(rho[i-1],rho[i+1]);
		}
		else if (MG_OOR)
		{
		    int j;
	            screen("ERROR in "
		           "MG_single_eos_load_pressure(), "
		           "interior region state out of range\n");
	            (void) printf("vsize = %d\n",vsize);
	            (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	            (void) printf("rho[%d] = %"FFMT"\n",i,rho[i]);
		    (void) printf("istart = %d, iend = %d\n",istart,iend);
	            for (j = 0; j < vsize; ++j)
	                (void) printf("    rho[%d] = %"FFMT"\n",j,rho[j]);
	            (void) printf("\n");
		    clean_up(ERROR);
		}
	    }
	    pr = P_ref(rho[i],mg);
	    er = E_ref(rho[i],mg);
	    Gm = Gamma(rho[i],mg);
	    p[i] = pr + rho[i]*Gm*(e[i] - er);
	}
	/*
	*   When the states at the end of the uni_arrays are out of range
	*   it is usually due to irregular stencil points that will be
	*   overwritten with their correct values in the irregular interior
	*   state sweep.  Copy out the adjacent interior values to these
	*   temporary values.
	*/
	for (i = 0; i < istart; ++i)
	    p[i] = p[istart];
	for (i = iend; i < vsize; ++i)
	    p[i] = p[iend-1];
}		/*end MG_single_eos_load_pressure*/

/*
*		MG_single_eos_load_sound_speed2():
*
*	Loads a vector of sound speeds into the	appropriate fields of the
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	MG_single_eos_load_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	double     *rho = vst->rho + offset;
	double     *c2 = vst->c2 + offset;
	double     *e = vst->e + offset;
	double     *FD;
	double     pr, dpr, d2pr, er, der, d2er, Gm, dGm, d2Gm;
	double     c2r, Gr;
	double     k, dk;
	double     Rho_max;
	int       i, istart, iend;
    	MG_params *mg;

	FD = (vst->FD != NULL) ? vst->FD + offset : NULL;

	mg = &MG_Eos(vst->state[offset])->MGparams;
	Rho_max = mg->Rho_max;
	for (istart = 0; istart < vsize; ++istart)
	    if (rho[istart] <= Rho_max)
		break;
	if (istart == vsize)
	{
	    if (MG_OOR)
	    {
	        screen("ERROR in "
		       "MG_single_eos_load_sound_speed2(), "
		       "all densities out of EOS range\n");
	        (void) printf("vsize = %d\n",vsize);
	        (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	        for (istart = 0; istart < vsize; ++istart)
	            (void) printf("    rho[%d] = %"FFMT"\n",istart,rho[istart]);
	        (void) printf("\n");
	        clean_up(ERROR);
	    }
	    else
		istart = 0;
	}
	for (iend = vsize-1; iend >= 0; --iend)
	    if (rho[iend] <= Rho_max)
		break;
	++iend;
	if (iend == 0)
	{
	    if (MG_OOR)
	    {
	        screen("ERROR in "
			"MG_single_eos_load_sound_speed2(), "
		       "all densities out of EOS range\n");
	        (void) printf("vsize = %d\n",vsize);
	        (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	        for (iend = vsize-1; iend >= 0; --iend)
	            (void) printf("    rho[%d] = %"FFMT"\n",iend,rho[iend]);
	        (void) printf("\n");
	        clean_up(ERROR);
	    }
	    else
		iend = vsize;
	}
	for (i = istart; i < iend; ++i)
	{
    	    if (rho[i] > Rho_max)
	    {                       
		if ((istart < i) && ((i+1) < iend))
		{
		    if ((rho[i-1] < Rho_max) && (rho[i+1] < Rho_max))
		    {
			rho[i] = 0.5*(rho[i-1] + rho[i+1]);
		    }
		    else
			rho[i] = min(rho[i-1],rho[i+1]);
		}
		else if (MG_OOR)
		{
		    int j;
	            screen("ERROR in "
		           "MG_single_eos_load_sound_speed2(), "
		           "interior region state out of range\n");
	            (void) printf("vsize = %d\n",vsize);
	            (void) printf("Rho_max = %"FFMT"\n",Rho_max);
	            (void) printf("rho[%d] = %"FFMT"\n",i,rho[i]);
		    (void) printf("istart = %d, iend = %d\n",istart,iend);
	            for (j = 0; j < vsize; ++j)
	                (void) printf("    rho[%d] = %"FFMT"\n",j,rho[j]);
	            (void) printf("\n");
		    clean_up(ERROR);
		}
	    }
	    pr = P_ref(rho[i],mg);
	    dpr = dP_ref(rho[i],mg);
	    er = E_ref(rho[i],mg);
	    der = dE_ref(rho[i],mg);
	    Gm = Gamma(rho[i],mg);
	    dGm = dGamma(rho[i],mg);
	    c2r = (Gm/rho[i])*(pr+der) - dpr/(rho[i]*rho[i]);
	    c2[i] = (e[i] - er)*(Gm*(Gm+1.0) - dGm/rho[i]) + c2r;
	    if (c2[i] < 0.0)
		c2[i] = -c2[i];
	    if (FD != NULL)
	    {
	        if (c2r > 0.0)
	        {
	            d2pr = d2P_ref(rho[i],mg);
	            d2er = d2E_ref(rho[i],mg);
	            d2Gm = d2Gamma(rho[i],mg);
	            k = pr + der;
	            dk = dpr + d2er;
	            Gr = 1.0 + 0.5*Gm +
			(0.5/(rho[i]*rho[i]*c2r))*
			(-2.0*dGm*k-Gm*dk+(Gm+2.0)*dpr + d2pr/rho[i]);

	            FD[i] = Gr + ((e[i]-er)/(rho[i]*rho[i]*2.0*c2[i]*c2r))*(
		        k*(2.0*Gm*dGm-2.0*dGm*dGm/rho[i] + Gm*d2Gm/rho[i]) +
		        dk*(Gm*Gm*(Gm + 1.0) - Gm*dGm/rho[i]) -
		        dpr*(Gm*(2.0 + Gm*(3.0 + Gm) - 3.0*dGm/rho[i]) -
		        2.0*dGm/rho[i] + d2Gm/(rho[i]*rho[i])) -
		        d2pr*((Gm*(Gm+1.0) - dGm/rho[i])/rho[i]));
	        }
	        else
	            FD[i] = 0.5*(Gm + 2.0);
	    }
	}
	/*
	*   When the states at the end of the uni_arrays are out of range
	*   it is usually due to irregular stencil points that will be
	*   overwritten with their correct values in the irregular interior
	*   state sweep.  Copy out the adjacent interior values to these
	*   temporary values.
	*/
	for (i = 0; i < istart; ++i)
	{
	    c2[i] = c2[istart];
	    if (FD != NULL)
		FD[i] = FD[istart];
	}
	for (i = iend; i < vsize; ++i)
	{
	    c2[i] = c2[iend-1];
	    if (FD != NULL)
		FD[i] = FD[iend-1];
	}
}		/*end MG_single_eos_load_sound_speed2*/

/***************END VECTORIZED THERMODYNAMIC FUNCTIONS *********************/

/***************RIEMANN SOLUTIONS UTILITY FUNCTIONS ************************/
/***************Purely Thermodynamic Hugoniot Functions*********************/

/*
*			MG_dens_Hugoniot():
*
*	Given the state state0 on one side of an oblique shock and the pressure
*	p1 on the other side, this function returns the density rho1 of the
*	state with pressure p1.  Rho1 is found by solving the Hugoniot relation
*
*		(p1 + p0)*(1/rho0 - 1/rho1) = 2*(e1 - e0)
*
*	where e0 and e1 are the specific internal energies of the two
*	respective states.  For a given equation of state the specific
*	internal energy can be expressed as a function of the
*	pressure and density.  Thus the above equation can be solved to
*	give rho1 as a function of state0 and p1.
*
*
*	Reference: Courant and Friedrichs page 302 ff.
*/


LOCAL	double	MG_dens_Hugoniot(
	double    p2,
	Locstate st1)
{
    	double       rho2, p, pmin, dp;
    	double       rho1, e1, p1, rho_last;
	double       Gm1, dGm1, Gm2, dGm2;
	double       pr2, dpr2, er2, der2;
	double       pr1, dpr1, der1;
	double       rho_l, rho_u;
	double       p_l, p_u, c2;
	double       D;
	double       rho_start, p_start;
	const double eps = MACH_EPS; /*TOLERANCE*/
    	MG_params   *mg;
	const int   MAX_ITER = 60;
	int         i, j;

	if (p2 == -HUGE_VAL)
	    return 0.0;

	mg = &MG_Eos(st1)->MGparams;
	rho1 = Dens(st1);
	test_out_of_range(rho1,st1,"MG_dens_Hugoniot, rho1",YES);
	p1 = MG_pressure(st1);
	if (fabs(p2/p1 - 1.0) < eps)
	    return rho1;
	e1 = MG_specific_internal_energy(st1);
	Gm1 = Gamma(rho1,mg);
	dGm1 = dGamma(rho1,mg);
	pr1 = P_ref(rho1,mg);
	dpr1 = dP_ref(rho1,mg);
	der1 = dE_ref(rho1,mg);

	c2 = Gm1*(pr1+der1)/rho1 - dpr1/(rho1*rho1) +
			     (Gm1 + 1.0 - dGm1/(Gm1*rho1))*(p1 - pr1)/rho1;
	if ((c2 <= 0.0) && (p2 < p1)) /* st1 is a vacuum */
	    return rho1;

	pmin = Min_pressure(st1);
	if (p2 < pmin)
	{
	    if (p1 <= pmin) /* st1 is a vacuum */
		return rho1;
	    p2 = pmin;
	}
	if (p2 < p1) /*Rarefaction shock*/
	{
	    double fac = 0.5;
	    rho_u = rho_start = rho1;
	    p_u = p1;
	    rho_l = rho1;
	    p_l = p_start = p1;
	    for (i = 0, j = 0; i < MAX_ITER; ++i)
	    {
		rho_l *= fac;
		if ((rho_l < eps) && ((1.0 - fac) < eps))
		{
		    rho2 = Params(st1)->vacuum_dens;
	            return rho2;
		}
		if (rho_l < eps)
		{
		    if (j == MAX_ITER)
		    {
			i = MAX_ITER;
			break;
		    }
		    ++j;
		    i = -1;
		    rho_l = rho_start;
		    p_l = p_start;
		    fac = 1.0 - 0.5*(1.0 - fac);
		}
		else
		{
	            Gm2 = Gamma(rho_l,mg);
	            D = rho1 - 0.5*Gm2*(rho_l - rho1);
	            if (D <= 0.0)
		    {
		        if (j == MAX_ITER)
		        {
			    i = MAX_ITER;
			    break;
		        }
		        ++j;
		        i = -1;
		        rho_l = rho_start;
		        p_l = p_start;
		        fac = 1.0 - 0.5*(1.0 - fac);
		    }
		    else
		    {
	                pr2 = P_ref(rho_l,mg);
	                dpr2 = dP_ref(rho_l,mg);
	                er2 = E_ref(rho_l,mg);
	                der2 = dE_ref(rho_l,mg);
	                dGm2 = dGamma(rho_l,mg);
	                p = pr2 + Gm2*(rho_l*rho1*(e1-er2) +
				    0.5*(p1+pr2)*(rho_l-rho1))/D;
	                c2 = Gm2*(pr2+der2)/rho_l - dpr2/(rho_l*rho_l) +
			     (Gm2 + 1.0 - dGm2/(Gm2*rho_l))*(p - pr2)/rho_l;
		        if ((c2 <= 0.0) || (p1 < p))
		        {
		            if (j == MAX_ITER)
		            {
			        i = MAX_ITER;
			        break;
		            }
		            ++j;
		            i = -1;
		            rho_l = rho_start;
		            p_l = p_start;
		            fac = 1.0 - 0.5*(1.0 - fac);
		        }
		        else
			{
		            p_l = p_start = p;
			    rho_start = rho_l;
			}
		    }
		}
	        if (p_l <= p2)
		    break;
		else
		{
		    p_u = p_l;
		    rho_u = rho_l;
		}
	    }
	    if (i == MAX_ITER)
	    {
		(void) printf("WARNING in MG_dens_Hugoniot(), "
		              "can't find lower limit for density\n");
		(void) printf("rho1 = %"FFMT", rho_l = %"FFMT", "
			      "rho1 - rho_l = %"FFMT", fac = %"FFMT"\n",
			      rho1,rho_l,rho1-rho_l,fac);
		(void) printf("p_l = %"FFMT", p_u = %"FFMT"\n",p_l,p_u);
		(void) printf("rho_l = %"FFMT", rho_u = %"FFMT", "
			      "rho1 = %"FFMT"\n",rho_l,rho_u,rho1);
		(void) printf("p1 = %"FFMT", p2 = %"FFMT", c2 = %"FFMT"\n",
			      p1,p2,c2);
		verbose_print_state("st1",st1);
		(void) printf("raw gas data\n");
		fprint_raw_gas_data(stdout,st1,Params(st1)->dim);
		return rho_l;
	    }
	}
	else if (p1 < p2) /*Shock*/
	{
	    double fac = 2.0;
	    double Rho_max = mg->Rho_max;

	    rho_l = rho1;
	    p_l = p1;
	    rho_u = rho_start = rho1;
	    p_u = p_start = p1;
	    for (i = 0, j = 0; i < MAX_ITER; ++i)
	    {
		if ((fac - 1.0) < eps)
	            return rho1;
		rho_u *= fac;
		if (rho_u >= Rho_max)
		    rho_u = 0.5*(rho_l + Rho_max);
		if (Rho_max - rho_u < eps)
	            return rho1;
	        Gm2 = Gamma(rho_u,mg);
	        D = rho1 - 0.5*Gm2*(rho_u - rho1);
	        if (D <= 0.0)
		{
		    if (j == MAX_ITER)
		    {
			i = MAX_ITER;
			break;
		    }
		    ++j;
		    i = -1;
		    rho_u = rho_start;
		    p_u = p_start;
		    fac = 1.0 + 0.5*(fac - 1.0);
		}
		else
		{
	            pr2 = P_ref(rho_u,mg);
	            dpr2 = dP_ref(rho_u,mg);
	            er2 = E_ref(rho_u,mg);
	            der2 = dE_ref(rho_u,mg);
	            dGm2 = dGamma(rho_u,mg);
	            p = pr2 + Gm2*(rho_u*rho1*(e1-er2) +
				    0.5*(p1+pr2)*(rho_u-rho1))/D;
	            c2 = Gm2*(pr2+der2)/rho_u - dpr2/(rho_u*rho_u) +
			 (Gm2 + 1.0 - dGm2/(Gm2*rho_u))*(p - pr2)/rho_u;
		    if ((c2 <= 0.0) || (p < p1))
		    {
		        if (j == MAX_ITER)
		        {
			    i = MAX_ITER;
			    break;
		        }
		        ++j;
		        i = -1;
		        rho_u = rho_start;
		        p_u = p_start;
		        fac = 1.0 + 0.5*(fac - 1.0);
		    }
		    else
		    {
	                p_u = p_start = p;
			rho_start = rho_u;
		    }
		}
		if (p2 <= p_u)
		    break;
		else
		{
		    p_l = p_u;
		    rho_l = rho_u;
		}
	    }
	    if (i == MAX_ITER)
	    {
		(void) printf("WARNING in MG_dens_Hugoniot(), "
		              "can't find upper limit for density\n");
		(void) printf("rho1 = %"FFMT", rho_u = %"FFMT", "
			      "rho_u - rho1 = %"FFMT", fac = %"FFMT"\n",
			      rho1,rho_u,rho_u-rho1,fac);
		(void) printf("p_l = %"FFMT", p_u = %"FFMT"\n",p_l,p_u);
		(void) printf("rho_l = %"FFMT", rho_u = %"FFMT"\n",rho_l,rho_u);
		(void) printf("p2 = %"FFMT", c2 = %"FFMT"\n",p2,c2);
		verbose_print_state("st1",st1);
		(void) printf("raw gas data\n");
		fprint_raw_gas_data(stdout,st1,Params(st1)->dim);
		return rho_u;
	    }
	}
	else
	    return rho2;
	rho2 = rho_u;
	for (i = 0; i < MAX_ITER; ++i)
	{
	    pr2 = P_ref(rho2,mg);
	    er2 = E_ref(rho2,mg);
	    dpr2 = dP_ref(rho2,mg);
	    der2 = dE_ref(rho2,mg);
	    Gm2 = Gamma(rho2,mg);
	    dGm2 = dGamma(rho2,mg);
	    D = rho1 - 0.5*Gm2*(rho2 - rho1);
	    if (D <= 0.0)
	    {
		screen("ERROR in MG_dens_Hugoniot(), unexpected case "
		       "interior point of full Hugoniot compression\n");
		clean_up(ERROR);
	    }
	    p = pr2 + Gm2*(rho2*rho1*(e1-er2) + 0.5*(p1+pr2)*(rho2-rho1))/D;
	    if (p < p2)
	    {
		rho_l = rho2;
		p_l = p;
	    }
	    else if (p2 < p)
	    {
	        rho_u = rho2;
		p_u = p;
	    }
	    else
		break;

	    if (fabs((p - p2)/p1) < eps)
		break;

	    dp = (rho1/D)*( -Gm2*rho2*der2 - ((Gm2*rho2 - dGm2)/Gm2)*(p - pr2)
		    + dpr2 -0.5*Gm2*rho2*(p + p1) );

	    rho_last = rho2;
	    rho2 = rho2*dp/(dp - (p - p2)*rho2);
	    if ((rho2 <= 0.0) || (rho_u < rho2) || (rho2 < rho_l))
		rho2 = 0.5*(rho_u + rho_l);
	    test_out_of_range(rho2,st1,"MG_dens_Hugoniot, rho2",YES);
	    if (fabs((rho2 - rho_last)/rho1) < eps)
		break;
	}
	return rho2;
}		/*end MG_dens_Hugoniot*/

/***************Functions to Compute Riemann Solutions**********************/


/*
*			MG_riemann_wave_curve():
*
*	Evalutes the forward wave family wave curve defined by
*
*		 _
*		|
*		|
*		|                                1/2
*               |   [ (Pstar  -  P1) * ( V1 - V) ]     if Pstar > P1
*		|
*		|
*	        / 
*	       /
*              \
*		\		
*		|
*               |        / Pstar     |
*               |       /            |
*               |       \      dP    |
*               |        \   ------  |		       if Pstar < P1
*               |         \   rho c  |
*               |         /          |
*               |        / P1        | S
*               |_
*
*/

LOCAL	double	MG_riemann_wave_curve(
	Locstate state1,
	double    pstar)
{
	double p1;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double pmin;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	p1 = pressure(state1);

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pmin = Min_pressure(state1);
	if (pstar < pmin)
	    pstar = pmin;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	if (pstar > p1)
	{
	    double V1, Vstar;
	    V1 = 1.0/Dens(state1);
	    Vstar = 1.0/dens_Hugoniot(pstar,state1);
	    return sqrt((pstar-p1)*(V1 - Vstar));
	}
	else
	    return (pstar - p1)*pr_Riemann_function(pstar,state1);
}		/*end MG_riemann_wave_curve*/

/*
*			pr_Riemann_function():
*
*	Evaluates the quantity
*
*
*			          / p           |
*		             1    \             |
*		RF(p) = ---------- \ dP/(rho*c) |
*		          p - p1   /            |
*			          /p1           | S = S1
*
*/

LOCAL	double	pr_Riemann_function(
	double    p,
	Locstate state1)
{
	double       rho1, p1, rho, c1, T1, c, c2;
	double       RF;
	const double rtol = 1.0e-6; /*TOLERANCE*/
	MG_params   *mg;

	mg = &MG_Eos(state1)->MGparams;
	rho1 = Dens(state1);
	p1 = pressure(state1);
	c1 = sound_speed(state1);
	if (p == p1)
	    return 1.0/(rho1*c1);
	T1 = temperature(state1);

	rho = MG_dens_rarefaction(p,state1);
	c2 = MG_sound_speed_squared_on_adiabat(rho,rho1,T1,mg);
	c = sqrt(c2);
	if ( (fabs(  c -   c1) < rtol*  c1) ||
	     (fabs(  p -   p1) < rtol*  p1) ||
	     (fabs(rho - rho1) < rtol*rho1) )
	    RF = 0.5*( 1.0/(rho1*c1) + 1.0/(rho*c) );
	else
	    RF = int_c_over_rho_drho(rho,state1)/(p - p1);
	return RF;
}		/*end pr_Riemann_function*/

/*
*			int_c_over_rho_drho():
*
*	Evaluates the quantity
*
*
*		          / rho2        |
*	                  \     c       |
*		           \   --- drho |
*		           /   rho      |
*		          / rho1        | S = S1
*
*/

LOCAL   double  dRF(double,POINTER);

struct _dRFprms {
    double     T1, rho1;
    MG_params *mg;
};

LOCAL	double	int_c_over_rho_drho(
	double    rho2,
	Locstate state1)
{
	struct _dRFprms   RF;
	double             rho1;
	double             epsabs, epsrel;
	double             I;
	double             abserr;
	int               neval;
	QUADRATURE_STATUS ier;

	RF.mg = &MG_Eos(state1)->MGparams;
	RF.rho1 = rho1 = Dens(state1);
	RF.T1 = MG_temperature(state1);

	epsrel = 1.0e-6;/*TOLERANCE*/
	epsabs = RF.mg->C0*epsrel;
	I = dqng(dRF,(POINTER)&RF,rho1,rho2,epsabs,epsrel,&abserr,&neval,&ier);
	switch (ier)
	{
	case INVALID_EPSILON:
	    screen("ERROR in int_c_over_rho_drho(), invalid epsilons\n"
		   "epsabs = %"FFMT", epsrel = %"FFMT"\n",epsabs,epsrel);
	    clean_up(ERROR);
	    break;
	case INACCURATE_INTEGRAL:
	    if ((fabs(abserr) > 20.0*epsabs) && (fabs(abserr) > 20.0*epsrel*I))
	    {
		/*
		 * Don't print a warning is we are close to satifying the
		 * error requirements
		 */
	        (void) printf("WARNING in int_c_over_rho_drho(), "
			      "inaccurate result\n \tneval = %d, "
			      "abserr = %"FFMT" result = %"FFMT"\n",
			      neval,abserr,I);
	    }
	    break;
	case ACCURATE_INTEGRAL:
	default:
	    break;
	}
	return I;
}		/*end int_c_over_rho_drho*/

/*ARGSUSED*/
LOCAL	double dRF(
	double	rho,
	POINTER prms)
{
	struct _dRFprms *RF = (struct _dRFprms   *)prms;
	double           c2;

	if (rho <= 0.0)
	    return 0.0;
	c2 = MG_sound_speed_squared_on_adiabat(rho,RF->rho1,RF->T1,RF->mg);
	if (c2 <= 0.0)
	{
	    (void) printf("WARNING in dRF(), imaginary sound speed\n");
	    return 0.0;
	}
	else
	    return sqrt(c2)/rho;
}		/*end dRF*/

/***************Purely Thermodynamic Hugoniot Functions*********************/


/***************End Purely Thermodynamic Hugoniot Functions*****************/
/***************Purely Thermodynamic Adiabatic Wave Curve Functions*********/

/*	
*			MG_dens_rarefaction():
*
*	Given the state st1 and the pressure p2 on the other side of
*	a simple wave in steady irrotational flow, this
* 	function returns the density on the other side.
*
*	The answer is give by the solution of the ordinary differential
*	equation
*
*		dh/dP = V,  h(p1) = h1;
*
*	where h is the specific enthalpy,  and the derivatives are taken
*	at constant entropy.
*/

LOCAL	double	MG_dens_rarefaction(
	double    p2,
	Locstate st1)
{
    	double rho2, p, pmin, dp;
    	double rho1, p1, T1, Cv, rho_last;
	double Gam2, dGam2;
	double pr, dpr, tr, der;
	double rho_l, rho_u;
	double p_l, p_u, c2;
	double rho_start, p_start;
	const double eps = MACH_EPS; /*TOLERANCE*/
    	MG_params *mg;
	const int MAX_ITER = 60;
	int i, j;

	if (p2 == -HUGE_VAL)
	    return 0.0;

	mg = &MG_Eos(st1)->MGparams;
	rho1 = Dens(st1);
	test_out_of_range(rho1,st1,"MG_dens_rarefaction, rho1",YES);
	p1 = MG_pressure(st1);
	if (p1 == p2)
	    return rho1;
	T1 = MG_temperature(st1);
	Cv = mg->CV;
	c2 = MG_sound_speed_squared_on_adiabat(rho1,rho1,T1,mg);
	if ((c2 <= 0.0) && (p2 < p1)) /* st1 is a vacuum */
	    return rho1;
	pmin = Min_pressure(st1);
	if (p2 < pmin)
	{
	    if (p1 <= pmin) /* st1 is a vacuum */
		return rho1;
	    p2 = pmin;
	}
	if (p2 < p1) /*Rarefaction*/
	{
	    double fac = 0.5;
	    rho_u = rho1;
	    p_u = p1;
	    rho_l = rho_start = rho1;
	    p_l = p_start = p1;
	    for (i = 0, j = 0; i < MAX_ITER; ++i)
	    {
		rho_l *= fac;
		if ((rho_l < eps) && ((1.0 - fac) < eps))
		{
		    rho2 = Params(st1)->vacuum_dens;
	            return MG_real_sound_speed_dens(rho2,rho1,T1,mg);
		}
		if (rho_l < eps)
		{
		    if (j == MAX_ITER)
		    {
			i = MAX_ITER;
			break;
		    }
		    ++j;
		    i = -1;
		    rho_l = rho_start;
		    p_l = p_start;
		    fac = 1.0 - 0.5*(1.0 - fac);
		}
		else
		{
	            c2 = MG_sound_speed_squared_on_adiabat(rho_l,rho1,T1,mg);
	            p = MG_pressure_on_adiabat(rho_l,rho1,T1,mg);
		    if ((c2 <= 0.0) || (p1 < p))
		    {
		        if (j == MAX_ITER)
		        {
			    i = MAX_ITER;
			    break;
		        }
		        ++j;
		        i = -1;
		        rho_l = rho_start;
		        p_l = p_start;
		        fac = 1.0 - 0.5*(1.0 - fac);
		    }
		    else
		    {
	                p_l = p;
			p_start = p_l;
			rho_start = rho_l;
		    }
		}
		if (p_l <= p2)
		    break;
		else
		{
		    p_u = p_l;
		    rho_u = rho_l;
		}
	    }
	    if (i == MAX_ITER)
	    {
		(void) printf("WARNING in MG_dens_rarefaction(), "
		              "can't find lower limit for density\n");
		(void) printf("rho1 = %"FFMT", rho_l = %"FFMT", "
			      "rho1 - rho_l = %"FFMT", fac = %"FFMT"\n",
			      rho1,rho_l,rho1-rho_l,fac);
		(void) printf("p_l = %"FFMT", p_u = %"FFMT"\n",p_l,p_u);
		(void) printf("rho_l = %"FFMT", rho_u = %"FFMT"\n",rho_l,rho_u);
		(void) printf("p2 = %"FFMT", c2 = %"FFMT"\n",p2,c2);
		verbose_print_state("st1",st1);
		(void) printf("raw gas data\n");
		fprint_raw_gas_data(stdout,st1,Params(st1)->dim);
		return rho_l;
	    }
	}
	else if (p1 < p2) /*Compression*/
	{
	    double fac = 2.0;
	    double Rho_max = mg->Rho_max;

	    rho_l = rho1;
	    p_l = p1;
	    rho_u = rho_start = rho1;
	    p_u = p_start = p1;
	    for (i = 0, j = 0; i < MAX_ITER; ++i)
	    {
		rho_u *= fac;
		if (rho_u >= Rho_max)
		    rho_u = 0.5*(rho_l + Rho_max);
	        c2 = MG_sound_speed_squared_on_adiabat(rho_u,rho1,T1,mg);
	        p = MG_pressure_on_adiabat(rho_u,rho1,T1,mg);
		if ((c2 <= 0.0))
		{
		    if (j == MAX_ITER)
		    {
			i = MAX_ITER;
			break;
		    }
		    ++j;
		    i = -1;
		    rho_u = rho_start;
		    p_u = p_start;
		    fac = 1.0 + 0.5*(fac - 1.0);
		}
		else
		{
	            p_u = p;
		    p_start = p_u;
		    rho_start = rho_u;
		}
	        if (p2 <= p_u)
		    break;
		else
		{
		    p_l = p_u;
		    rho_l = rho_u;
		}
	    }
	    if (i == MAX_ITER)
	    {
		(void) printf("WARNING in MG_dens_rarefaction(), "
		              "can't find upper limit for density\n");
		(void) printf("rho1 = %"FFMT", rho_u = %"FFMT", "
			      "rho_u - rho1 = %"FFMT", "
			      "fac = %"FFMT"\n",rho1,rho_u,rho_u-rho1,fac);
		(void) printf("p_l = %"FFMT", p_u = %"FFMT"\n",p_l,p_u);
		(void) printf("rho_l = %"FFMT", rho_u = %"FFMT"\n",rho_l,rho_u);
		(void) printf("p2 = %"FFMT", c2 = %"FFMT"\n",p2,c2);
		verbose_print_state("st1",st1);
		(void) printf("raw gas data\n");
		fprint_raw_gas_data(stdout,st1,Params(st1)->dim);
		return rho_u;
	    }
	}
	else
	    return rho1;
	rho2 = rho_u;
	for (i = 0; i < MAX_ITER; ++i)
	{
	    pr = P_ref(rho2,mg);
	    dpr = dP_ref(rho2,mg);
	    der = dE_ref(rho2,mg);
	    tr = T_ref(rho2,mg);
	    Gam2 = Gamma(rho2,mg);
	    dGam2 = dGamma(rho2,mg);
	    p = pr + Cv*Gam2*rho2*(T1*IrhoG(rho2,rho1,mg) - tr);
	    if (p < p2)
	    {
		rho_l = rho2;
		p_l = p;
	    }
	    else if (p2 < p)
	    {
	        rho_u = rho2;
		p_u = p;
	    }
	    else
		break;

	    if (fabs((p - p2)/p1) < eps)
		break;

	    dp = dpr +
		(dGam2/Gam2 - (Gam2+1.0)*rho2)*(p - pr) - Gam2*rho2*(pr + der);

	    rho_last = rho2;
	    rho2 = rho2*dp/(dp - (p - p2)*rho2);
	    if ((rho2 <= 0.0) || (rho_u < rho2) || (rho2 < rho_l))
		rho2 = 0.5*(rho_u + rho_l);
	    test_out_of_range(rho2,st1,"MG_dens_rarefaction, rho2",YES);
	    if (fabs((rho2 - rho_last)/rho1) < eps)
		break;
	}
	return MG_real_sound_speed_dens(rho2,rho1,T1,mg);
}		/*end MG_dens_rarefaction*/

/*
*			MG_real_sound_speed_dens():
*
*	Given a state with density rho on the isentrope through the state with
*	density rho1 and temperature T1,  determine the sound speed at
*	that point,  and if this quantity if imaginary adiabatically raise
*	the density to the point where the sound speed becomes real.
*/

LOCAL	double MG_real_sound_speed_dens(
	double     rho_in,
	double     rho1,
	double     T1,
    	MG_params *mg)
{
	double c2;
	double rho_out, rho_h, rho_l, rho_m;
	double rho_start;
	double fac = 2.0;
	const double eps_l = MACH_EPS; /*TOLERANCE*/
	const double eps_h = 100.0*MACH_EPS; /*TOLERANCE*/
	const int MAX_ITER = 60;
	int i, j;

	c2 = MG_sound_speed_squared_on_adiabat(rho_in,rho1,T1,mg);
	if (c2 >= eps_l)
	    return rho_in;

	c2 = MG_sound_speed_squared_on_adiabat(rho1,rho1,T1,mg);
	if (c2 < eps_l) /*st1 is already a vacuum */
	    return rho1;

	rho_h = rho_start = rho_in;
	for (i = 0, j = 0; i < MAX_ITER; ++i)
	{
	    rho_l = rho_h;
	    rho_h *= fac;
	    rho_out = rho_h;
	    c2 = MG_sound_speed_squared_on_adiabat(rho_out,rho1,T1,mg);
	    if (c2 <= 0.0)
	    {
		if (j == MAX_ITER)
		{
	            i = MAX_ITER;
	            break;
		}
		++j;
		i = -1;
		rho_l = rho_start;
		rho_h = rho_start;
	        rho_out = rho_h;
		fac = 1.0 + 0.5*(fac - 1.0);
	        c2 = MG_sound_speed_squared_on_adiabat(rho_out,rho1,T1,mg);
	    }
	    if (eps_l <= c2)
	        break;
	    rho_start = rho_out;
	}
	if (i == MAX_ITER)
	{
	    (void) printf("WARNING in MG_real_sound_speed_dens(), "
		          "can't sound speed >= eps_l = %"FFMT", "
			  "c2 = %"FFMT"\n",eps_l,c2);
	    return rho_out;
	}
	if (c2 > eps_h)
	{
	    for (i = 0; i < MAX_ITER; ++i)
	    {
		rho_m = 0.5*(rho_l + rho_h);
	        rho_out = rho_m;
	        c2 = MG_sound_speed_squared_on_adiabat(rho_out,rho1,T1,mg);
		if (c2 < eps_l)
		    rho_l = rho_m;
		if (c2 > eps_h)
		    rho_h = rho_m;
	        if ( (eps_l <= c2) && (c2 <= eps_h) )
		    break;
	    }
	    if (i == MAX_ITER)
	    {
	        (void) printf("WARNING in MG_real_sound_speed_dens(), "
		       "can't sound speed betwen eps_l = %"FFMT" and "
		       "ups_h = %"FFMT", c2 = %"FFMT"\n",eps_l,eps_h,c2);
		return rho_out;
	    }
	}
	return rho_out;
}

/*	
*			MG_pressure_rarefaction():
*
*	Given the state st1 and the density on the other side of
*	a simple wave in steady irrotational flow, this
* 	function returns the pressure on the other side.
*
*	The answer is give by the solution of the ordinary differential
*	equation
*
*		de/dV = -P,  e(V1) = e1;
*
*	where e is the specific internal energy,  and the derivatives are taken
*	at constant entropy.
*/

LOCAL	double	MG_pressure_rarefaction(
	double    rho2,
	Locstate st1)
{
    	double p2;
    	double rho1, T1;
    	MG_params *mg;

	mg = &MG_Eos(st1)->MGparams;
	rho1 = Dens(st1);
	test_out_of_range(rho1,st1,"MG_pressure_rarefaction, rho1",YES);
	if (test_out_of_range(rho2,st1,"MG_pressure_rarefaction, rho2",NO))
	{
	    /* Just try to survive */
	    p2 = pressure(st1) + sound_speed_squared(st1)*(rho2 - rho1);
	}
	else
	{
	    T1 = MG_temperature(st1);
	    p2 = MG_pressure_on_adiabat(rho2,rho1,T1,mg);
	}
	return  p2;
}		/*end MG_pressure_rarefaction*/


/***************End Purely Thermodynamic Adiabatic Wave Curve Functions*****/
/***************General Wave Curve Functions********************************/

/*
*			MG_mass_flux():
*
*	Returns the mass flux across a wave.
*
*				
*		     | (P - P0) |
*		m  = | -------  |
*		     | (U - U0) |
*
*	Where 
*		P0 = pressure ahead of the shock
*		U0 = velocity ahead of the shock
*		P = pressure behind the shock
*		U = velocity behind the shock
*
*/

LOCAL	double	MG_mass_flux(
	double    p,
	Locstate state0)
{
	return sqrt(MG_mass_flux_squared(p,state0));
}		/*end MG_mass_flux*/

/*
*			MG_mass_flux_squared():
*
*	Returns the square of the mass flux across a wave.
*
*				 2
*		 2   | (P - P0) |
*		m  = | -------  |
*		     | (U - U0) |
*
*	Where 
*		P0 = pressure ahead of the shock
*		U0 = velocity ahead of the shock
*		P = pressure behind the shock
*		U = velocity behind the shock
*
*/

LOCAL	double	MG_mass_flux_squared(
	double p,
	Locstate state0)
{
	double p0 = pressure(state0);
	double m2, V, V0;

	if (fabs(p - p0) < p0*EPSILON)
	    m2 = acoustic_impedance_squared(state0);
	else if (p < p0)
	{
	    double mi;
	    mi = pr_Riemann_function(p,state0);
	    m2 = sqr(1.0/mi);
	}
	else
	{
	    const double epsilon = 10.0*MACH_EPS; /*TOLERANCE*/
	    V0 = 1.0/Dens(state0);
	    V = 1.0/dens_Hugoniot(p,state0);
	    if (fabs(V0-V) < epsilon*V0)
	        m2 = acoustic_impedance_squared(state0);
	    else
	        m2 = (p-p0)/(V0-V);
	    if (m2 < 0.0)
	    {
	        screen("ERROR in MG_mass_flux_squared(), "
			"imaginary mass flux\n");
		(void) printf("m2 = %"FFMT", p = %"FFMT", "
			      "dens_Hugoniot(p,state0) = %"FFMT"\n",
			      m2,p,dens_Hugoniot(p,state0));
	        verbose_print_state("state0",state0);
		(void) printf("raw gas data\n");
		fprint_raw_gas_data(stdout,state0,Params(state0)->dim);
	        clean_up(ERROR);
	    }
	}
	return m2;
}		/*end MG_mass_flux_squared*/

/***************End General Wave Curve Functions****************************/

/***************Functions for the Evaluation of Riemann Solutions***********/

/*
*			MG_oned_fan_state():
*
*	This is a utility function provided for the evaluation of states
*	in a simple wave.   Given sta, it solves for stm using the
*	equation:
*
*	                     / p_m        |            / p_m        |
*	                    /             |           /             |
*	                    \       dP    |           \       G dP  |
*	    w = c_m - c_a +  \    -----   |         =  \     ------ |
*	                      \   rho c   |             \    rho c  |
*	                      /           |             /           |
*	                     /p_a         | S = S_a    / p_a        | S = S_a
*
*	                                               / c_m        |
*	                                              /             |
*	                                              \        dc   |
*	                                            =  \     ------ |
*	                                                \     mu^2  |
*	                                                /           |
*	                                               / c_a        | S = S_a
*
*
*	here c is the sound speed,  rho the density,  S the specific entropy,
*	p the pressure,  and mu^2 = (G - 1)/G,  where G is the fundamental
*	derivative of gas dynamics.  The returned state1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	state1 can be used to evaluate the pressure, density, and sound speed
*	of the state inside the rarefaction fan.
*	
*	Input data:
*		w = value of w as defined above
*		sta = state ahead of fan
*		stb = state behind fan
*
*	Output data:
*		stm = state inside fan
*		vacuum = 1 if stm is a vacuum,  0 otherwise
*
*	Returns the sound speed of the answer state stm.
*/

struct FAN_AUX	{
	Locstate sta;
	Locstate stm;
	double	ca, cm;
	double   pa;
};
LOCAL	boolean	oned_fan_aux(double,double*,POINTER);

LOCAL	double	MG_oned_fan_state(
	double    w,
	Locstate sta,
	Locstate stb,
	Locstate stm,
	int      stype_m,
	boolean     *vacuum)
{
	double   rhoa, rhob;
	double	rhom, ca;
	double   wb;
	double   eps, delta;
	struct FAN_AUX	Fprms;

	*vacuum = NO;
	Fprms.sta = sta;
	Fprms.stm = stm;
	Fprms.ca = ca = sound_speed(sta);
	eps = EPSILON*ca;
	Fprms.pa = MG_pressure(sta);
	rhoa = Dens(sta);
	rhob = Dens(stb);
	delta = EPSILON*fabs(rhob-rhoa);
	Set_params(stm,sta);

	if (!oned_fan_aux(Dens(stb),&wb,(POINTER)&Fprms))
	{
	    screen("ERROR in MG_oned_fan_state(), can't evalute "
		   "state behind rarefaction fan\n");
	    clean_up(ERROR);
	    return Fprms.cm;
	}
	if (((w <= wb) && (wb <= 0.0)) || ((0.0 <= wb) && (wb <= w)))
	{
	    set_state(stm,stype_m,stb);
            Fprms.cm = sound_speed(stm);
	}
	else if (((wb <= 0.0) && (0.0 <= w)) || ((w <= 0.0) && (0.0 <= wb)))
	{
	    set_state(stm,stype_m,sta);
            Fprms.cm = sound_speed(stm);
	}
	else if (find_root(oned_fan_aux,(POINTER)&Fprms,w,&rhom,rhoa,rhob,
		           eps,delta)) 
	{
	    set_state(stm,stype_m,stm);
	}
	else
	{
	    double wm;
	    *vacuum = YES;
	    screen("ERROR in MG_oned_fan_state(), can't find root\n");
	    (void) printf("w = %"FFMT"\n",w);
	    (void) printf("w(stb) = %"FFMT"\n",wb);
	    verbose_print_state("sta",sta);
	    (void) printf("raw gas data\n");
	    fprint_raw_gas_data(stdout,sta,Params(sta)->dim);
	    verbose_print_state("stb",stb);
	    (void) printf("raw gas data\n");
	    fprint_raw_gas_data(stdout,stb,Params(stb)->dim);
	    (void) printf("Current stm\n");
	    verbose_print_state("stm",stm);
	    (void) printf("raw gas data\n");
	    fprint_raw_gas_data(stdout,stm,Params(stm)->dim);
	    (void) oned_fan_aux(Dens(stm),&wm,(POINTER)&Fprms);
	    (void) printf("w(stm) = %"FFMT",w(stm) - w = %"FFMT"\n",wm,wm - w);
	    (void) printf("Convergence tolerances\n");
	    (void) printf("Solution value = %"FFMT"\n",eps);
	    (void) printf("Solution interval = %"FFMT"\n",delta);
	    clean_up(ERROR);
	    return Fprms.cm;
	}

	return Fprms.cm;
}		/* MG_end oned_fan_state*/

LOCAL	boolean	oned_fan_aux(
	double	rhom,
	double	*w,
	POINTER	prms)
{
	struct FAN_AUX	*fprms = (struct FAN_AUX  *)prms;
	Locstate sta = fprms->sta;
	Locstate stm = fprms->stm;
	double    pm, cm;

	pm = MG_pressure_rarefaction(rhom,sta);
	*w = int_c_over_rho_drho(rhom,sta);
	Dens(stm) = rhom;
	Press(stm) = pm;
	set_type_of_state(stm,TGAS_STATE);
	fprms->cm = cm = sound_speed(stm);
	*w += cm - fprms->ca;
	return FUNCTION_SUCCEEDED;
}		/*end oned_fan_aux*/
/***************End Functions for the Evaluation of Riemann Solutions********/

/***************END RIEMANN SOLUTIONS UTILITY FUNCTIONS ********************/

/***************INITIALIZATION UTILITY FUNCTIONS****************************/

/*
*			MG_fprint_EOS_params():
*
*	Prints the parameters that define the given equation of state.
*	NOTE:  This is not really an initialization function,  but it is
*	convenient to locate it next the the corresponding read function.
*/

LOCAL	void	MG_fprint_EOS_params(
	FILE      *file,
	Gas_param *params)
{
	MG_EOS    *mgeos = (MG_EOS*) params->eos;
	MG_params *mg = &mgeos->MGparams;

	(void) fprintf(file,"\tEquation of state = %d MIE_GRUNEISEN\n",
		       MIE_GRUNEISEN);
	(void) fprintf(file,"Material name = %s\n",mg->name);
	if (is_binary_output())
	{
	    (void) fprintf(file,"\f%c",15);
	    (void) fwrite((const void *) &mg->Rho0,FLOAT,15,file);
	}
	else
	{
	    (void) fprintf(file,"\tRho0, P0, E0, T0 = "
			        "%"FFMT" %"FFMT" %"FFMT" %"FFMT"\n",
			   mg->Rho0,mg->P0,mg->E0,mg->T0);
	    (void) fprintf(file,"\tC0, S1, S2, S3 = "
			        "%"FFMT" %"FFMT" %"FFMT" %"FFMT"\n",
	                   mg->C0,mg->S1,mg->S2,mg->S3);
	    (void) fprintf(file,"\tgamma0, b = %"FFMT" %"FFMT"\n",
			   mg->gamma0,mg->b);
	    (void) fprintf(file,"\tCV = %"FFMT"\n",mg->CV);
	    (void) fprintf(file,"\tRho_max = %"FFMT"\n",mg->Rho_max);
	    (void) fprintf(file,"\tV_min = %"FFMT"\n",mg->V_min);
	    (void) fprintf(file,"\tP_inf = %"FFMT"\n",mg->P_inf);
	    (void) fprintf(file,"\tE_inf = %"FFMT"\n",mg->E_inf);
	}
	(void) fprintf(file,"\n\n");
}		/*end MG_fprint_EOS_params */

LOCAL	void	print_MG_params(
	MG_params *mg)
{
	(void) printf("Material name = %s\n",mg->name);
	(void) printf("\tRho0, P0, E0, T0 = "
		      "%"FFMT" %"FFMT" %"FFMT" %"FFMT"\n",
		      mg->Rho0,mg->P0,mg->E0,mg->T0);
	(void) printf("\tC0, S1, S2, S3 = "
		      "%"FFMT" %"FFMT" %"FFMT" %"FFMT"\n",
	              mg->C0,mg->S1,mg->S2,mg->S3);
	(void) printf("\tgamma0, b = %"FFMT" %"FFMT"\n",mg->gamma0,mg->b);
	(void) printf("\tCV = %"FFMT"\n",mg->CV);
	(void) printf("\tRho_max = %"FFMT"\n",mg->Rho_max);
	(void) printf("\tV_min = %"FFMT"\n",mg->V_min);
	(void) printf("\tP_inf = %"FFMT"\n",mg->P_inf);
	(void) printf("\tE_inf = %"FFMT"\n",mg->E_inf);
}		/*end print_MG_params*/

/*
*			MG_read_print_EOS_params():
*
*	Reads the equation of state data as printed by MG_fprint_EOS_params.
*	This is restart function.
*/

/*ARGSUSED*/
LOCAL	void	MG_read_print_EOS_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
    	FILE   *file = io_type->file;
	MG_EOS *mgeos = (MG_EOS*) params->eos;
	MG_params *mg = &mgeos->MGparams;
	char s[Gets_BUF_SIZE];
	int  c;

	if (!fgetstring(file,"Material name = "))
	{
	    screen("ERROR in MG_read_print_EOS_params(), "
		   "can't find EOS data\n");
	    clean_up(ERROR);
	}
	(void) fgets(s,Gets_BUF_SIZE-2,file);
	s[strlen(s)-1] = '\0';
	mg->name = strdup(s);
	if ((c = getc(file)) != '\f') /*NOBINARY*/
	{
	    (void) ungetc(c,file);
	    mg->Rho0 = fread_float("Rho0, P0, E0, T0 = ",io_type);
	    mg->P0 = fread_float(NULL,io_type);
	    mg->E0 = fread_float(NULL,io_type);
	    mg->T0 = fread_float(NULL,io_type);
	    mg->C0 = fread_float("C0, S1, S2, S3 = ",io_type);
	    mg->S1 = fread_float(NULL,io_type);
	    mg->S2 = fread_float(NULL,io_type);
	    mg->S3 = fread_float(NULL,io_type);
	    mg->gamma0 = fread_float("gamma0, b = ",io_type);
	    mg->b = fread_float(NULL,io_type);
	    mg->CV = fread_float("CV = ",io_type);
	    mg->Rho_max = fread_float("Rho_max = ",io_type);
	    mg->V_min = fread_float("V_min = ",io_type);
	    mg->P_inf = fread_float("P_inf = ",io_type);
	    mg->E_inf = fread_float("E_inf = ",io_type);
	}
	else
	{
	    (void) getc(file);
	    (void) read_binary_real_array(&mg->Rho0,15,io_type);
	}
	initialize_Tref(mg);
}		/*end MG_read_print_EOS_params*/

/*
*			MG_prompt_for_EOS_params():
*
*	Prompts for equation of state parameters.
*/

/*ARGSUSED*/
LOCAL	void	MG_prompt_for_EOS_params(
	INIT_DATA  *init,
	Gas_param  *params,
	const char *message1,
	const char *message2)
{
	boolean found = NO;
	char s[Gets_BUF_SIZE];
	double xi_max;
	MG_EOS *mgeos = (MG_EOS*) params->eos;
	MG_params *mg = &mgeos->MGparams;

	screen("Input the Mie-Gruneisen EOS parameters for the%s gas%s\n",
	       message1,message2);
	screen("Enter the name of a material or custom to input the "
	       "parameters directly\n\t(default = custom): ");
	(void) Gets(s);
	if ((s[0] != '\0') && strcasecmp(s,"custom"))
	{
	    int i;

	    for (i = 0; MGMaterial[i].name != NULL; ++i)
	    {
		if (strcasecmp(s,MGMaterial[i].name)==0)
		{
		    *mg = MGMaterial[i];
		    mg->P0 = 1.0e-6;/*1 bar*/
		    mg->T0 = 300.0;/*Room temperature*/
		    found = YES;
		    break;
		}
	    }
	    if (!found)
		screen("Material %s not found in the database, "
		       "enter the EOS parameters directly\n");
	}
	if (!found)
	{
	    if (s[0] != '\0')
		mg->name = strdup(s);
	    else
		mg->name = strdup("CUSTOM");
	    screen("Enter the reference state\n\t"
	           "density, pressure, and positive temperature: ");
	    (void) Scanf("%f %f %f\n",&mg->Rho0,
			              &mg->P0,
				      &mg->T0);
	    screen("Enter the reference state sound speed: ");
	    (void) Scanf("%f\n",&mg->C0);
	    screen("Enter the Us-Up slope coefficients S1, S2, S3: ");
	    (void) Scanf("%f %f %f\n",&mg->S1,
			              &mg->S2,
				      &mg->S3);
	    screen("Enter the Gruneisen coefficient at the reference state: ");
	    (void) Scanf("%f\n",&mg->gamma0);
	    screen("Enter the Gruneisen coefficient at infinite compression: ");
	    (void) Scanf("%f\n",&mg->b);
	    screen("Enter the specific heat at constant volume: ");
	    (void) Scanf("%f\n",&mg->CV);
	}
	mg->E0 = sqr(mg->C0)/mg->gamma0 - mg->P0/mg->Rho0;
	mg->P_inf = mg->Rho0*sqr(mg->C0)/(1.0 + mg->gamma0) - mg->P0;
	mg->E_inf = (mg->P0 + (mg->gamma0 + 1.0)*mg->P_inf)/
	             (mg->Rho0*mg->gamma0) - mg->E0;
	params->vacuum_dens = max(
		mg->Rho0*pow(MACH_EPS/sqr(mg->C0),1.0/mg->gamma0),
		params->vacuum_dens);
	params->min_pressure = -0.95*mg->P_inf;/*TOLERANCE*/
	params->min_energy = (mg->E_inf > 0.0) ? -HUGE_VAL :
	    (params->min_pressure + (mg->gamma0+1.0)*mg->P_inf)/mg->gamma0;
	if (mg->S1 <= 0.0)
	{
	    screen("ERROR in MG_prompt_for_EOS_params, "
		   "nonpositive S1 = %"FFMT"\n",mg->S1);
	    clean_up(ERROR);
	}
	xi_max = max_compression(mg);
	mg->V_min = (1.0 - xi_max)/mg->Rho0;
	mg->Rho_max = (xi_max < 1.0) ? mg->Rho0/(1.0 - xi_max) : HUGE_VAL;
	initialize_Tref(mg);

	if (debugging("eosref") || debugging("beosref"))
	{
	    char s[80];
	    FILE *file;
	    int  i, m;
	    double rho, V;
	    double Gm, pr, dpr, er, der, Tr, c2r;
	    boolean  binary = NO;
	    binary = debugging("beosref") ? YES : NO;
	    (void) sprintf(s,"%s-EosReferenceCurves",mg->name);
	    if ((file = fopen(s,"w")) == NULL)
		return;
	    foutput(file);
	    (void) fprintf(file,"Equation of State Reference Curves for %s\n",
			   mg->name);
	    (void) fprintf(file,"\n\n");
	    print_machine_parameters(file);
	    (void) fprintf(file,"\n\n");
	    foutput(file);
	    (void) fprintf(file,"%-25s %-26s %-26s %-26s %-26s %-26s\n",
			         "rho","V",  "P",  "e",  "T",  "c2");
	    m = mg->TrefSpline->m;
	    for (i = 0; i < m; ++i)
	    {
		rho = mg->TrefSpline->rho[i];
		V = 1.0/rho;
		Gm = Gamma(rho,mg);
		pr = P_ref(rho,mg);
		dpr = dP_ref(rho,mg);
		er = E_ref(rho,mg);
		der = dE_ref(rho,mg);
		Tr = T_ref(rho,mg);
	        c2r = (Gm/rho)*(pr+der) - dpr/(rho*rho);
		if (binary)
		{
		    (void) fprintf(file,"\f%c",6);
		    (void) fwrite((const void *)&rho,sizeof(double),1,file);
		    (void) fwrite((const void *)&V,sizeof(double),1,file);
		    (void) fwrite((const void *)&pr,sizeof(double),1,file);
		    (void) fwrite((const void *)&er,sizeof(double),1,file);
		    (void) fwrite((const void *)&Tr,sizeof(double),1,file);
		    (void) fwrite((const void *)&c2r,sizeof(double),1,file);
		}
		else
	            (void) fprintf(file,"%-"FFMT" %-"FFMT" %-"FFMT" "
			                "%-"FFMT" %-"FFMT" %-"FFMT"\n",
				        rho,V,pr,er,Tr,c2r);
	    }
	    (void) fprintf(file,"\n\n");
	    trace_foutput(file);
	    (void) fclose(file);
	}
}		/*end MG_prompt_for_EOS_params*/

/*
*			max_compression():
*
*	Computes the maximum value of the compression x = (V0-V)/V0
*	for which the reference curves are finite and the sound speed
*	on the reference curve is real.  The value is determined by finding
*	the minimum root of the two polynomials 
*
*		p1(x) = 1 - (S1 + S2*x + S3*x*x)*x
*
*		p2(x) = 1 + (S1 - 1)*x + (S1*(gamma0-1) + 3*S2)*x*x +
*                        (S1*(b-gamma0) + S2*(2*gamma0 - 3)+5*S3)*x*x*x +
*                        (2*S2*(b-gamma0) + S3*(3*gamma0-5))*x*x*x*x +
*                        (3*S3*(b-gamma0))*x*x*x*x*x +
*
*	on the interval 0 <= xi <= 1.0. Powers of p1(x) appear as the
*	denominator of the reference curves and is constrained to be positive.
*	Thus 0 <= p1(x),  for 0 <= x <= x_max. The numerator of the square
*	of the sound speed along the reference curve is proportional to p2(x)
*	so we require 0 <= p2(x) for 0 <= x <= x_max.
*
*	Note: We assume S1 > 0.
*/

LOCAL	boolean p2_aux(double,double*,POINTER);

LOCAL	double	max_compression(
	MG_params *mg)
{
    	double x_max, c2x_max, rho, c2ref;
	double S1 = mg->S1, S2 = mg->S2, S3 = mg->S3;
	double gamma0 = mg->gamma0, b = mg->b;
	double p2[6];
	double d, x, r, q;
	double p1r, p1rp;
	const double epsilon = 10.0*MACH_EPS; /*TOLERANCE*/
	double xp2min, p2min;
	const int num_iter = 5;
	int i;

	debug_print("eosref","Entered max_compression\n");

	/* Find smallest root of p1(x) in [0, 1] */
	x_max = 1.0;
	if (S2 == 0.0 && S3 == 0.0) /* Common case,  handle directly */
	{
	    x_max = min(x_max,1.0/S1);
	}
	else if (S3 == 0.0)
	{
	    if ((S1*S1 + 4.0*S2) > 0.0)
	    {
		d = sqrt(S1*S1 + 4.0*S2);
		x = (S2 > 0.0) ?  2.0/(d + S1) : -0.5*(S1 - d)/S2;
		x_max = min(x_max,x);
	    }
	}
	else
	{
	    double x1, x2;
	    boolean  do_newton;
	    do_newton = NO;
	    r = 0.0; /*Newton method start point, or smallest root*/
	    if (S2*S2 > 3.0*S1*S3) /*Local max and mins exist*/
	    {
		d = sqrt(S2*S2 - 3.0*S1*S3);
		x1 = -(S2 + d)/(3.0*S3);
		x2 = -(S2 - d)/(3.0*S3);
		if (x2 < x1)
		{
		    x1 = -(S2 - d)/(3.0*S3);
		    x2 = -(S2 + d)/(3.0*S3);
		}
		if ((1.0 < x1) || (x2 < 0.0)) /* No max/min in [0, 1] */
		{
		    q = 1.0 - S1 - S2 - S3; /* p1(1) */
		    if (q < 0.0)
			do_newton = YES;
		    else if (q == 0.0)
			r = 1.0;
		}
		else if ((0.0 <= x1) && (1.0 <= x2)) /* One max/min in [0, 1]*/
		{
		    q = 1.0 - x1*(S1 + x1*(S2 + S3*x1)); /* p1(x1) */
		    if (q < 0.0)
			do_newton = YES;
		    else if (q == 0.0)
			r = x1;
		}
		else if ((x1 <= 0.0) && (x2 <= 1.0)) /* One max/min in [0, 1]*/
		{
		    q = 1.0 - x2*(S1 + x2*(S2 + S3*x2)); /* p1(x2) */
		    if (q < 0.0)
			do_newton = YES;
		    else if (q == 0.0)
			r = x2;
		}
		else /* Two max/min in [0,1] */
		{
		    q = 1.0 - x1*(S1 + x1*(S2 + S3*x1)); /* p1(x1) */
		    if (q < 0.0)
			do_newton = YES;
		    else if (q == 0.0)
			r = x1;
		    else
		    {
		        q = 1.0 - S1 - S2 - S3; /* p1(1) */
			r = 1.0;
		        if (q < 0.0)
			{
			    do_newton = YES;
		            for (i = 0; i < num_iter; ++i)
		            {
		                p1r = 1.0 - r*(S1 + r*(S2 + S3*r));
		                p1rp = -(S1 + r*(2.0*S2 + 3.0*r*S3));
		                r -= p1r/p1rp;
		            }
		            q = 1.0 - r*(S1 + r*(S2 + S3*r)); /* p1(r) */
			    if (q < 0.0)
				r = 0.5*(x2 + r);
			    else if (q == 0.0)
				do_newton = NO;
			}
		    }
		}
	    }
	    else
	    {
		q = 1.0 - S1 - S2 - S3; /* p1(1) */
		if (q < 0.0)
		    do_newton = YES;
		else if (q == 0.0)
		    r = 1.0;
	    }
	    if (do_newton)
	    {
		for (i = 0; i < num_iter; ++i)
		{
		    p1r = 1.0 - r*(S1 + r*(S2 + S3*r));
		    p1rp = -(S1 + r*(2.0*S2 + 3.0*r*S3));
		    r -= p1r/p1rp;
		}
	    }
	    x_max = min(r,1.0);
	}
	(void) printf("Reference pressure denominator vanishes at\n"
		      "    compression x_max = %"FFMT",\n"
		      "    density = %"FFMT"\n",x_max,
		      (x_max < 1.0) ? mg->Rho0/(1.0 - x_max) : HUGE_VAL);
	
	/* Now determine interval for real sound speeds */

	p2[0] = 1.0;
	p2[1] = S1 - 1.0;
	p2[2] = 3.0*S2 - (gamma0 + 1.0)*S1;
	p2[3] = (gamma0-b)*S1 - (2.0*gamma0 + 3.0)*S2 + 5.0*S3;
	p2[4] = 2.0*(b-gamma0)*S2 + (5.0 + 3.0*gamma0)*S3;
	p2[5] = 3.0*(gamma0-b)*S3;
	c2x_max = x_max;
	for (i = 0; i < 5; ++i)
	{
	    double a[5];
	    int   k;
	    (void) p2_aux(c2x_max,&q,(POINTER)p2);
	    if (debugging("eosref"))
	    {
		(void) printf("\np2 = (\n"
			      "      %"FFMT"\n"
			      "      %"FFMT"\n"
			      "      %"FFMT"\n"
			      "      %"FFMT"\n"
			      "      %"FFMT"\n"
			      "      %"FFMT"\n"
			      "     )\n",
			      p2[0],p2[1],p2[2],p2[3],p2[4],p2[5]);
		(void) printf("p2_aux(%"FFMT") = %"FFMT"\n",c2x_max,q);
	    }
	    if (q < 0.0)
	    {
	        if (!find_root(p2_aux,(POINTER)p2,0.0,&r,0.0,c2x_max,
				 epsilon,epsilon))
		{
		    screen("ERROR in max_compression(), find_root() failed "
			   "for q < 0.0\n");
		    clean_up(ERROR);
		}
	        c2x_max = r;
	        if (debugging("eosref"))
		{
		    (void) printf("root c2x_max = %"FFMT"\n",c2x_max);
	            (void) p2_aux(c2x_max,&q,(POINTER)p2);
		    (void) printf("p2_aux(%"FFMT") = %"FFMT"\n",c2x_max,q);
		}
	    }
	    else if (find_separation_point(p2_aux,(POINTER)p2,0.0,&r,0.0,
			                   c2x_max,&xp2min,&p2min,epsilon))
	    {
	        if (!find_root(p2_aux,(POINTER)p2,0.0,&r,
				 0.0,r,epsilon,epsilon))
		{
		    screen("ERROR in max_compression(), find_root() failed "
			   "for q >= 0.0\n");
		    clean_up(ERROR);
		}
	        c2x_max = r;
	    }
	    else
		break;

	    a[4] = p2[5];
	    for (k = 3; k >= 0; --k)
	        a[k] = (p2[k+1] + r*a[k+1]);
	    p2[5] = 0.0;
	    for (k = 0; k < 5; ++k)
	        p2[k] = -a[k];
	}
	if (c2x_max < 1.0)
	{
	    double Gm, pr, dpr, der;
	    rho = mg->Rho0/(1.0 - c2x_max);
	    Gm = Gamma(rho,mg);
	    pr = P_ref(rho,mg);
	    dpr = dP_ref(rho,mg);
	    der = dE_ref(rho,mg);
	    c2ref = (Gm/rho)*(pr+der) - dpr/(rho*rho);
	}
	else
	{
	    rho = HUGE_VAL;
    	    c2ref = 0.0;
	}
	(void) printf("Reference sound speed vanishes at\n"
		      "    compression x_max = %"FFMT",\n"
		      "    density = %"FFMT",\n"
		      "    c2ref = %"FFMT"\n",
		      c2x_max,rho,c2ref);
	debug_print("eosref","Left max_compression\n");
	return min(x_max,c2x_max);
}		/*end max_compression*/

LOCAL	boolean p2_aux(
	double x,
	double *p2x,
	POINTER prms)
{
    	double *p2 = (double*)prms;
	*p2x = p2[0] + x*(p2[1] + x*(p2[2] + x*(p2[3] + x*(p2[4] + x*p2[5]))));
	return FUNCTION_SUCCEEDED;
}		/*end p2_aux*/

/***************END INITIALIZATION UTILITY FUNCTIONS************************/

/***************EQUATION OF STATE DOMAIN FUNCTIONS**************************/

LOCAL	double	MG_Min_energy(
	Locstate	state)
{
	double emin = -HUGE_VAL;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double rho, Gm, dGm, pr, dpr, er, der;
	double Cv, Tr, Tmin;
	double c2ref;
	double min_energy;
    	MG_params *mg;
	mg = &MG_Eos(state)->MGparams;
	rho = Dens(state);
	Gm = Gamma(rho,mg);
	dGm = dGamma(rho,mg);
	pr = P_ref(rho,mg);
	dpr = dP_ref(rho,mg);
	er = E_ref(rho,mg);
	der = dE_ref(rho,mg);
	Tr = T_ref(rho,mg);
	Cv = mg->CV;
	c2ref = (Gm/rho)*(pr+der) - dpr/(rho*rho);
	emin = er + (MACH_EPS - c2ref)/(Gm*(Gm+1.0) - dGm/rho);
	Tmin = Tr + (emin - er)/Cv;
	if (Tmin < MACH_EPS)
	    emin = er + Cv*(MACH_EPS - Tr);
	emin *= rho;
	min_energy = Params(state)->min_energy;
	if (emin < min_energy)
	    emin = min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return emin;
}	/*end MG_Min_energy*/

/*
*			MG_Min_pressure():
*
*	Selects the minimum pressure at the specified density of the state so
*	that both the temperature is positive and the sound
*	speed is real.
*/

LOCAL	double	MG_Min_pressure(
	Locstate	state)
{
	double pmin = -HUGE_VAL;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double rho, Gm, dGm, pr, dpr, der, Tr, Tmin, Cv;
	double c2ref;
	double min_pressure;
    	MG_params *mg;
	mg = &MG_Eos(state)->MGparams;
	rho = Dens(state);
	Gm = Gamma(rho,mg);
	dGm = dGamma(rho,mg);
	pr = P_ref(rho,mg);
	dpr = dP_ref(rho,mg);
	der = dE_ref(rho,mg);
	Tr = T_ref(rho,mg);
	Cv = mg->CV;
	c2ref = (Gm/rho)*(pr+der) - dpr/(rho*rho);
	pmin = pr+rho*Gm*(MACH_EPS - c2ref)/(Gm*(Gm+1.0)-dGm/rho);
	Tmin = Tr + (pmin - pr)/(Cv*rho*Gm);
	if (Tmin < MACH_EPS)
	    pmin = pr + rho*Gm*Cv*(MACH_EPS - Tr);
	min_pressure = Params(state)->min_pressure;
	if (pmin < min_pressure)
	    pmin = min_pressure;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pmin;
}	/*end MG_Min_pressure*/

LOCAL	double	MG_Min_temperature(
	Locstate	state)
{
	double Tmin = -HUGE_VAL;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	double rho, Gm, dGm, pr, dpr, der;
	double Cv, Tr;
	double c2ref;
    	MG_params *mg;
	mg = &MG_Eos(state)->MGparams;
	rho = Dens(state);
	Gm = Gamma(rho,mg);
	dGm = dGamma(rho,mg);
	pr = P_ref(rho,mg);
	dpr = dP_ref(rho,mg);
	der = dE_ref(rho,mg);
	Tr = T_ref(rho,mg);
	Cv = mg->CV;
	c2ref = (Gm/rho)*(pr+der) - dpr/(rho*rho);
	Tmin = Tr + (MACH_EPS - c2ref)/(Cv*(Gm*(Gm+1.0) - dGm/rho));
	if (Tmin < MACH_EPS)
	    Tmin = MACH_EPS;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return Tmin;
}	/*end MG_Min_temperature*/

/***************END EQUATION OF STATE DOMAIN FUNCTIONS*************************/


/*
*		MG_pressure_on_adiabat():
*
*	Returns the pressure of the state with density rho on the adiabat
*	through the state with density rho1 and temperature T1.
*/

LOCAL double MG_pressure_on_adiabat(
	double     rho,
	double     rho1,
	double     T1,
	MG_params *mg)
{
    	double Cv, Gm, pr, tr;

	Cv = mg->CV;
	pr = P_ref(rho,mg);
	tr = T_ref(rho,mg);
	Gm = Gamma(rho,mg);
	return pr + Cv*Gm*rho*(T1*IrhoG(rho,rho1,mg) - tr);
}		/*end MG_pressure_on_adiabat*/

/*
*		MG_sound_speed_squared_on_adiabat():
*
*	Returns the square of the sound speed of the state with
*	density rho on the adiabat through the state with density rho1
*	and temperature T1.
*/


LOCAL double MG_sound_speed_squared_on_adiabat(
	double     rho,
	double     rho1,
	double     T1,
	MG_params *mg)
{
	double Cv, Gm, dGm, pr, dpr, der, Tr, c2r;

	if (rho <= 0.0)
	    return 0.0;

	Cv = mg->CV;
	Gm = Gamma(rho,mg);
	dGm = dGamma(rho,mg);
	pr = P_ref(rho,mg);
	dpr = dP_ref(rho,mg);
	der = dE_ref(rho,mg);
	Tr = T_ref(rho,mg);
	c2r = Gm*(pr + der)/rho - dpr/(rho*rho);
    
	return  fabs(c2r +
		     (Gm+1.0-dGm/(Gm*rho))*Gm*Cv*(IrhoG(rho,rho1,mg)*T1 - Tr));
}		/*end MG_sound_speed_squared_on_adiabat*/

LOCAL	boolean test_out_of_range(
	double      rho,
	Locstate   state,
	const char *fname,
	boolean       fatal)
{
    	boolean out_of_range;
	MG_params *mg;
	mg = &MG_Eos(state)->MGparams;
	out_of_range = (rho > mg->Rho_max) ? YES : NO;
	if (out_of_range && fatal && MG_OOR)
	{
	    boolean bin_out;
	    screen("ERROR in %s, density out of range of valid EOS\n",fname);
	    screen("rho = %"FFMT", Rho_max = %"FFMT"\n",rho,mg->Rho_max);
	    bin_out = is_binary_output();
	    set_binary_output(NO);
	    fprint_gas_data(stdout,state);
	    fprint_Gas_param(stdout,Params(state));
	    set_binary_output(bin_out);
	    clean_up(ERROR);
	}
	return out_of_range;
}		/*end test_out_of_range*/

/* All functions that depend explicitly on the reference curves follow */

/*
*			Gamma():
*
*	Return the Gruneisen exponent at the given density
*/

LOCAL	double Gamma(
	double     rho,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double gamma0 = mg->gamma0;
	double b = mg->b;
	double xi;

	if (rho >= Rho0)
	{
	    xi = (1.0 - Rho0/rho);
	    return gamma0*(1.0 - xi) + b*xi;
	}
	else if (rho >= 0.0)
	    return gamma0;
	else
	{
	    screen("ERROR in Gamma(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end Gamma*/

/*
*			dGamma():
*
*	Return the derivative of the Gruneisen exponent with respect
*	to specific volume at the given density.
*/

LOCAL	double dGamma(
	double     rho,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double gamma0 = mg->gamma0;
	double b = mg->b;

	if (rho >= Rho0)
	    return Rho0*(gamma0 - b);
	else if (rho >= 0.0)
	    return 0.0;
	else
	{
	    screen("ERROR in dGamma(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end dGamma*/

/*
*			d2Gamma():
*
*	Return the second derivative of the Gruneisen exponent with respect
*	to specific volume at the given density.
*
*	Note: for the specific Mie-Gruneisen equation state used here
*	this quantity is identically zero.  However it is included as a
*	separate function since all of the functions that use d2Gamma()
*	are written in a form that only assumes the Gruneisen exponent is a
*	given function of specific volume.
*/

/*ARGSUSED*/
LOCAL	double d2Gamma(
	double     rho,
	MG_params *mg)
{
	if (rho >= 0.0)
	    return 0.0;
	else
	{
	    screen("ERROR in d2Gamma(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end d2Gamma*/

/*
*			P_ref():
*
*	Evaluates pressure along the reference curve.
*/

LOCAL	double P_ref(
	double     rho,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double xi, S, C0, P_inf, P0, D;

	P0 = mg->P0;
	if (rho >= Rho0)
	{
	    C0 = mg->C0;
	    xi = (1.0 - Rho0/rho);
	    S = mg->S1 + mg->S2*xi + mg->S3*xi*xi;
	    D = 1.0 - S*xi;
	    return P0 + Rho0*C0*C0*xi/(D*D);
	}
	else if (rho >= 0.0)
	{
	    P_inf = mg->P_inf;
	    return P0 + (P0 + P_inf)*(pow(rho/Rho0,mg->gamma0+1.0) - 1.0);
	}
	else
	{
	    screen("ERROR in P_ref(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end P_ref*/

/*
*			dP_ref():
*
*	Evaluates the derivative of pressure with respect
*	to specific volume along the reference curve.
*/

LOCAL	double dP_ref(
	double     rho,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double xi, S, dS, C0, D;

	C0 = mg->C0;
	if (rho >= Rho0)
	{
	    xi = (1.0 - Rho0/rho);
	    S  = mg->S1 + mg->S2*xi +     mg->S3*xi*xi;
	    dS =          mg->S2    + 2.0*mg->S3*xi;
	    D = 1.0 - S*xi;
	    return -Rho0*Rho0*C0*C0*(1.0 + xi*(S + 2.0*xi*dS))/(D*D*D);
	}
	else if (rho >= 0.0)
	{
	    return -sqr(Rho0*C0)*pow(rho/Rho0,mg->gamma0+2.0);
	}
	else
	{
	    screen("ERROR in dP_ref(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end dP_ref*/

/*
*			d2P_ref():
*
*	Evaluates the second derivative of pressure with respect
*	to specific volume along the reference curve.
*/

LOCAL	double d2P_ref(
	double     rho,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double xi, S, dS, d2S, C0, D;

	C0 = mg->C0;
	if (rho >= Rho0)
	{
	    xi = (1.0 - Rho0/rho);
	    S   = mg->S1 + mg->S2*xi +     mg->S3*xi*xi;
	    dS  =          mg->S2    + 2.0*mg->S3*xi;
	    d2S =                      2.0*mg->S3;
	    D = 1.0 - S*xi;
	    return 2.0*Rho0*Rho0*Rho0*C0*C0*
	    (2.0*S+xi*(S*S+4.0*dS+xi*(d2S+2.0*S*dS+xi*(3.0*dS*dS-S*d2S))))/
		   (D*D*D*D);
	}
	else if (rho >= 0.0)
	{
	    double gamma0 = mg->gamma0;
	    return Rho0*sqr(Rho0*C0)*(gamma0+2.0)*pow(rho/Rho0,gamma0+3.0);
	}
	else
	{
	    screen("ERROR in d2P_ref(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end d2P_ref*/

/*
*			E_ref():
*
*	Evaluates specific internal energy along the reference curve.
*/

LOCAL	double E_ref(
	double     rho,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double xi, S, C0, P0, P_inf, E0, D;

	C0 = mg->C0;
	P0 = mg->P0;
	E0 = mg->E0;
	if (rho >= Rho0)
	{
	    xi = (1.0 - Rho0/rho);
	    S = mg->S1 + mg->S2*xi + mg->S3*xi*xi;
	    D = 1.0 - S*xi;
	    return E0 + (P0/Rho0)*xi + 0.5*C0*C0*xi*xi/(D*D);
	}
	else if (rho > 0.0)
	{
	    double gamma0 = mg->gamma0;
	    P_inf = mg->P_inf;
	    return  E0 + ((P0+P_inf)/(Rho0*gamma0))*(pow(rho/Rho0,gamma0) - 1.0)
		       + P_inf*(1.0/rho - 1.0/Rho0);
	}
	else if (rho == 0.0)
	{
	    if (mg->P_inf > 0.0)
	        return HUGE_VAL;
	    else if (mg->P_inf < 0.0)
	        return -HUGE_VAL;
	    else
		return 0.0;
	}
	else
	{
	    screen("ERROR in E_ref(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end E_ref*/

/*
*			dE_ref():
*
*	Evaluates the derivative of specific internal energy with respect
*	to specific volume along the reference curve.
*/

LOCAL	double dE_ref(
	double     rho,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double xi, S, dS, C0, P_inf, P0, D;

	P0 = mg->P0;
	if (rho >= Rho0)
	{
	    C0 = mg->C0;
	    xi = (1.0 - Rho0/rho);
	    S  = mg->S1 + mg->S2*xi +     mg->S3*xi*xi;
	    dS =          mg->S2    + 2.0*mg->S3*xi;
	    D = 1.0 - S*xi;
	    return -P0 - Rho0*C0*C0*xi*(1.0 + xi*xi*dS)/(D*D*D);
	}
	else if (rho >= 0.0)
	{
	    P_inf = mg->P_inf;
	    return -P0 - (P0 + P_inf)*(pow(rho/Rho0,mg->gamma0+1.0) - 1.0);
	}
	else
	{
	    screen("ERROR in dE_ref(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end dE_ref*/

/*
*			d2E_ref():
*
*	Evaluates the second derivative of specific internal energy with
*	respect to specific volume along the reference curve.
*/

LOCAL	double d2E_ref(
	double     rho,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double xi, S, dS, d2S, C0, D;

	C0 = mg->C0;
	if (rho >= Rho0)
	{
	    xi = (1.0 - Rho0/rho);
	    S  = mg->S1 + mg->S2*xi +     mg->S3*xi*xi;
	    dS =          mg->S2    + 2.0*mg->S3*xi;
	    d2S =                     2.0*mg->S2;
	    D = 1.0 - S*xi;
	    return Rho0*Rho0*C0*C0*
		   (1.0+xi*(2.0*S+xi*(6.0*dS+xi*(d2S+(3.0*dS*dS-S*d2S)*xi))))/
		          (D*D*D*D);
	}
	else if (rho >= 0.0)
	{
	    return sqr(Rho0*C0)*pow(rho/Rho0,mg->gamma0+2.0);
	}
	else
	{
	    screen("ERROR in d2E_ref(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    return -HUGE_VAL;
	}
}		/*end d2E_ref*/

/*
*			T_ref():
*
*	Evaluates temperature along the reference curve.
*
*	T_ref() is computed as the solution evaluated at
*	V = 1/rho of the differential equation
*
*	dTr/dV + (Gamma(V)/V)*Tr = (P_ref(V) + dE_ref(V))/CV
*/

#if defined(__cplusplus)
extern "C" {    
#endif /* defined(__cplusplus) */
    LOCAL   LSODE_FUNC  Trefrhs;
    LOCAL   LSODE_JAC 	jacTrefrhs;
#if defined(__cplusplus)
}                          
#endif /* defined(__cplusplus) */

LOCAL	MG_params *mgTref;

LOCAL	double T_ref(
	double     rho,
	MG_params *mg)
{
	double T;
	double Rho0 = mg->Rho0;
	double g0 = mg->gamma0;

	if (rho < 0.0)
	{
	    screen("ERROR in T_ref(), negative density %"FFMT"\n",rho);
	    print_MG_params(mg);
	    clean_up(ERROR);
	    T = -HUGE_VAL;
	}
	else if (rho < Rho0)
	    T = mg->T0*pow(rho/Rho0,g0);
	else
	{
	    Tref_Spline *TrefSpline = mg->TrefSpline;
	    double  rho_max = TrefSpline->rho[TrefSpline->m-1];

	    if (rho <= rho_max)
	    {
	        int   one = 1, ier;
	        FORTRAN_NAME(splev)(TrefSpline->t,&TrefSpline->n,TrefSpline->c,
		                    &TrefSpline->k,&rho,&T,&one,&ier);
	        if (ier != 0)
	        {
		    int i;
		    screen("Error in T_ref(), invalid spline evaluation, "
		           "rho = %"FFMT", ier = %d\n",rho,ier);
		    (void) printf("Rho0 = %"FFMT", Rho_max = %"FFMT"\n",
				  Rho0,mg->Rho_max);
		    (void) printf("TrefSpline->n = %d, TrefSpline->k = %d\n",
			          TrefSpline->n,TrefSpline->k);
		    output();
	            (void) printf("%-30s %-30s Spline_Coef\n",
			          "knots","B-Spline_Coefs");
	            for (i = 0; i < TrefSpline->n; ++i)
	            (void) printf(" %-30.24g %-30.24g\n",
			          TrefSpline->t[i],TrefSpline->c[i]);
	            (void) printf("\n");

		    clean_up(ERROR);
	        }
	    }
	    else if ((rho < mg->Rho_max) && (rho < HUGE_VAL))
	    {
		extend_spline_Tref_curve(rho,mg);
		return T_ref(rho,mg);
	    }
	    else
		return HUGE_VAL;
	}
	return T;
}		/*end T_ref*/

LOCAL	void extend_spline_Tref_curve(
	double     up_to_rho,
	MG_params *mg)
{
	Tref_Spline *TrefSpline = mg->TrefSpline;
    	double       rho_max, old_rho_max;
	double       Rho0 = mg->Rho0;
	double       *rho, *old_rho;
	double       *T, *old_T;
	double       drho;
	double       rtol, atol, terr;
	double       *w, *t, *c;
	double       *wrk;
	double       s, fp;
	double       xb, xe;
	double       N;
	int         iopt;
	int         nest;
	int         *iwrk;
	int   	    lwrk;
	int         old_m, m;
	int         i, k, n, istate, ier;

	old_m = TrefSpline->m;
	old_rho = TrefSpline->rho;
	old_T = TrefSpline->T;
	drho = old_rho[1] - old_rho[0];
	old_rho_max = TrefSpline->rho[old_m-1];

	N = ceil(log((up_to_rho - Rho0)/(old_rho_max - Rho0))/log(2.0));
	rho_max = old_rho_max + pow(2.0,N)*(old_rho_max - Rho0);
	rho_max = min(mg->Rho_max,rho_max);
	m = (int) ceil((rho_max - Rho0)/drho);
	if ((Rho0 + (m-1)*drho) < up_to_rho)
	{
	    screen("ERROR in extend_spline_Tref_curve(), "
		   "can't extend curve to rho = %"FFMT" due to limit %"FFMT" "
		   "on the reference curves\n",up_to_rho,mg->Rho_max);
	    clean_up(ERROR);
	}

	uni_array(&rho,m,FLOAT);
	uni_array(&T,m,FLOAT);
	for (i = 0; i < old_m; ++i)
	{
	    rho[i] = old_rho[i];
	    T[i] = old_T[i];
	}
	free_these(2,TrefSpline->rho,TrefSpline->T);
	TrefSpline->rho = rho;
	TrefSpline->T = T;
	for (i = old_m; i < m; ++i)
	    rho[i] = Rho0 + i*drho;
	rho[m-1] = rho_max;

	atol = mg->TrefSpline->Tref_atol;
	rtol = mg->TrefSpline->Tref_rtol;
	terr = mg->TrefSpline->Tref_terr;
	mgTref = mg;
	istate = 1;
	if ((!ode_solver(Trefrhs,T+old_m-1,rho[old_m-1],drho,m-old_m+1,rtol,atol,
		         terr,&istate,jacTrefrhs)) || (istate != 2))
	{
	    screen("ERROR in extend_spline_Tref_curve() can't solve ODE\n"
	           "error flag = %d\n",istate);
	    clean_up(ERROR);
	}

	for (i = old_m; i < m; ++i)
	{
	    if (log(T[i]) > 2*log(T[i-1]))
	    {
		(void) printf("WARNING in extend_spline_Tref_curve(), "
			      "resetting Rho_max from %"FFMT" to %"FFMT"\n",
			      mg->Rho_max,rho[i-1]);
		(void) printf("T_ref(%"FFMT") = %"FFMT", "
			      "T_ref(%"FFMT") = %"FFMT"\n",
			      rho[i-1],T[i-1],rho[i],T[i]);
	        m = i;
	        mg->Rho_max = rho_max = rho[m-1];
	        mg->V_min = 1.0/mg->Rho_max;
		break;
	    }
	}
	if (debugging("mg-eos"))
	{
	    output();
	    (void) printf("%-30s %-30s %-30s Temperature_Reference_Curve\n",
			  "Density","Temperature","logTemperature");
	    for (i = 0; i < m; ++i)
	        (void) printf(" %-30.24g %-30.24g %-30.24g\n",
			      rho[i],T[i],log(T[i]));
	    (void) printf("\n");
	}

	k = TrefSpline->k;
	TrefSpline->m = m;
	TrefSpline->nest = nest = m + k + 1;

	free_these(2,TrefSpline->t,TrefSpline->c);
	uni_array(&TrefSpline->t,nest,FLOAT);
	uni_array(&TrefSpline->c,nest,FLOAT);
	t = TrefSpline->t;
	c = TrefSpline->c;
	/* Work space storage */
	uni_array(&w,m,FLOAT);
	for (i = 0; i < m; ++i)
	    w[i] = 1.0;
	lwrk = m*(k+1)+nest*(7+3*k);
	uni_array(&wrk,lwrk,FLOAT);
	uni_array(&iwrk,nest,FLOAT);

	xb = Rho0; xe = rho_max;
	for (i = 0; i < k; ++i)
	    t[i] = xb;
	for (i = k; i < (nest-k); ++i)
	    t[i] = rho[i-k];
	for (i = nest-k; i < nest; ++i)
	    t[i] = xe;

	s = 0.0; /* s = 0 implies interpolation */
	iopt = 0;/* Smoothing spline */
	FORTRAN_NAME(curfit)(&iopt,&m,rho,T,w,&xb,&xe,&k,&s,&nest,&n,t,c,&fp,
		             wrk,&lwrk,iwrk,&ier);

	switch (ier)
	{
	case 0:
	    break;
	case -1:
	    break;
	case -2:
	    break;
	case 1:
	    screen("ERROR in extend_spline_Tref_curve() spline fit failed\n"
	           "error flag = %d, not enough storage\n",ier);
	    clean_up(ERROR);
	    break;
	case 2:
	    screen("ERROR in extend_spline_Tref_curve() spline fit failed\n"
	           "error flag = %d, s too small?\n",ier);
	    clean_up(ERROR);
	    break;
	case 3:
	    screen("ERROR in extend_spline_Tref_curve() spline fit failed\n"
	           "error flag = %d, max interations exceeded "
		   "s too small?\n",ier);
	    clean_up(ERROR);
	    break;
	case 10:
	    screen("ERROR in extend_spline_Tref_curve() spline fit failed\n"
	           "error flag = %d, invalid input data\n",ier);
	    (void) printf("m = %d, nest = %d\n",m,nest);
	    (void) printf("drho = %"FFMT"\n",drho);
	    (void) printf("%-10s %-"SFMT" %-"SFMT"\n","i","rho","T");
	    for (i = 0; i < m; ++i)
	        (void) printf("%-10d %-"FFMT" %-"FFMT"\n",i,rho[i],T[i]);
	    clean_up(ERROR);
	    break;
	default:
	    screen("ERROR in extend_spline_Tref_curve() spline fit failed\n"
	           "error flag = %d, unknown error\n",ier);
	    clean_up(ERROR);
	    break;
	}

	TrefSpline->n = n;

	free_these(3,w,wrk,iwrk);
}		/*end extend_spline_Tref_curve*/

LOCAL	void initialize_Tref(
	MG_params *mg)
{
	Tref_Spline *TrefSpline;
	double       Rho0, rho_max;
	double       *rho, *T, drho;
	double       rtol, atol, terr;
	double       *w, *t, *c;
	double       xb, xe;
	double       s, fp;
	double       *wrk;
	int         *iwrk;
	int         lwrk;
	int         i, m, n, istate;
	int         iopt, k, nest;
	int         ier;

	Rho0 = mg->Rho0;
	rho_max = min(mg->Rho_max,11.0*Rho0);
	xb = Rho0; xe = rho_max;

	m = 1001; /*TOLERANCE*/
	iopt = 0;/* Smoothing spline */
	scalar(&TrefSpline,sizeof(Tref_Spline));
	mg->TrefSpline = TrefSpline;
	TrefSpline->m = m;
	uni_array(&TrefSpline->rho,m,FLOAT);
	uni_array(&TrefSpline->T,m,FLOAT);
	TrefSpline->k = k = 3; /* Cubic spline approximates T_ref */
	TrefSpline->nest = nest = m + k + 1;
	uni_array(&TrefSpline->t,nest,FLOAT);
	uni_array(&TrefSpline->c,nest,FLOAT);


	/* Work space storage */
	uni_array(&w,m,FLOAT);
	for (i = 0; i < m; ++i)
	    w[i] = 1.0;
	lwrk = m*(k+1)+nest*(7+3*k);
	uni_array(&wrk,lwrk,FLOAT);
	uni_array(&iwrk,nest,FLOAT);

	rho = TrefSpline->rho;
	T = TrefSpline->T;
	t = TrefSpline->t;
	c = TrefSpline->c;

	drho = (rho_max - Rho0)/(m - 1);
	rho[0] = Rho0;
	for (i = 1; i < m; ++i)
	    rho[i] = Rho0 + i*drho;
	rho[m-1] = rho_max;

	for (i = 0; i < k; ++i)
	    t[i] = xb;
	for (i = k; i < (nest-k); ++i)
	    t[i] = rho[i-k];
	for (i = nest-k; i < nest; ++i)
	    t[i] = xe;

	T[0] = mg->T0;
	mgTref = mg;
	mg->TrefSpline->Tref_rtol = rtol = 1.0e-6;/*TOLERANCE*/
	mg->TrefSpline->Tref_atol = atol = rtol*mg->T0;/*TOLERANCE*/
	mg->TrefSpline->Tref_terr = terr = 1.0e-3*drho;/*TOLERANCE*/
	istate = 1;
	if ((!ode_solver(Trefrhs,T,Rho0,drho,m,rtol,atol,
		         terr,&istate,jacTrefrhs)) || (istate != 2))
	{
	    screen("ERROR in initialize_Tref() can't solve ODE\n"
	           "error flag = %d\n",istate);
	    clean_up(ERROR);
	}


	for (i = 1; i < m; ++i)
	{
	    if (log(T[i]) > 2*log(T[i-1]))
	    {
		(void) printf("WARNING in initialize_Tref(), "
			      "resetting Rho_max from %"FFMT" to %"FFMT"\n",
			      mg->Rho_max,rho[i-1]);
		(void) printf("T_ref(%"FFMT") = %"FFMT", "
			      "T_ref(%"FFMT") = %"FFMT"\n",
			      rho[i-1],T[i-1],rho[i],T[i]);
	        TrefSpline->m = m = i;
	        mg->Rho_max = rho_max = rho[m-1];
	        mg->V_min = 1.0/mg->Rho_max;
	        TrefSpline->nest = nest = m + k + 1;
		break;
	    }
	}
	if (debugging("mg-eos"))
	{
	    output();
	    (void) printf("%-30s %-30s %-30s Temperature_Reference_Curve\n",
			  "Density","Temperature","logTemperature");
	    for (i = 0; i < m; ++i)
	        (void) printf(" %-30.24g %-30.24g %-30.24g\n",
			      rho[i],T[i],log(T[i]));
	    (void) printf("\n");
	}

	s = 0.0; /* s = 0 implies interpolation */
	FORTRAN_NAME(curfit)(&iopt,&m,rho,T,w,&xb,&xe,&k,&s,&nest,&n,t,c,&fp,
		             wrk,&lwrk,iwrk,&ier);

	switch (ier)
	{
	case 0:
	    break;
	case -1:
	    break;
	case -2:
	    break;
	case 1:
	    screen("ERROR in initialize_Tref() spline fit failed\n"
	           "error flag = %d, not enough storage\n",ier);
	    clean_up(ERROR);
	    break;
	case 2:
	    screen("ERROR in initialize_Tref() spline fit failed\n"
	           "error flag = %d, s too small?\n",ier);
	    clean_up(ERROR);
	    break;
	case 3:
	    screen("ERROR in initialize_Tref() spline fit failed\n"
	           "error flag = %d, max interations exceeded "
		   "s too small?\n",ier);
	    clean_up(ERROR);
	    break;
	case 10:
	    screen("ERROR in initialize_Tref() spline fit failed\n"
	           "error flag = %d, invalid input data\n");
	    clean_up(ERROR);
	    break;
	default:
	    screen("ERROR in initialize_Tref() spline fit failed\n"
	           "error flag = %d, unknown error\n",ier);
	    clean_up(ERROR);
	    break;
	}

	TrefSpline->n = n;
	free_these(3,w,wrk,iwrk);
}		/*end initialize_Tref*/

#if defined(__cplusplus)
extern "C" {    
#endif /* defined(__cplusplus) */

/*ARGSUSED*/
LOCAL	void Trefrhs(
	int		*neq,
	double		*x,
	double		*y,
	double		*yp)
{
    	double rho = x[0];
	double T = y[0];
	double Gm;
	double CV = mgTref->CV;

	Gm = Gamma(rho,mgTref);
	yp[0] = (Gm/rho)*T -
	    (1.0/(CV*sqr(rho)))*(P_ref(rho,mgTref)+dE_ref(rho,mgTref));
}		/*end Trefrhs*/

/*ARGSUSED*/
LOCAL	void	jacTrefrhs(
	int	*neq,
	double	*x,
	double	*y,
	int	*ml,
	int	*mu,
	double	*pd,
	int	*nrowpd)
{
    	double rho = x[0];
	double Gm;

	Gm = Gamma(rho,mgTref);
	*pd = Gm/rho;
}	/*end jacTrefrhs*/

#if defined(__cplusplus)
}                          
#endif /* defined(__cplusplus) */

/*
*			IrhoG():
*
*	Returns the integral from V0 to V of the Gruneisent exponent divided
*	by specific volume with respect to V,
*
*				--                  --
*			        | / V2               |
*	IrhoG(rho1,rho2,mg) =exp| \    Gamma(v)/v dv |, V1 = 1/rho1, V2 = 1/rho2
*			        | /V1                |
*				--                  --
*
*	V0 = 1/rho0,  the reference state density and V = 1/rho
*/

LOCAL	double IrhoG(
	double     rho1,
	double     rho2,
	MG_params *mg)
{
	double Rho0 = mg->Rho0;
	double b = mg->b;
	double g0 = mg->gamma0;

	if ((rho1 < 0.0) || (rho2 < 0.0))
	{
	    screen("ERROR in IrhoG(), negative density\n");
	    (void) printf("rho1 = %"FFMT", rho2 = %"FFMT"\n",rho1,rho2);
	    print_MG_params(mg);
	    clean_up(ERROR);
	}
	if ((rho1 >= Rho0) && (rho2 >= Rho0))
	{
	    return pow(rho1/rho2,b)*exp((g0-b)*(Rho0/rho2 - Rho0/rho1));
	}
	else if (rho1 < Rho0)
	{
	    if (rho2 >= Rho0)
	        return pow(Rho0/rho2,b)*exp((g0-b)*(Rho0/rho2 - 1.0))*pow(rho1/Rho0,g0);
	    else if (g0 == 0.0)
		return 1.0;
	    else if ((rho2 > 0.0) && (rho1 > 0.0))
		return pow(rho1/rho2,g0);
	    else if (rho2 > 0.0) /*rho1 == 0*/
		return (g0 > 0.0) ? 0.0 : HUGE_VAL;
	    else if (rho1 > 0.0) /*rho2 == 0*/
		return (g0 > 0.0) ? HUGE_VAL : 0.0;
	    else
		return 1.0;
	}
	else
	{
	    return pow(rho1/Rho0,b)*exp((g0-b)*(1.0 - Rho0/rho1))*pow(Rho0/rho2,g0);
	}
}		/*end IrhoG*/
