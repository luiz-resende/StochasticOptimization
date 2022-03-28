#ifndef A4_MAIN_CPP
#define A4_MAIN_CPP

#include "A0_definit.h"


clock_t		start;


int main()
{
	matrix_flt	xi_volum = {			// Possible precipitation volumes for all months
		{ -50.0, 100.0, 250.0 },
		{ -50.0, 100.0, 250.0 },
		{ 50.0, 150.0, 350.0 },
		{ 50.0, 150.0, 350.0 }
	};	

	matrix_flt	xi_probs = {			// Precipitation probabilities for all months
		{ 0.25, 0.50, 0.25 },
		{ 0.25, 0.50, 0.25 },
		{ 0.25, 0.50, 0.25 },
		{ 0.25, 0.50, 0.25 }
	};

	try {
		double probability_b1 = accumulate(xi_probs[0].begin(), xi_probs[0].end(), 0.0);
		double probability_b2 = accumulate(xi_probs[2].begin(), xi_probs[2].end(), 0.0);
		string except = "IOException: ";

		if ((probability_b1 != 1.0) && (probability_b2 != 1.0)) {
			except += "The probabilities of precipitation for months in the first and second blocks do not add to 1.0...\n\n";
			throw except;
		}
		else if ((probability_b1 != 1.0) && (probability_b2 == 1.0)) {
			except += "The probabilities of precipitation for months in the first block do not add to 1.0...\n\n";
			throw except;
		}
		else if ((probability_b1 == 1.0) && (probability_b2 != 1.0)) {
			except += "The probabilities of precipitation for months in the second block do not add to 1.0...\n\n";
			throw except;
		}
	}
	catch (string e) {
		printf("%s", e.c_str());
		system("PAUSE");
		return 1;
	}

	/* INITIALIZING THE GLOBAL VARIABLES */
	int		T = 4;						// The total number of periods t = 1, ..., T
	int		S = 3;						// The number of possible weather realizations per period t

	double	R_max = 200.0;				// The maximum possible level the dam can release in a given period
	double	L_max = 250.0;				// The maximum level below flood stage the dam level could be
	double	L_0 = (L_max - 150.0);		// Initial level below flood stage

	double	cost_lowering = 0.0;		// The cost per millimetre of dam level lowered
	double	cost_import = 5000.0;		// The cost per millimetre of water imported
	double	cost_flood = 10000.0;		// The cost per millimetre of flood

	int		objct_type = 1;

	/* BLOCK FOR SOLVING PROBLEM USING CPLEX MODEL */
	printf("Solving problem with CPLEX...\n\n");
	double opt_value_RP, opt_value_expect, opt_value_determ, opt_value_evpi;
	try {
		opt_value_RP = Dam_CPLEX_RP(
			T,				// The total number of periods
			S,				// The number of scenarios in each period
			R_max,			// The maximum possible level the dam can release in a given period
			L_max,			// The maximum level below flood stage the dam level could be lowered in a given period
			L_0,			// The initial/current level below flood stage the dam's water level is in
			cost_lowering,	// The cost per millimetre of dam level lowered
			cost_flood,		// The cost per millimetre of water flood
			cost_import,	// The cost per millimetre of water imported
			xi_volum,		// The vector with possible precipitation volumes for months in the first block
			xi_probs,		// The vector with precipitation probabilities for months in the first block
			objct_type,		// The objective function type (1-minimize expected cost; 2-minimize probability)
			false,			// Whether or not show the display
			true,			// Whether or not to save the result to files
			2,				// Information display intensity (default is 2)
			-1.0,			// Percentage tolerance gap (default is -1.0)
			-1,				// Time limit in seconds (default is -1)
			-1.0,			// Cuts factor (default is -1.0)
			0,				// Gomory cuts tolerance (default is 0)
			-1.0,			// Memory size limit for Branch and Bound tree (default is -1.0)
			0,				// MIP emphasis (default is 0)
			0				// Limit count on the threads used (default is 0)
		);

		if (objct_type == 1) {
			printf("The optimal objective value for the recourse problem (RP) is: $%.2f\n\n", opt_value_RP);
		}
		else {
			printf("The optimal objective value for the recourse problem (RP) is: %.2f%%\n\n", (100.0 * opt_value_RP));
		}

		tie(opt_value_expect, opt_value_determ) = Dam_CPLEX_expected(
			T,				// The total number of periods
			S,				// The number of scenarios in each period
			R_max,			// The maximum possible level the dam can release in a given period
			L_max,			// The maximum level below flood stage the dam level could be lowered in a given period
			L_0,			// The initial/current level below flood stage the dam's water level is in
			cost_lowering,	// The cost per millimetre of dam level lowered
			cost_flood,		// The cost per millimetre of water flood
			cost_import,	// The cost per millimetre of water imported
			xi_volum,		// The vector with possible precipitation volumes for months in the first block
			xi_probs,		// The vector with precipitation probabilities for months in the first block
			objct_type,		// The objective function type (1-minimize expected cost; 2-minimize probability)
			false,			// Whether or not show the display
			true,			// Whether or not to save the result to files
			2,				// Information display intensity (default is 2)
			-1.0,			// Percentage tolerance gap (default is -1.0)
			-1,				// Time limit in seconds (default is -1)
			-1.0,			// Cuts factor (default is -1.0)
			0,				// Gomory cuts tolerance (default is 0)
			-1.0,			// Memory size limit for Branch and Bound tree (default is -1.0)
			0,				// MIP emphasis (default is 0)
			0				// Limit count on the threads used (default is 0)
		);

		if (objct_type == 1) {
			printf("The optimal objective value for the expected value problem is: $%.2f\n", opt_value_expect);
			printf("\tThe optimal objective value for the deterministic problem is: $%.2f\n\n", opt_value_determ);

			printf("The value of stochastic solution (VSS) is: $%.2f\n\n", (opt_value_expect - opt_value_RP));
		}
		else {
			printf("The optimal objective value for the expected value problem is: %.2f%%\n", (opt_value_expect * 100.0));
			printf("\tThe optimal objective value for the deterministic problem is: %.2f%%\n\n", (opt_value_determ * 100.0));

			printf("The value of stochastic solution (VSS) is: %.2f%%\n\n", ((opt_value_expect - opt_value_RP) * 100.0));
		}

		opt_value_evpi = Dam_CPLEX_WS(
			T,				// The total number of periods
			S,				// The number of scenarios in each period
			R_max,			// The maximum possible level the dam can release in a given period
			L_max,			// The maximum level below flood stage the dam level could be lowered in a given period
			L_0,			// The initial/current level below flood stage the dam's water level is in
			cost_lowering,	// The cost per millimetre of dam level lowered
			cost_flood,		// The cost per millimetre of water flood
			cost_import,	// The cost per millimetre of water imported
			xi_volum,		// The vector with possible precipitation volumes for months in the first block
			xi_probs,		// The vector with precipitation probabilities for months in the first block
			objct_type,		// The objective function type (1-minimize expected cost; 2-minimize probability)
			false,			// Whether or not show the display
			true,			// Whether or not to save the result to files
			2,				// Information display intensity (default is 2)
			-1.0,			// Percentage tolerance gap (default is -1.0)
			-1,				// Time limit in seconds (default is -1)
			-1.0,			// Cuts factor (default is -1.0)
			0,				// Gomory cuts tolerance (default is 0)
			-1.0,			// Memory size limit for Branch and Bound tree (default is -1.0)
			0,				// MIP emphasis (default is 0)
			0				// Limit count on the threads used (default is 0)
		);

		if (objct_type == 1) {
			printf("The optimal objective value for the wait-and-see problem (WS) is: $%.2f\n", opt_value_evpi);
			printf("\tThe expected value of perfect information (EVPI) is: $%.2f\n\n", (opt_value_RP - opt_value_evpi));
		}
		else {
			printf("The optimal objective value for the wait-and-see problem (WS) is: %.2f%%\n", (opt_value_evpi * 100.0));
			printf("\tThe expected value of perfect information (EVPI) is: %.2f%%\n\n", ((opt_value_RP - opt_value_evpi) * 100.0));
		}
	}
	catch (const exception& e) {
		printf("\nERROR: Failed to solve problem using CPLEX... \n");
		printf("\t%s\n", e.what());
		exit(1);
	}

	system("PAUSE");
	return 0;
}



#endif