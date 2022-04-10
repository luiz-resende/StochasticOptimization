#ifndef A3_PROBLEM_I_CPP
#define A3_PROBLEM_I_CPP

#include "A0_definit.h"


void solve_problem_1(double	ticket_economy_price) {
	matrix_flt	xi_demand = {			// The matrix with possible demands for different seat types.
		{ 200.0, 50.0, 20.0 },
		{ 175.0, 25.0, 10.0 },
		{ 150.0, 10.0, 5.0 },
	};	

	double prob = (1.0 / 3.0);
	vector_flt	xi_probs = {			// The vector with the probabilities of each scenario.
		prob, prob, prob
	};

	try {
		double probability = accumulate(xi_probs.begin(), xi_probs.end(), 0.0);
		string except = "IOException: ";

		if ((probability != 1.0)) {
			except += "The probabilities of demand for seats do not add up to 1.0...\n\n";
			throw except;
		}
	}
	catch (string e) {
		printf("%s", e.c_str());
		system("PAUSE");
	}


	/* INITIALIZING THE GLOBAL VARIABLES */
	int			S = 3;					// The integer for the number of realizations.
	int			I = 3;					// The number of different seat types.
	double		Max_seats = 200.0;		// The maximum number of seats in a plane.
	vector_flt	A = {					// Seat partition relation between different types.
		1.0 , 1.5 , 2.0
	};
	vector_flt	c_x = {					// The for partitioning seats.
		0.0 , 0.0 , 0.0
	};
	vector_flt	r_y = {					// The revenue per seat of type `i` sold.
		(ticket_economy_price) , (2 * ticket_economy_price) , (3 * ticket_economy_price)
	};


	/* BLOCK FOR SOLVING PROBLEM USING CPLEX MODEL */
	printf("Solving Problem #1 with CPLEX...\n\n");
	double opt_value_RP, opt_value_expect, opt_value_determ, opt_value_evpi;
	try {
		opt_value_RP = Airline_CPLEX_RP(
			S,				// The integer for the number of realizations.
			I,				// The number of different seat types.
			Max_seats,		// The maximum number of seats in a plane.
			A,				// Seat partition relation between different types.
			c_x,			// The for partitioning seats.
			r_y,			// The revenue per seat of type `i` sold.
			xi_demand,		// The matrix with possible demands for different seat types.
			xi_probs,		// The vector with the probabilities of each scenario.
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
		printf("The optimal objective value for the recourse problem (RP) is: $%.2f\n\n", opt_value_RP);

		tie(opt_value_expect, opt_value_determ) = Airline_CPLEX_expected(
			S,				// The integer for the number of realizations.
			I,				// The number of different seat types.
			Max_seats,		// The maximum number of seats in a plane.
			A,				// Seat partition relation between different types.
			c_x,			// The for partitioning seats.
			r_y,			// The revenue per seat of type `i` sold.
			xi_demand,		// The matrix with possible demands for different seat types.
			xi_probs,		// The vector with the probabilities of each scenario.
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
		printf("The optimal objective value for the expected value problem is: $%.2f\n", opt_value_expect);
		printf("\tThe optimal objective value for the deterministic problem is: $%.2f\n\n", opt_value_determ);

		printf("The value of stochastic solution (VSS) is: $%.2f\n\n", abs(opt_value_expect - opt_value_RP));

		opt_value_evpi = Airline_CPLEX_WS(
			S,				// The integer for the number of realizations.
			I,				// The number of different seat types.
			Max_seats,		// The maximum number of seats in a plane.
			A,				// Seat partition relation between different types.
			c_x,			// The for partitioning seats.
			r_y,			// The revenue per seat of type `i` sold.
			xi_demand,		// The matrix with possible demands for different seat types.
			xi_probs,		// The vector with the probabilities of each scenario.
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
		printf("The optimal objective value for the wait-and-see problem (WS) is: $%.2f\n", opt_value_evpi);
		printf("\tThe expected value of perfect information (EVPI) is: $%.2f\n\n", abs(opt_value_evpi - opt_value_RP));
	}
	catch (const exception& e) {
		printf("\nERROR: Failed to solve problem using CPLEX... \n");
		printf("\t%s\n", e.what());
		system("PAUSE");
	}
}



#endif