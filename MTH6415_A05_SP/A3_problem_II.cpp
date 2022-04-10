#ifndef A3_PROBLEM_II_CPP
#define A3_PROBLEM_II_CPP

#include "A0_definit.h"


void solve_problem_2() {
	vector<matrix_flt>	Omega = {		// The matrix with possible scenarios.
		{
			{200 , 1.0},
			{50 , 1.0},
			{20 , 1.0}
		} ,
		{
			{200 , 0.9},
			{50 , 0.9},
			{20 , 0.9}
		} ,
		{
			{200 , 0.8},
			{50 , 0.8},
			{20 , 0.8}
		} ,
		{
			{175 , 1.0},
			{25 , 1.0},
			{10 , 1.0}
		} ,
		{
			{175 , 0.9},
			{25 , 0.9},
			{10 , 0.9},
		} ,
		{
			{175 , 0.8},
			{25 , 0.8},
			{10 , 0.8}
		} ,
		{
			{150 , 1.0},
			{10 , 1.0},
			{5 , 1.0},
		} ,
		{
			{150 , 0.9},
			{10 , 0.9},
			{5 , 0.9},
		} ,
		{
			{150 , 0.8},
			{10 , 0.8},
			{5 , 0.8}
		}
	};	

	double prob_1 = ((1.0 / 3.0) * (1.0 / 2.0));
	double prob_2 = ((1.0 / 3.0) * (1.0 / 4.0));
	vector_flt	Omega_probs = {			// The vector with the probabilities of each scenario.
		prob_1 , prob_2 , prob_2 , prob_1 , prob_2 , prob_2 , prob_1 , prob_2 , prob_2
	};

	try {
		double probability = accumulate(Omega_probs.begin(), Omega_probs.end(), 0.0);
		string except = "IOException: ";

		if ((probability != 1.0)) {
			except += "The probabilities of scenarios do not add up to 1.0...\n\n";
			throw except;
		}
	}
	catch (string e) {
		printf("%s", e.c_str());
		system("PAUSE");
	}


	/* INITIALIZING THE GLOBAL VARIABLES */
	int			I = 3;					// The number of different seat types.
	double		b = 200;				// The maximum number of Economy seats in a plane.
	vector_flt	A = {					// Seat partition relation between different types.
		1.0 , 1.5 , 2.0
	};
	vector_flt	c_part = {				// The cost for partitioning seats.
		0.0 , 0.0 , 0.0
	};
	vector_flt	c_deny = {				// The cost for denying reservations.
		1.5 , 3.0 , 6.0
	};
	vector_flt	r_book = {				// The revenue per seat of type `i` booked.
		0.0 , 0.0 , 0.0
	};
	vector_flt	r_show = {				// The revenue per seat of type `i` actually occupied.
		1.0 , 2.0 , 4.0
	};
	vector_flt	x = {					// The fixed partition of seats
		140.0 , 24.0 , 12.0
	};


	/* BLOCK FOR SOLVING PROBLEM USING CPLEX MODEL */
	printf("Solving Problem #2 with CPLEX...\n\n");
	double opt_value_RP_3S, opt_value_expect_3S;
	try {
		opt_value_RP_3S = Airline_CPLEX_RP_3S(
			I,				// The number of different seat types.
			b,				// The maximum number of seats in a plane.
			A,				// Seat partition relation between different types.
			Omega,			// The matrix with possible scenarios.
			Omega_probs,	// The vector with the probabilities of each scenario.
			c_part,			// The cost for partitioning seats.
			c_deny,			// The cost per seat of type `i` whose reservation is denied.
			r_book,			// The revenue per seat of type `i` booked.
			r_show,			// The revenue per seat of type `i` actually filled.
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
		printf("The optimal objective value for the 3-stage recourse problem is: $%.2f\n\n", opt_value_RP_3S);

		opt_value_expect_3S = Airline_CPLEX_RP_3S_fixed(
			I,				// The number of different seat types.
			Omega,			// The matrix with possible scenarios.
			Omega_probs,	// The vector with the probabilities of each scenario.
			x,				// The seat partitioning.
			c_part,			// The cost for partitioning seats.
			c_deny,			// The cost per seat of type `i` whose reservation is denied.
			r_book,			// The revenue per seat of type `i` booked.
			r_show,			// The revenue per seat of type `i` actually filled.
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
		printf("The optimal objective value for the 3-stage recourse problem with fixed partition is: $%.2f\n\n", opt_value_expect_3S);
	}
	catch (const exception& e) {
		printf("\nERROR: Failed to solve problem using CPLEX... \n");
		printf("\t%s\n", e.what());
		system("PAUSE");
	}
}



#endif