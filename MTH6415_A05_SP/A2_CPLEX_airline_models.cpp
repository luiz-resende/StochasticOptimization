#ifndef A2_CPLEX_AIRLINE_MODELS_CPP
#define A2_CPLEX_AIRLINE_MODELS_CPP

#include "A0_definit.h"


/**
* Function solves the Northam Airlines two-stage stochastic problem (Exercise 1, pp. 50, Birge & Louveaux, 2011) using CPLEX library.
* 
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
* 
* @param S	The integer for the number of realizations.
* @param I	The number of different seat types.
* @param Max_seats	The maximum number of seats in a plane.
* @param A	Seat partition relation between different types.
* @param c_x	The for partitioning seats.
* @param r_y	The revenue per seat of type `i` sold.
* @param xi_demand	The matrix with possible demands for different seat types.
* @param xi_probs	The vector with the probabilities of each scenario.
* @param show_cplex_display	Whether or not to show CPLEX solution display. The default is `true`.
* @param save_results Whether or not to save the results to files. The default is `true`.
* @param display_intensity	Sets different levels of output display. The default is `2`. (0-No display , 1-Display integer feasible solutions , 2-Default , 3-Number cuts , 4-More information , 5-All information)
* @param gap_tolerance	Gap tolerance for stoping criteria. The default is `-1.0` to keep the CPLEX default.
* @param time_tolerance	The maximum limit time, in seconds, for running CPLEX. The defautl is `-1.0` to keep CPLEX default value.
* @param cuts_factor	Limits the number of cuts that can be added. The default is `-1.0`. (-1.0-CPLEX dynamically adjusts the limit, [0.0,1.0]-For values between the range [0.0, 1.0] CPLEX generates no cuts, >1.0-CPLEX limits the number of rows in the model with cuts added)
* @param gomory_cuts_tolerance	Limits the number of cuts added by CPLEX. The default is `0`. (-1-Do not limit, 0-Default, 1-Moderate, 2-Aggressively)
* @param tree_size_limit	Absolute upper limit on the size (in megabytes, uncompressed) of the branch-and-cut tree. The default is `-1.0` to keep CPLEX default value.
* @param mip_emphasis	Controls trade-offs between speed, feasibility, optimality, and moving bounds in MIP. The default is `0`. (0-Balance optimality and feasibility, 1-Feasibility over optimality, 2-Optimality over feasibility, 3-Moving best bound, 4-Finding hidden feasible solutions)
* @param thread_limit	Sets the default maximal number of parallel threads that will be invoked by any CPLEX parallel optimizer. The default is `0`. (0-Automatic, 1-Single thread, int N>1-Uses up to N threads limited by the available processors)
* @return objective_fnc_value	The objective function value for the RP problem.
*/
double Airline_CPLEX_RP(
	int			S,
	int			I,
	double		Max_seats,
	vector_flt	A,
	vector_flt	c_x,
	vector_flt	r_y,
	matrix_flt	xi_demand,
	vector_flt	xi_probs,
	bool show_cplex_display,//=true,
	bool save_results,//=false
	int display_intensity,//=2,
	double gap_tolerance,//=-1.0,
	int time_tolerance,//=-1.0,
	double cuts_factor,//=-1.0,
	int gomory_cuts_tolerance,//=0,
	double tree_size_limit,//=-1.0,
	int mip_emphasis,//=0,
	int thread_limit//=0
	) {

	/* CREATE ENVIRONMENT AND MODELLING OBJECTS */
	IloEnv		env;							// CPLEX environment
	IloModel	model(env, "RP");				// CPLEX model to store problem and naming it "RP"
	IloNum		lb_x = 0.0;						// x variables lower bound
	IloNum		ub_x = IloInfinity;				// x variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		rhs_partition_ct = 0.0;			// Right-hand-side of partition capacity constraints
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// The CPLEX status of the solution found for the problem


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	IloNumVarArray				x(env, I);					// Variables for partition decision of seat type `i` (first-stage)
	IloArray<IloNumVarArray>	y(env, S);					// Decision variables for seats of type `i` sold under scenario `s`
	
	for (int i = 0; i < I; i++) {
		x[i] = IloNumVar(env, lb_x, ub_x, ILOINT);			// Defining x variables
	}

	for (int s = 0; s < S; s++) {
		y[s] = IloNumVarArray(env, I);							// Assigning a array with the number of types
		
		for (int i = 0; i < I; i++) {
			y[s][i] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);		// Defining y variables
		}
	}


	/* FIRST-STAGE CONSTRAINTS */
	IloExpr first_stage(env);

	for (int i = 0; i < I; i++) {
		model.add(x[i] >= rhs_nonneg_ct);						// Adding 'x_t >= 0.0' constraint to model
		
		first_stage += (A[i] * x[i]);
	}
	model.add((first_stage - Max_seats) <= rhs_partition_ct);	// Adding 'aTx <= b' constraint to model
	first_stage.end();


	/* SECOND-STAGE CONSTRAINTS */

	for (int s = 0; s < S; s++) {

		for (int i = 0; i < I; i++) {

			model.add(y[s][i] >= rhs_nonneg_ct);						// Adding 'y_si >= 0.0' constraint to model

			model.add((y[s][i] - x[i]) <= rhs_nonneg_ct);				// Adding 'y_si <= x_i' constraint to model

			model.add((y[s][i] - xi_demand[s][i]) <= rhs_nonneg_ct);	// Adding 'y_si <= xi_si' constraint to model
		}
	}

	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);										// Expression for the objective function first-stage.
	IloExpr	obj_fnc_2(env);										// Expression for the objective function second-stage.


	for (int i = 0; i < I; i++) {
		obj_fnc_1 -= (c_x[i] * x[i]);						// Partition cost
	}

	for (int s = 0; s < S; s++) {
		IloExpr	temp_expr(env);

		for (int i = 0; i < I; i++) {
			temp_expr += (r_y[i] * y[s][i]);
		}
		obj_fnc_2 += (xi_probs[s] * temp_expr);				// Expected profit from selling different seats under scenario `s`
		temp_expr.end();
	}

	IloObjective objective_function = IloMaximize(		// Defining an objective function to be maximized
		env,
		(obj_fnc_1 + obj_fnc_2),
		"Objective_Function");
	model.add(objective_function);						// Adding objective function expression to the model
	obj_fnc_1.end();									// Ending objective fnc expression 1 to prevent memory leaks.
	obj_fnc_2.end();									// Ending objective fnc expression 2 to prevent memory leaks.


	/* CREATE SOLUTION OBJECT, WHICH WILL BE RESPONSIBLE TO CALL SOLVING FUNCTION AND RESULTS */
	IloCplex mdl_solution(model);

	/* CPLEX PARAMETERS TUNING */
	if (!show_cplex_display) {		// Turn-off CPLEX logging screen
		mdl_solution.setOut(env.getNullStream());
	}
	if (show_cplex_display) {
		mdl_solution.setParam(IloCplex::Param::MIP::Display, display_intensity);
	}
	if (time_tolerance >= 0) {
		mdl_solution.setParam(IloCplex::Param::TimeLimit, time_tolerance);
	}
	if (gap_tolerance >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, gap_tolerance);
	}
	if (cuts_factor >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::CutsFactor, cuts_factor);
	}
	if (gomory_cuts_tolerance != 0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Cuts::Gomory, gomory_cuts_tolerance);
	}
	if (tree_size_limit >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::TreeMemory, tree_size_limit);
	}
	if (mip_emphasis > 0) {
		mdl_solution.setParam(IloCplex::Param::Emphasis::MIP, mip_emphasis);
	}
	if (thread_limit > 0) {
		mdl_solution.setParam(IloCplex::Param::Threads, thread_limit);
	}
	//mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1e-10);	// e-optimal solution (absolute value)


	/* CALLING SOLVING FUNCTION */
	mdl_solution.solve();


	/* GETTING SOLUTION STATUS */
	solution_status = mdl_solution.getStatus();
	if (solution_status != 3) {
		objective_fnc_value = (double)mdl_solution.getObjValue();
	}
	if (show_cplex_display) {
		env.out() << "\n###########################################################\n";
		env.out() << "CPLEX Solution status = " << mdl_solution.getStatus() << endl;
		if (solution_status != 3) {
			env.out() << "SP objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}

	/* RETRIEVING SOLUTION */
	vector_flt	x_solution(I, 0.0);
	matrix_flt	y_solution(S, vector_flt(I, 0.0));

	if (solution_status != 3) {
		for (int i = 0; i < I; i++) {
			x_solution[i] = (double)mdl_solution.getValue(x[i]);

			for (int s = 0; s < S; s++) {
				y_solution[s][i] = (double)mdl_solution.getValue(y[s][i]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_vector(x_solution, false, true, "solution_x_i_CPLEX_RP.txt", "Outputs", "\t");
			print_matrix(y_solution, false, true, "solution_y_si_CPLEX_RP.txt", "Outputs", "\t");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	x.end();
	y.end();
	mdl_solution.end();
	env.end();
	//model.end();
	x_solution.clear();
	x_solution.shrink_to_fit();
	y_solution.clear();
	y_solution.shrink_to_fit();

	return objective_fnc_value;
}


/**
* Function solves the Northam Airlines expected value solution problem (Exercise 1, pp. 50, Birge & Louveaux, 2011) using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param S	The integer for the number of realizations.
* @param I	The number of different seat types.
* @param Max_seats	The maximum number of seats in a plane.
* @param A	Seat partition relation between different types.
* @param c_x	The for partitioning seats.
* @param r_y	The revenue per seat of type `i` sold.
* @param xi_demand	The matrix with possible demands for different seat types.
* @param xi_probs	The vector with the probabilities of each scenario.
* @param show_cplex_display	Whether or not to show CPLEX solution display. The default is `true`.
* @param save_results Whether or not to save the results to files. The default is `true`.
* @param display_intensity	Sets different levels of output display. The default is `2`. (0-No display , 1-Display integer feasible solutions , 2-Default , 3-Number cuts , 4-More information , 5-All information)
* @param gap_tolerance	Gap tolerance for stoping criteria. The default is `-1.0` to keep the CPLEX default.
* @param time_tolerance	The maximum limit time, in seconds, for running CPLEX. The defautl is `-1.0` to keep CPLEX default value.
* @param cuts_factor	Limits the number of cuts that can be added. The default is `-1.0`. (-1.0-CPLEX dynamically adjusts the limit, [0.0,1.0]-For values between the range [0.0, 1.0] CPLEX generates no cuts, >1.0-CPLEX limits the number of rows in the model with cuts added)
* @param gomory_cuts_tolerance	Limits the number of cuts added by CPLEX. The default is `0`. (-1-Do not limit, 0-Default, 1-Moderate, 2-Aggressively)
* @param tree_size_limit	Absolute upper limit on the size (in megabytes, uncompressed) of the branch-and-cut tree. The default is `-1.0` to keep CPLEX default value.
* @param mip_emphasis	Controls trade-offs between speed, feasibility, optimality, and moving bounds in MIP. The default is `0`. (0-Balance optimality and feasibility, 1-Feasibility over optimality, 2-Optimality over feasibility, 3-Moving best bound, 4-Finding hidden feasible solutions)
* @param thread_limit	Sets the default maximal number of parallel threads that will be invoked by any CPLEX parallel optimizer. The default is `0`. (0-Automatic, 1-Single thread, int N>1-Uses up to N threads limited by the available processors)
* @param (objective_fnc_value, objective_fnc_determ)	Tuple with the objective function values for the expected problem and deterministic problem.
*/
tuple<double, double> Airline_CPLEX_expected(
	int			S,
	int			I,
	double		Max_seats,
	vector_flt	A,
	vector_flt	c_x,
	vector_flt	r_y,
	matrix_flt	xi_demand,
	vector_flt	xi_probs,
	bool show_cplex_display,
	bool save_results,
	int display_intensity,
	double gap_tolerance,
	int time_tolerance,
	double cuts_factor,
	int gomory_cuts_tolerance,
	double tree_size_limit,
	int mip_emphasis,
	int thread_limit
) {

	/* CREATE ENVIRONMENT AND MODELLING OBJECTS */
	IloEnv		env;							// CPLEX environment
	IloModel	model(env, "RP_expected");		// CPLEX model to store problem and naming it "RP_expected"
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	double		objective_fnc_deter = -999.0;	// Variable to store the objective function value of the deterministic model
	int			solution_status;				// The CPLEX status of the solution found for the problem


	/* CREATING VECTOR OF SIZE T WITH MEAN PRECIPITATION/EVAPORATION VOLUMES PER MONTH */
	vector_flt xi_mean; // = { 175.0 , 28.0 , 12.0 };
	for (int i = 0; i < I; i++) {
		double sum = 0.0;

		for (int s = 0; s < S; s++) {
			sum += (xi_demand[s][i] * xi_probs[s]);
		}
		xi_mean.push_back(sum);
	}
	xi_mean = round_values(xi_mean);


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	vector_flt					x;								// The release decisions for mean realiations 
	IloArray<IloNumVarArray>	y(env, S);						// Variables for seats sold of type `i` under scenario `s`

	for (int s = 0; s < S; s++) {
		y[s] = IloNumVarArray(env, I);							// Assigning a vector with the number of types
		
		for (int i = 0; i < I; i++) {
			y[s][i] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);		// Defining y variables
		}
	}


	/* GETTING THE x VARIABLE VALUES FOR THE PROBLEM USING MEAN PRECIPITATION/EVAPORATION VALUES */
	tie(x, objective_fnc_deter) = Airline_CPLEX_deterministic(S, I, Max_seats, A, c_x, r_y, xi_mean, false, true);		// Saves the deterministic result


	/* SECOND-STAGE CONSTRAINTS */

	for (int s = 0; s < S; s++) {

		for (int i = 0; i < I; i++) {
			model.add(y[s][i] >= rhs_nonneg_ct);						// Adding 'y_si >= 0.0' constraint to model
			model.add((y[s][i] - x[i]) <= rhs_nonneg_ct);				// Adding 'y_si <= x_i' constraint to model
			model.add((y[s][i] - xi_demand[s][i]) <= rhs_nonneg_ct);	// Adding 'y_si <= xi_si' constraint to model
		}
	}

	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);							// Expression for the objective function water level lowering part (first-stage).
	IloExpr	obj_fnc_2(env);							// Expression for the objective function second-stage part.

	for (int i = 0; i < I; i++) {
		obj_fnc_1 -= (c_x[i] * x[i]);						// Cost for partition
	}

	for (int s = 0; s < S; s++) {
		IloExpr	temp_expr(env);

		for (int i = 0; i < I; i++) {
			temp_expr += (r_y[i] * y[s][i]);
		}
		obj_fnc_2 += (xi_probs[s] * temp_expr);				// Expected revenue from selling seats
		temp_expr.end();
	}

	IloObjective objective_function = IloMaximize(env,
		(obj_fnc_1 + obj_fnc_2),
		"Objective_Function");								// Defining an objective function to be maximized
	model.add(objective_function);							// Adding objective function expression to the model
	obj_fnc_1.end();										// Ending objective fnc expression 1 to prevent memory leaks.
	obj_fnc_2.end();										// Ending objective fnc expression 2 to prevent memory leaks.


	/* CREATE SOLUTION OBJECT, WHICH WILL BE RESPONSIBLE TO CALL SOLVING FUNCTION AND RESULTS */
	IloCplex mdl_solution(model);


	/* CPLEX PARAMETERS TUNING */
	if (!show_cplex_display) {		// Turn-off CPLEX logging screen
		mdl_solution.setOut(env.getNullStream());
	}
	if (show_cplex_display) {
		mdl_solution.setParam(IloCplex::Param::MIP::Display, display_intensity);
	}
	if (time_tolerance >= 0) {
		mdl_solution.setParam(IloCplex::Param::TimeLimit, time_tolerance);
	}
	if (gap_tolerance >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, gap_tolerance);
	}
	if (cuts_factor >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::CutsFactor, cuts_factor);
	}
	if (gomory_cuts_tolerance != 0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Cuts::Gomory, gomory_cuts_tolerance);
	}
	if (tree_size_limit >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::TreeMemory, tree_size_limit);
	}
	if (mip_emphasis > 0) {
		mdl_solution.setParam(IloCplex::Param::Emphasis::MIP, mip_emphasis);
	}
	if (thread_limit > 0) {
		mdl_solution.setParam(IloCplex::Param::Threads, thread_limit);
	}
	//mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1e-10);	// e-optimal solution (absolute value)


	/* CALLING SOLVING FUNCTION */
	mdl_solution.solve();


	/* GETTING SOLUTION STATUS */
	solution_status = mdl_solution.getStatus();
	if (solution_status != 3) {
		objective_fnc_value = (double)mdl_solution.getObjValue();
	}
	if (show_cplex_display) {
		env.out() << "\n###########################################################\n";
		env.out() << "CPLEX Solution status = " << mdl_solution.getStatus() << endl;
		if (solution_status != 3) {
			env.out() << "SP objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}


	/* RETRIEVING SOLUTION */
	matrix_flt	y_solution(S, vector_flt(I, 0.0));

	if (solution_status != 3) {

		for (int s = 0; s < S; s++) {

			for (int i = 0; i < I; i++) {
				y_solution[s][i] = (double)mdl_solution.getValue(y[s][i]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_vector(x, false, true, "solution_x_i_CPLEX_expected.txt", "Outputs", "\t");
			print_matrix(y_solution, false, true, "solution_y_si_CPLEX_expected.txt", "Outputs", "\t");
		}
	}


	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	y.end();
	mdl_solution.end();
	env.end();
	//model.end();
	x.clear();
	x.shrink_to_fit();
	y_solution.clear();
	y_solution.shrink_to_fit();

	return make_tuple(objective_fnc_value, objective_fnc_deter);
}


/**
* Function solves the Northam Airlines mean precipitation problem (Exercise 1, pp. 50, Birge & Louveaux, 2011) using CPLEX library to get solution
* for calculating the deterministic problem.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param S	The integer for the number of realizations.
* @param I	The number of different seat types.
* @param Max_seats	The maximum number of seats in a plane.
* @param A	Seat partition relation between different types.
* @param c_x	The for partitioning seats.
* @param r_y	The revenue per seat of type `i` sold.
* @param xi_mean	The vector with the mean demands for the seats.
* @param show_cplex_display	Whether or not to show CPLEX solution display. The default is `true`.
* @param save_results Whether or not to save the results to files. The default is `true`.
* @param display_intensity	Sets different levels of output display. The default is `2`. (0-No display , 1-Display integer feasible solutions , 2-Default , 3-Number cuts , 4-More information , 5-All information)
* @param gap_tolerance	Gap tolerance for stoping criteria. The default is `-1.0` to keep the CPLEX default.
* @param time_tolerance	The maximum limit time, in seconds, for running CPLEX. The defautl is `-1.0` to keep CPLEX default value.
* @param cuts_factor	Limits the number of cuts that can be added. The default is `-1.0`. (-1.0-CPLEX dynamically adjusts the limit, [0.0,1.0]-For values between the range [0.0, 1.0] CPLEX generates no cuts, >1.0-CPLEX limits the number of rows in the model with cuts added)
* @param gomory_cuts_tolerance	Limits the number of cuts added by CPLEX. The default is `0`. (-1-Do not limit, 0-Default, 1-Moderate, 2-Aggressively)
* @param tree_size_limit	Absolute upper limit on the size (in megabytes, uncompressed) of the branch-and-cut tree. The default is `-1.0` to keep CPLEX default value.
* @param mip_emphasis	Controls trade-offs between speed, feasibility, optimality, and moving bounds in MIP. The default is `0`. (0-Balance optimality and feasibility, 1-Feasibility over optimality, 2-Optimality over feasibility, 3-Moving best bound, 4-Finding hidden feasible solutions)
* @param thread_limit	Sets the default maximal number of parallel threads that will be invoked by any CPLEX parallel optimizer. The default is `0`. (0-Automatic, 1-Single thread, int N>1-Uses up to N threads limited by the available processors)
* @param (x_solution, objective_fnc_value)	Tuple with the solution vector `x` and the objective function value for the deterministic problem.
*/
tuple<vector_flt, double> Airline_CPLEX_deterministic(
	int			S,
	int			I,
	double		Max_seats,
	vector_flt	A,
	vector_flt	c_x,
	vector_flt	r_y,
	vector_flt	xi_mean,
	bool show_cplex_display,
	bool save_results,
	int display_intensity,
	double gap_tolerance,
	int time_tolerance,
	double cuts_factor,
	int gomory_cuts_tolerance,
	double tree_size_limit,
	int mip_emphasis,
	int thread_limit
) {

	/* CREATE ENVIRONMENT AND MODELLING OBJECTS */
	IloEnv		env;							// CPLEX environment
	IloModel	model(env, "Deterministic");	// CPLEX model to store problem and naming it "Deterministic"
	IloNum		lb_x = 0.0;						// x variables lower bound
	IloNum		ub_x = IloInfinity;				// x variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		rhs_partition_ct = 0.0;			// Right-hand-side of partition capacity
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// The status of the solution for the problem


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	IloNumVarArray			x(env, I);				// First-stage decision variables for partition
	IloNumVarArray			y(env, I);				// Second-stage varaibles - seats sold

	for (int i = 0; i < I; i++) {

		x[i] = IloNumVar(env, lb_x, ub_x, ILOINT);				// Defining x variables
		y[i] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);				// Defining y variables
	}


	/* PARTITION CAPACITY CONSTRAINTS */
	IloExpr	partition(env);

	for (int i = 0; i < I; i++) {
		model.add(x[i] >= rhs_nonneg_ct);						// Adding 'x_i >= 0.0' constraint to model

		partition += (A[i] * x[i]);
	}
	model.add((partition - Max_seats) <= rhs_partition_ct);		// Adding 'aTx <= b' constraint to model
	partition.end();


	/* SEATS SELLING CAPACITY CONSTRAINTS */

	for (int i = 0; i < I; i++) {
		model.add(y[i] >= rhs_nonneg_ct);					// Adding 'y_i >= 0.0' constraint to model

		model.add((y[i] - x[i]) <= rhs_nonneg_ct);			// Adding 'y_i <= x_i' constraint to model
		model.add((y[i] - xi_mean[i]) <= rhs_nonneg_ct);	// Adding 'y_i <= xi_i' constraint to model
	}


	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);									// Expression for cost of partitioning.
	IloExpr	obj_fnc_2(env);									// Expression for the revenue from selling.

	for (int i = 0; i < I; i++) {

		obj_fnc_1 -= (c_x[i] * x[i]);						// Cost for partitioning
		obj_fnc_2 += (r_y[i] * y[i]);						// Revenue from selling
	}

	IloObjective objective_function = IloMaximize(				// Defining an objective function to be maximized
		env,
		(obj_fnc_1 + obj_fnc_2),
		"Objective_Function");
	model.add(objective_function);								// Adding objective function expression to the model
	obj_fnc_1.end();											// Ending objective fnc expression 1 to prevent memory leaks.
	obj_fnc_2.end();											// Ending objective fnc expression 2 to prevent memory leaks.


	/* CREATE SOLUTION OBJECT, WHICH WILL BE RESPONSIBLE TO CALL SOLVING FUNCTION AND RESULTS */
	IloCplex mdl_solution(model);


	/* CPLEX PARAMETERS TUNING */
	if (!show_cplex_display) {		// Turn-off CPLEX logging screen
		mdl_solution.setOut(env.getNullStream());
	}
	if (show_cplex_display) {
		mdl_solution.setParam(IloCplex::Param::MIP::Display, display_intensity);
	}
	if (time_tolerance >= 0) {
		mdl_solution.setParam(IloCplex::Param::TimeLimit, time_tolerance);
	}
	if (gap_tolerance >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, gap_tolerance);
	}
	if (cuts_factor >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::CutsFactor, cuts_factor);
	}
	if (gomory_cuts_tolerance != 0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Cuts::Gomory, gomory_cuts_tolerance);
	}
	if (tree_size_limit >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::TreeMemory, tree_size_limit);
	}
	if (mip_emphasis > 0) {
		mdl_solution.setParam(IloCplex::Param::Emphasis::MIP, mip_emphasis);
	}
	if (thread_limit > 0) {
		mdl_solution.setParam(IloCplex::Param::Threads, thread_limit);
	}
	//mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1e-10);	// e-optimal solution (absolute value)


	/* CALLING SOLVING FUNCTION */
	mdl_solution.solve();


	/* GETTING SOLUTION STATUS */
	solution_status = mdl_solution.getStatus();
	if (solution_status != 3) {
		objective_fnc_value = (double)mdl_solution.getObjValue();
	}
	if (show_cplex_display) {
		env.out() << "\n###########################################################\n";
		env.out() << "CPLEX Solution status = " << mdl_solution.getStatus() << endl;
		if (solution_status != 3) {																// If the solution is not infeasible
			env.out() << "Deterministic problem objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}


	/* RETRIEVING SOLUTION */
	vector_flt	x_solution(I, -999.0);
	vector_flt	y_solution(I, -999.0);

	if (solution_status != 3) {											// If problem is not infeasible

		for (int i = 0; i < I; i++) {
			x_solution[i] = (double)mdl_solution.getValue(x[i]);
			y_solution[i] = (double)mdl_solution.getValue(y[i]);
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_vector(xi_mean, false, true, "solution_mean_xi_i_CPLEX_deterministic.txt", "Outputs", "\t");
			print_vector(x_solution, false, true, "solution_x_i_CPLEX_deterministic.txt", "Outputs", "\t");
			print_vector(y_solution, false, true, "solution_y_i_CPLEX_deterministic.txt", "Outputs", "\t");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	x.end();
	y.end();
	mdl_solution.end();
	env.end();
	y_solution.clear();
	y_solution.shrink_to_fit();
	//model.end();

	return make_tuple(x_solution, objective_fnc_value);
}


/**
* Function solves the Northam Airlines wait-and-see stochastic problem (Exercise 1, pp. 50, Birge & Louveaux, 2011) using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param S	The integer for the number of realizations.
* @param I	The number of different seat types.
* @param Max_seats	The maximum number of seats in a plane.
* @param A	Seat partition relation between different types.
* @param c_x	The for partitioning seats.
* @param r_y	The revenue per seat of type `i` sold.
* @param xi_demand	The matrix with possible demands for different seat types.
* @param xi_probs	The vector with the probabilities of each scenario.
* @param show_cplex_display	Whether or not to show CPLEX solution display. The default is `true`.
* @param save_results Whether or not to save the results to files. The default is `true`.
* @param display_intensity	Sets different levels of output display. The default is `2`. (0-No display , 1-Display integer feasible solutions , 2-Default , 3-Number cuts , 4-More information , 5-All information)
* @param gap_tolerance	Gap tolerance for stoping criteria. The default is `-1.0` to keep the CPLEX default.
* @param time_tolerance	The maximum limit time, in seconds, for running CPLEX. The defautl is `-1.0` to keep CPLEX default value.
* @param cuts_factor	Limits the number of cuts that can be added. The default is `-1.0`. (-1.0-CPLEX dynamically adjusts the limit, [0.0,1.0]-For values between the range [0.0, 1.0] CPLEX generates no cuts, >1.0-CPLEX limits the number of rows in the model with cuts added)
* @param gomory_cuts_tolerance	Limits the number of cuts added by CPLEX. The default is `0`. (-1-Do not limit, 0-Default, 1-Moderate, 2-Aggressively)
* @param tree_size_limit	Absolute upper limit on the size (in megabytes, uncompressed) of the branch-and-cut tree. The default is `-1.0` to keep CPLEX default value.
* @param mip_emphasis	Controls trade-offs between speed, feasibility, optimality, and moving bounds in MIP. The default is `0`. (0-Balance optimality and feasibility, 1-Feasibility over optimality, 2-Optimality over feasibility, 3-Moving best bound, 4-Finding hidden feasible solutions)
* @param thread_limit	Sets the default maximal number of parallel threads that will be invoked by any CPLEX parallel optimizer. The default is `0`. (0-Automatic, 1-Single thread, int N>1-Uses up to N threads limited by the available processors)
* @return objective_fnc_value	The objective function value for the RP problem.
*/
double Airline_CPLEX_WS(
	int			S,
	int			I,
	double		Max_seats,
	vector_flt	A,
	vector_flt	c_x,
	vector_flt	r_y,
	matrix_flt	xi_demand,
	vector_flt	xi_probs,
	bool show_cplex_display,
	bool save_results,
	int display_intensity,
	double gap_tolerance,
	int time_tolerance,
	double cuts_factor,
	int gomory_cuts_tolerance,
	double tree_size_limit,
	int mip_emphasis,
	int thread_limit
) {

	/* CREATE ENVIRONMENT AND MODELLING OBJECTS */
	IloEnv		env;							// CPLEX environment
	IloModel	model(env, "WS");				// CPLEX model to store problem and naming it "WS"
	IloNum		lb_x = 0.0;						// x variables lower bound
	IloNum		ub_x = IloInfinity;				// x variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		rhs_patition_ct = 0.0;			// Right-hand-side of partition capacity constraints
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// The CPLEX status of the solution found for the problem


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	IloArray<IloNumVarArray>	x(env, S);					// Variables for seat partition of type `i` under scenario `s`
	IloArray<IloNumVarArray>	y(env, S);					// Decision variables for seats sold of type `i` under scenario `s`

	for (int s = 0; s < S; s++) {
		x[s] = IloNumVarArray(env, I);						// Assigning a array with the number of types
		y[s] = IloNumVarArray(env, I);						// Assigning a array with the number of types

		for (int i = 0; i < I; i++) {
			x[s][i] = IloNumVar(env, lb_x, ub_x, ILOINT);		// Defining x variables
			y[s][i] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);		// Defining y variables
		}
	}


	/* LINKING LEVEL OF THE DAM TO y AND w VARIABLES */

	for (int s = 0; s < S; s++) {
		IloExpr	partition_s(env);

		for (int i = 0; i < I; i++) {

			model.add(x[s][i] >= rhs_nonneg_ct);						// Adding 'x_i >= 0.0' constraint to model
			model.add(y[s][i] >= rhs_nonneg_ct);						// Adding 'y_i >= 0.0' constraint to model

			model.add((y[s][i] - x[s][i]) <= rhs_nonneg_ct);			// Adding 'y_si <= x_si' constraint to model
			model.add((y[s][i] - xi_demand[s][i]) <= rhs_nonneg_ct);	// Adding 'y_si <= x_si' constraint to model

			partition_s += (A[i] * x[s][i]);
		}
		model.add((partition_s - Max_seats) <= rhs_patition_ct);		// Adding 'aTx_s <= b'
		partition_s.end();
	}

	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc(env);											// Expression for the objective function.

	for (int s = 0; s < S; s++) {
		IloExpr	temp_expr(env);

		for (int i = 0; i < I; i++) {
			temp_expr += (( - (c_x[i] * x[s][i])) + (r_y[i] * y[s][i]));
		}
		obj_fnc += (xi_probs[s] * temp_expr);			// Expected profit over types `i` under scenario `s`
		temp_expr.end();
	}

	IloObjective objective_function = IloMaximize(		// Defining an objective function to be maximized
		env,
		obj_fnc,
		"Objective_Function");
	model.add(objective_function);						// Adding objective function expression to the model
	obj_fnc.end();										// Ending objective fnc expression to prevent memory leaks.


	/* CREATE SOLUTION OBJECT, WHICH WILL BE RESPONSIBLE TO CALL SOLVING FUNCTION AND RESULTS */
	IloCplex mdl_solution(model);

	/* CPLEX PARAMETERS TUNING */
	if (!show_cplex_display) {		// Turn-off CPLEX logging screen
		mdl_solution.setOut(env.getNullStream());
	}
	if (show_cplex_display) {
		mdl_solution.setParam(IloCplex::Param::MIP::Display, display_intensity);
	}
	if (time_tolerance >= 0) {
		mdl_solution.setParam(IloCplex::Param::TimeLimit, time_tolerance);
	}
	if (gap_tolerance >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, gap_tolerance);
	}
	if (cuts_factor >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::CutsFactor, cuts_factor);
	}
	if (gomory_cuts_tolerance != 0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Cuts::Gomory, gomory_cuts_tolerance);
	}
	if (tree_size_limit >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::TreeMemory, tree_size_limit);
	}
	if (mip_emphasis > 0) {
		mdl_solution.setParam(IloCplex::Param::Emphasis::MIP, mip_emphasis);
	}
	if (thread_limit > 0) {
		mdl_solution.setParam(IloCplex::Param::Threads, thread_limit);
	}
	//mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1e-10);	// e-optimal solution (absolute value)


	/* CALLING SOLVING FUNCTION */
	mdl_solution.solve();


	/* GETTING SOLUTION STATUS */
	solution_status = mdl_solution.getStatus();
	if (solution_status != 3) {
		objective_fnc_value = (double)mdl_solution.getObjValue();
	}
	if (show_cplex_display) {
		env.out() << "\n###########################################################\n";
		env.out() << "CPLEX Solution status = " << mdl_solution.getStatus() << endl;
		if (solution_status != 3) {
			env.out() << "SP objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}

	/* RETRIEVING SOLUTION */
	matrix_flt	x_solution(S, vector_flt(I, 0.0));
	matrix_flt	y_solution(S, vector_flt(I, 0.0));

	if (solution_status != 3) {

		for (int i = 0; i < I; i++) {

			for (int s = 0; s < S; s++) {
				x_solution[s][i] = (double)mdl_solution.getValue(x[s][i]);
				y_solution[s][i] = (double)mdl_solution.getValue(y[s][i]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_matrix(x_solution, false, true, "solution_x_si_CPLEX_EVPI.txt", "Outputs", "\t");
			print_matrix(y_solution, false, true, "solution_y_si_CPLEX_EVPI.txt", "Outputs", "\t");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	x.end();
	y.end();
	mdl_solution.end();
	env.end();
	//model.end();
	x_solution.clear();
	x_solution.shrink_to_fit();
	y_solution.clear();
	y_solution.shrink_to_fit();

	return objective_fnc_value;
}


/**
* Function solves the Northam Airlines three-stage stochastic problem (Exercise 1, pp. 83, Birge & Louveaux, 2011) using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param I	The number of different seat types.
* @param b	The maximum number of economy seats in the plane.
* @param a_part	The vector with the seat partition proportions.
* @param Omega	The demand and show-up scenarios for each seat type.
* @param Omega_probs	The probability of each scenario `s`.
* @param c_part	The vector with cost for partitioning seats.
* @param c_deny	The vector with the costs for denying reservations.
* @param r_book	The vector with revenue from booking seats of type `i`.
* @param r_show	The vector with revenue per reservation show-up of type `i`.
* @param show_cplex_display	Whether or not to show CPLEX solution display. The default is `true`.
* @param save_results Whether or not to save the results to files. The default is `true`.
* @param display_intensity	Sets different levels of output display. The default is `2`. (0-No display , 1-Display integer feasible solutions , 2-Default , 3-Number cuts , 4-More information , 5-All information)
* @param gap_tolerance	Gap tolerance for stoping criteria. The default is `-1.0` to keep the CPLEX default.
* @param time_tolerance	The maximum limit time, in seconds, for running CPLEX. The defautl is `-1.0` to keep CPLEX default value.
* @param cuts_factor	Limits the number of cuts that can be added. The default is `-1.0`. (-1.0-CPLEX dynamically adjusts the limit, [0.0,1.0]-For values between the range [0.0, 1.0] CPLEX generates no cuts, >1.0-CPLEX limits the number of rows in the model with cuts added)
* @param gomory_cuts_tolerance	Limits the number of cuts added by CPLEX. The default is `0`. (-1-Do not limit, 0-Default, 1-Moderate, 2-Aggressively)
* @param tree_size_limit	Absolute upper limit on the size (in megabytes, uncompressed) of the branch-and-cut tree. The default is `-1.0` to keep CPLEX default value.
* @param mip_emphasis	Controls trade-offs between speed, feasibility, optimality, and moving bounds in MIP. The default is `0`. (0-Balance optimality and feasibility, 1-Feasibility over optimality, 2-Optimality over feasibility, 3-Moving best bound, 4-Finding hidden feasible solutions)
* @param thread_limit	Sets the default maximal number of parallel threads that will be invoked by any CPLEX parallel optimizer. The default is `0`. (0-Automatic, 1-Single thread, int N>1-Uses up to N threads limited by the available processors)
* @return objective_fnc_value	The objective function value for the RP problem.
*/
double Airline_CPLEX_RP_3S(
	int					I,
	int					b,
	vector_flt			a_part,
	vector<matrix_flt>	Omega,
	vector_flt			Omega_probs,
	vector_flt			c_part,
	vector_flt			c_deny,
	vector_flt			r_book,
	vector_flt			r_show,
	bool show_cplex_display,//=true,
	bool save_results,//=false
	int display_intensity,//=2,
	double gap_tolerance,//=-1.0,
	int time_tolerance,//=-1.0,
	double cuts_factor,//=-1.0,
	int gomory_cuts_tolerance,//=0,
	double tree_size_limit,//=-1.0,
	int mip_emphasis,//=0,
	int thread_limit//=0
) {

	/* CREATE ENVIRONMENT AND MODELLING OBJECTS */
	IloEnv		env;							// CPLEX environment
	IloModel	model(env, "3S");				// CPLEX model to store problem and naming it "3S"
	IloNum		lb_x = 0.0;						// x variables lower bound
	IloNum		ub_x = IloInfinity;				// x variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	IloNum		lb_u = 0.0;						// u variables lower bound
	IloNum		ub_u = IloInfinity;				// u variables upper bound
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// The CPLEX status of the solution found for the problem


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	int							O = (int)Omega.size();		// Number of scenarios
	IloNumVarArray				x(env, I);					// Decision variables for seats of type `i` partition
	IloArray<IloNumVarArray>	y(env, O);					// Decision variables for seats of type `i` reserved under scenario `s`
	IloArray<IloNumVarArray>	w(env, O);					// Decision variables for seats of type `i` show under scenario `s`
	IloArray<IloNumVarArray>	u(env, O);					// Decision variables for seats of type `i` denied under scenario `s`

	for (int i = 0; i < I; i++) {
		x[i] = IloNumVar(env, lb_x, ub_x, ILOINT);			// Defining x variables
	}

	for (int s = 0; s < O; s++) {
		y[s] = IloNumVarArray(env, I);							// Assigning a array with the number of types
		w[s] = IloNumVarArray(env, I);							// Assigning a array with the number of types
		u[s] = IloNumVarArray(env, I);							// Assigning a array with the number of types

		for (int i = 0; i < I; i++) {
			y[s][i] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);		// Defining y variables
			w[s][i] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);		// Defining w variables
			u[s][i] = IloNumVar(env, lb_u, ub_u, ILOFLOAT);		// Defining u variables
		}
	}


	/* FIRST-STAGE CONSTRAINTS */
	IloExpr first_stage(env);

	for (int i = 0; i < I; i++) {
		model.add(x[i] >= rhs_nonneg_ct);						// Adding 'x_t >= 0.0' constraint to model

		first_stage += (a_part[i] * x[i]);
	}
	model.add((first_stage - b) <= rhs_nonneg_ct);				// Adding 'aTx <= b' constraint to model
	first_stage.end();


	/* SECOND-STAGE CONSTRAINTS */

	for (int s = 0; s < O; s++) {
		int s_prime = ((int)s - 1);

		for (int i = 0; i < I; i++) {
			model.add(y[s][i] <= Omega[s][i][0]);							// Adding 'y_si <= xi_si' constraint to model

			if ((s > 0) && (Omega[s][i][0] == Omega[s_prime][i][0])) {
				model.add(y[s][i] == y[s_prime][i]);						// Adding 'y_si - y_s'i == 0' -> nonanticipativity
			}

			model.add(y[s][i] >= rhs_nonneg_ct);							// Adding 'y_si >= 0.0' constraint to model
		}
	}


	/* THIRD-STAGE CONSTRAINTS */

	for (int s = 0; s < O; s++) {

		for (int i = 0; i < I; i++) {

			model.add(w[s][i] == ((Omega[s][i][1] * y[s][i]) - u[s][i]));	// Adding 'w_si == (o_si * y_si) - u_si' constraint to model
			model.add(w[s][i] <= x[i]);										// Adding 'w_si <= x_i' constraint to model

			model.add(w[s][i] >= rhs_nonneg_ct);							// Adding 'w_si >= 0.0' constraint to model
			model.add(u[s][i] >= rhs_nonneg_ct);							// Adding 'u_si >= 0.0' constraint to model
		}
	}

	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr obj_fnc_1(env);										// Expression for first-stage part of the objective function.
	IloExpr	obj_fnc_2(env);										// Expression for second-stage part of the objective function.
	IloExpr	obj_fnc_3(env);										// Expression for third-stage part of the objective function.

	for (int i = 0; i < I; i++) {								// First-stage part of objective function
		obj_fnc_1 -= (c_part[i] * x[i]);
	}

	for (int s = 0; s < O; s++) {								// Second-stage part of objective function
		IloExpr	temp_expr(env);

		for (int i = 0; i < I; i++) {
			temp_expr += (r_book[i] * y[s][i]);
		}
		obj_fnc_2 += (Omega_probs[s] * temp_expr);				// Expected profit from booking seats under scenario `s`
		temp_expr.end();
	}

	for (int s = 0; s < O; s++) {								// Third-stage part of objective function
		IloExpr	temp_expr(env);

		for (int i = 0; i < I; i++) {
			temp_expr += ((r_show[i] * w[s][i]) - (c_deny[i] * u[s][i]));
		}
		obj_fnc_3 += (Omega_probs[s] * temp_expr);				// Expected profit/penalty from show-up/no-show seats under scenario `s`
		temp_expr.end();
	}

	IloObjective objective_function = IloMaximize(		// Defining an objective function to be maximized
		env,
		(obj_fnc_1 + obj_fnc_2 + obj_fnc_3),
		"Objective_Function");
	model.add(objective_function);						// Adding objective function expression to the model
	obj_fnc_1.end();									// Ending objective fnc expression 1 to prevent memory leaks.
	obj_fnc_2.end();									// Ending objective fnc expression 2 to prevent memory leaks.
	obj_fnc_3.end();									// Ending objective fnc expression 3 to prevent memory leaks.


	/* CREATE SOLUTION OBJECT, WHICH WILL BE RESPONSIBLE TO CALL SOLVING FUNCTION AND RESULTS */
	IloCplex mdl_solution(model);

	/* CPLEX PARAMETERS TUNING */
	if (!show_cplex_display) {		// Turn-off CPLEX logging screen
		mdl_solution.setOut(env.getNullStream());
	}
	if (show_cplex_display) {
		mdl_solution.setParam(IloCplex::Param::MIP::Display, display_intensity);
	}
	if (time_tolerance >= 0) {
		mdl_solution.setParam(IloCplex::Param::TimeLimit, time_tolerance);
	}
	if (gap_tolerance >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, gap_tolerance);
	}
	if (cuts_factor >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::CutsFactor, cuts_factor);
	}
	if (gomory_cuts_tolerance != 0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Cuts::Gomory, gomory_cuts_tolerance);
	}
	if (tree_size_limit >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::TreeMemory, tree_size_limit);
	}
	if (mip_emphasis > 0) {
		mdl_solution.setParam(IloCplex::Param::Emphasis::MIP, mip_emphasis);
	}
	if (thread_limit > 0) {
		mdl_solution.setParam(IloCplex::Param::Threads, thread_limit);
	}
	//mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1e-10);	// e-optimal solution (absolute value)


	/* CALLING SOLVING FUNCTION */
	mdl_solution.solve();


	/* GETTING SOLUTION STATUS */
	solution_status = mdl_solution.getStatus();
	if (solution_status != 3) {
		objective_fnc_value = (double)mdl_solution.getObjValue();
	}
	if (show_cplex_display) {
		env.out() << "\n###########################################################\n";
		env.out() << "CPLEX Solution status = " << mdl_solution.getStatus() << endl;
		if (solution_status != 3) {
			env.out() << "Three-Stage SP objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}

	/* RETRIEVING SOLUTION */
	vector_flt	x_solution(I, -1.0);
	matrix_flt	y_solution(O, vector_flt(I, -1.0));
	matrix_flt	w_solution(O, vector_flt(I, -1.0));
	matrix_flt	u_solution(O, vector_flt(I, -1.0));

	if (solution_status != 3) {

		for (int i = 0; i < I; i++) {
			x_solution[i] = (double)mdl_solution.getValue(x[i]);

			for (int s = 0; s < O; s++) {
				y_solution[s][i] = (double)mdl_solution.getValue(y[s][i]);
				w_solution[s][i] = (double)mdl_solution.getValue(w[s][i]);
				u_solution[s][i] = (double)mdl_solution.getValue(u[s][i]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_vector(x_solution, false, true, "solution_x_i_CPLEX_RP_3S.txt", "Outputs", "\t");
			print_matrix(y_solution, false, true, "solution_y_si_CPLEX_RP_3S.txt", "Outputs", "\t");
			print_matrix(w_solution, false, true, "solution_w_si_CPLEX_RP_3S.txt", "Outputs", "\t");
			print_matrix(u_solution, false, true, "solution_u_si_CPLEX_RP_3S.txt", "Outputs", "\t");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	x.end();
	y.end();
	w.end();
	u.end();
	mdl_solution.end();
	env.end();
	//model.end();
	x_solution.clear();
	x_solution.shrink_to_fit();
	y_solution.clear();
	y_solution.shrink_to_fit();
	w_solution.clear();
	w_solution.shrink_to_fit();
	u_solution.clear();
	u_solution.shrink_to_fit();

	return objective_fnc_value;
}


/**
* Function solves the Northam Airlines three-stage stochastic problem with fixed partitioning (Exercise 1, pp. 83, Birge & Louveaux, 2011) using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param I	The number of different seat types.
* @param Omega	The demand and show-up scenarios for each seat type.
* @param Omega_probs	The probability of each scenario `s`.
* @param x	The for partitioning seats.
* @param c_part	The vector with costs for partitioning.
* @param c_deny	The vector with penalty for denying reservation.
* @param r_book	The vector with revenue per reservation booked of type `i`.
* @param r_show	The vector with revenue per show-up reservation of type `i`.
* @param show_cplex_display	Whether or not to show CPLEX solution display. The default is `true`.
* @param save_results Whether or not to save the results to files. The default is `true`.
* @param display_intensity	Sets different levels of output display. The default is `2`. (0-No display , 1-Display integer feasible solutions , 2-Default , 3-Number cuts , 4-More information , 5-All information)
* @param gap_tolerance	Gap tolerance for stoping criteria. The default is `-1.0` to keep the CPLEX default.
* @param time_tolerance	The maximum limit time, in seconds, for running CPLEX. The defautl is `-1.0` to keep CPLEX default value.
* @param cuts_factor	Limits the number of cuts that can be added. The default is `-1.0`. (-1.0-CPLEX dynamically adjusts the limit, [0.0,1.0]-For values between the range [0.0, 1.0] CPLEX generates no cuts, >1.0-CPLEX limits the number of rows in the model with cuts added)
* @param gomory_cuts_tolerance	Limits the number of cuts added by CPLEX. The default is `0`. (-1-Do not limit, 0-Default, 1-Moderate, 2-Aggressively)
* @param tree_size_limit	Absolute upper limit on the size (in megabytes, uncompressed) of the branch-and-cut tree. The default is `-1.0` to keep CPLEX default value.
* @param mip_emphasis	Controls trade-offs between speed, feasibility, optimality, and moving bounds in MIP. The default is `0`. (0-Balance optimality and feasibility, 1-Feasibility over optimality, 2-Optimality over feasibility, 3-Moving best bound, 4-Finding hidden feasible solutions)
* @param thread_limit	Sets the default maximal number of parallel threads that will be invoked by any CPLEX parallel optimizer. The default is `0`. (0-Automatic, 1-Single thread, int N>1-Uses up to N threads limited by the available processors)
* @return objective_fnc_value	The objective function value for the RP problem.
*/
double Airline_CPLEX_RP_3S_fixed(
	int					I,
	vector<matrix_flt>	Omega,
	vector_flt			Omega_probs,
	vector_flt			x,
	vector_flt			c_part,
	vector_flt			c_deny,
	vector_flt			r_book,
	vector_flt			r_show,
	bool show_cplex_display,//=true,
	bool save_results,//=false
	int display_intensity,//=2,
	double gap_tolerance,//=-1.0,
	int time_tolerance,//=-1.0,
	double cuts_factor,//=-1.0,
	int gomory_cuts_tolerance,//=0,
	double tree_size_limit,//=-1.0,
	int mip_emphasis,//=0,
	int thread_limit//=0
) {

	/* CREATE ENVIRONMENT AND MODELLING OBJECTS */
	IloEnv		env;							// CPLEX environment
	IloModel	model(env, "3S");				// CPLEX model to store problem and naming it "3S"
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	IloNum		lb_u = 0.0;						// u variables lower bound
	IloNum		ub_u = IloInfinity;				// u variables upper bound
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// The CPLEX status of the solution found for the problem


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	int							O = (int)Omega.size();		// Number of scenarios
	IloArray<IloNumVarArray>	y(env, O);					// Decision variables for seats of type `i` reserved under scenario `s`
	IloArray<IloNumVarArray>	w(env, O);					// Decision variables for seats of type `i` show under scenario `s`
	IloArray<IloNumVarArray>	u(env, O);					// Decision variables for seats of type `i` denied under scenario `s`

	for (int s = 0; s < O; s++) {
		y[s] = IloNumVarArray(env, I);							// Assigning a array with the number of types
		w[s] = IloNumVarArray(env, I);							// Assigning a array with the number of types
		u[s] = IloNumVarArray(env, I);							// Assigning a array with the number of types

		for (int i = 0; i < I; i++) {
			y[s][i] = IloNumVar(env, lb_y, ub_y, ILOINT);		// Defining y variables
			w[s][i] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);		// Defining w variables
			u[s][i] = IloNumVar(env, lb_u, ub_u, ILOFLOAT);		// Defining u variables
		}
	}


	/* SECOND-STAGE CONSTRAINTS */

	for (int s = 0; s < O; s++) {
		int s_prime = ((int)s - 1);

		for (int i = 0; i < I; i++) {
			model.add(y[s][i] <= Omega[s][i][0]);							// Adding 'y_si <= xi_si' constraint to model

			if ((s > 0) && (Omega[s][i][0] == Omega[s_prime][i][0])) {
				model.add(y[s][i] == y[s_prime][i]);						// Adding 'y_si - y_s'i == 0' -> nonanticipativity
			}

			model.add(y[s][i] >= rhs_nonneg_ct);							// Adding 'y_si >= 0.0' constraint to model
		}
	}


	/* THIRD-STAGE CONSTRAINTS */

	for (int s = 0; s < O; s++) {

		for (int i = 0; i < I; i++) {

			model.add(w[s][i] == ((Omega[s][i][1] * y[s][i]) - u[s][i]));	// Adding 'w_si == (o_si * y_si) - u_si' constraint to model
			model.add(w[s][i] <= x[i]);										// Adding 'w_si <= x_i' constraint to model

			model.add(w[s][i] >= rhs_nonneg_ct);							// Adding 'w_si >= 0.0' constraint to model
			model.add(u[s][i] >= rhs_nonneg_ct);							// Adding 'u_si >= 0.0' constraint to model
		}
	}

	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr obj_fnc_1(env);										// Expression for first-stage part of the objective function.
	IloExpr	obj_fnc_2(env);										// Expression for second-stage part of the objective function.
	IloExpr	obj_fnc_3(env);										// Expression for third-stage part of the objective function.

	for (int i = 0; i < I; i++) {								// First-stage part of objective function
		obj_fnc_1 -= (c_part[i] * x[i]);
	}

	for (int s = 0; s < O; s++) {								// Second-stage part of objective function
		IloExpr	temp_expr(env);

		for (int i = 0; i < I; i++) {
			temp_expr += (r_book[i] * y[s][i]);
		}
		obj_fnc_2 += (Omega_probs[s] * temp_expr);				// Expected profit from booking seats under scenario `s`
		temp_expr.end();
	}

	for (int s = 0; s < O; s++) {								// Third-stage part of objective function
		IloExpr	temp_expr(env);

		for (int i = 0; i < I; i++) {
			temp_expr += ((r_show[i] * w[s][i]) - (c_deny[i] * u[s][i]));
		}
		obj_fnc_3 += (Omega_probs[s] * temp_expr);				// Expected profit/penalty from show-up/no-show seats under scenario `s`
		temp_expr.end();
	}

	IloObjective objective_function = IloMaximize(		// Defining an objective function to be maximized
		env,
		(obj_fnc_1 + obj_fnc_2 + obj_fnc_3),
		"Objective_Function");
	model.add(objective_function);						// Adding objective function expression to the model
	obj_fnc_1.end();									// Ending objective fnc expression 1 to prevent memory leaks.
	obj_fnc_2.end();									// Ending objective fnc expression 2 to prevent memory leaks.
	obj_fnc_3.end();									// Ending objective fnc expression 3 to prevent memory leaks.


	/* CREATE SOLUTION OBJECT, WHICH WILL BE RESPONSIBLE TO CALL SOLVING FUNCTION AND RESULTS */
	IloCplex mdl_solution(model);

	/* CPLEX PARAMETERS TUNING */
	if (!show_cplex_display) {		// Turn-off CPLEX logging screen
		mdl_solution.setOut(env.getNullStream());
	}
	if (show_cplex_display) {
		mdl_solution.setParam(IloCplex::Param::MIP::Display, display_intensity);
	}
	if (time_tolerance >= 0) {
		mdl_solution.setParam(IloCplex::Param::TimeLimit, time_tolerance);
	}
	if (gap_tolerance >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, gap_tolerance);
	}
	if (cuts_factor >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::CutsFactor, cuts_factor);
	}
	if (gomory_cuts_tolerance != 0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Cuts::Gomory, gomory_cuts_tolerance);
	}
	if (tree_size_limit >= 0.0) {
		mdl_solution.setParam(IloCplex::Param::MIP::Limits::TreeMemory, tree_size_limit);
	}
	if (mip_emphasis > 0) {
		mdl_solution.setParam(IloCplex::Param::Emphasis::MIP, mip_emphasis);
	}
	if (thread_limit > 0) {
		mdl_solution.setParam(IloCplex::Param::Threads, thread_limit);
	}
	//mdl_solution.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1e-10);	// e-optimal solution (absolute value)


	/* CALLING SOLVING FUNCTION */
	mdl_solution.solve();


	/* GETTING SOLUTION STATUS */
	solution_status = mdl_solution.getStatus();
	if (solution_status != 3) {
		objective_fnc_value = (double)mdl_solution.getObjValue();
	}
	if (show_cplex_display) {
		env.out() << "\n###########################################################\n";
		env.out() << "CPLEX Solution status = " << mdl_solution.getStatus() << endl;
		if (solution_status != 3) {
			env.out() << "Three-Stage SP objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}

	/* RETRIEVING SOLUTION */
	matrix_flt	y_solution(O, vector_flt(I, -1.0));
	matrix_flt	w_solution(O, vector_flt(I, -1.0));
	matrix_flt	u_solution(O, vector_flt(I, -1.0));

	if (solution_status != 3) {

		for (int i = 0; i < I; i++) {

			for (int s = 0; s < O; s++) {
				y_solution[s][i] = (double)mdl_solution.getValue(y[s][i]);
				w_solution[s][i] = (double)mdl_solution.getValue(w[s][i]);
				u_solution[s][i] = (double)mdl_solution.getValue(u[s][i]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_vector(x, false, true, "solution_x_i_CPLEX_RP_3S_fixed.txt", "Outputs", "\t");
			print_matrix(y_solution, false, true, "solution_y_si_CPLEX_RP_3S_fixed.txt", "Outputs", "\t");
			print_matrix(w_solution, false, true, "solution_w_si_CPLEX_RP_3S_fixed.txt", "Outputs", "\t");
			print_matrix(u_solution, false, true, "solution_u_si_CPLEX_RP_3S_fixed.txt", "Outputs", "\t");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	y.end();
	w.end();
	u.end();
	mdl_solution.end();
	env.end();
	//model.end();
	y_solution.clear();
	y_solution.shrink_to_fit();
	w_solution.clear();
	w_solution.shrink_to_fit();
	u_solution.clear();
	u_solution.shrink_to_fit();

	return objective_fnc_value;
}


/**
* Function rounds the fractional values in a vector.
* 
* Values are rounded up if the first decimal is >= 5, and rounded down otherwise.
* 
* @param vec	The vector with values to round.
* @return temp	The vector with rounded values.
*/
vector_flt round_values(vector_flt vec) {
	vector_flt temp;

	for (int i = 0; i < vec.size(); i++) {
		int int_part = (int)vec[i];
		double flt_part = (vec[i] - (double)int_part);

		if (flt_part == 0.0) {
			temp.push_back(vec[i]);
		}
		else if ((flt_part * 10.0) >= 5.0) {
			temp.push_back(((double)int_part + 1.0));
		}
		else {
			temp.push_back((double)int_part);
		}
	}

	return temp;
}


/**
* Function cumulates the values of a vector from time t = 1 to time t = `tau`.
* 
* @param arr	The array of IloCplex variables with values to cumulate.
* @param tau	The time until cumulate to.
* @param env	The problem environment.
* @return cumulative	The cumulative value over `t = 1, ..., tau`.
*/
IloExpr cumulate(IloNumVarArray arr, int tau, IloEnv& env) {
	IloExpr cumulative(env, 0.0);

	for (int t = 0; t <= tau; t++) {
		cumulative += arr[t];
	}

	return cumulative;
}



/**
* Function cumulates the values of a vector from time t = 1 to time t = `tau`.
*
* @param vec	The vector with values to cumulate.
* @param tau	The time until cumulate to.
* @return cumulative	The cumulative value over `t = 1, ..., tau`.
*/
IloNum cumulate(vector_flt vec, int tau) {
	IloNum cumulative = 0.0;

	for (int t = 0; t <= tau; t++) {
		cumulative += (IloNum)vec[t];
	}

	return cumulative;
}


/**
* Function calculates the total number of possible scenarios.
*
* @param R	The integer for the number of possible realizations `\xi` at a given time `t`.
* @param H	The integer total time horizon.
* @return O	The number of all possible scenarios.
*/
int num_scenarios(int R, int H) {

	int O = (int)pow((double)R, (double)H);	// Calculating the total number of possible scenarios

	return O;
}


/**
* Function generates the set of all possible scenarios.
* 
* @param R	The integer for the number of possible realizations `\xi` at a given time `t`.
* @param H	The integer total time horizon.
* @return mtrx	The matrix with all possible scenarios
*/
matrix_int get_scenarios_set(int R, int H) {
	int O = (int)pow((double)R, (double)H);								// Calculating the total number of possible scenarios
	matrix_int mtrx(O, vector_int(H, -1));								// Creating an empty matrix to store scenarios

	int s_p = 0;														// Initializing first possible realization index

	for (int t = 0; t < H; t++) {										// Looping through time periods

		int lim = (int)pow((double)R, ((double)H - (double)t - 1.0));	// Calculating, for a given time period, when the realization should change
		
		for (int s = 0; s < O; s++) {									// Looping through scenarios
			
			if ((s > 0) && (s % lim == 0)) {							// Increasing or reseting the realization value
				s_p += 1;

				if (s_p == R) {
					s_p = 0;
				}
			}

			mtrx[s][t] = s_p;											// Assigning realization `\xi` for scenario `s` at period `t`
		}
		s_p = 0;
	}

	return mtrx;
}


/**
* Function calculates the probability of a given scenario `s` in the set of all possible scenarios.
* 
* @param H	The integer total time horizon.
* @param Omega	The matrix with the set of all possible scenarios.
* @param xi_probs	The matrix with the probability of a given realization `\xi` at time `t`.
* @return probs	The vector with the probability of a given scenario `s` in the set of scenarios `\Omega`.
*/
vector_flt get_scenarios_probs(int H, matrix_int Omega, matrix_flt xi_probs) {
	vector_flt probs((int)Omega.size(), 1.0);

	for (int s = 0; s < (int)Omega.size(); s++) {		// Looping through scenarios

		for (int t = 0; t < H; t++) {					// Looping through time periods

			probs[s] *= xi_probs[t][Omega[s][t]];		// Calculating probability of scenario `s` from the independent probabilities of `\xi`
		}
	}

	return probs;
}


/**
* Function creates the matrix with the given weather at time `t` under scenario `s`.
* 
* @param H	The integer total time horizon.
* @param xi	The matrix with the possible weathers per period.
* @param Omega	The matrix with the set/configuration of all possible scenarios.
*/
matrix_flt create_scenarios(int H, matrix_flt xi, matrix_int Omega) {
	matrix_flt	xi_s;

	for (int s = 0; s < (int)Omega.size(); s++) {
		vector_flt	temp_vec;

		for (int t = 0; t < H; t++) {

			temp_vec.push_back(xi[t][Omega[s][t]]);
		}
		xi_s.push_back(temp_vec);
		temp_vec.clear();
	}

	return xi_s;
}



#endif