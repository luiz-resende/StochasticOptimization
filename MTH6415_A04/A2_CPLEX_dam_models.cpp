#ifndef A2_CPLEX_DAM_MODELS_CPP
#define A2_CPLEX_DAM_MODELS_CPP

#include "A0_definit.h"


/**
* Function solves the Clear Lake Dam two-stage stochastic problem (Exercise 3, pp. 50-51, Birge & Louveaux, 2011) using CPLEX library.
* 
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
* 
* @param T	The integer for the time horizon.
* @param S	The integer for the number of realizations in a period.
* @param R_max	The maximum possible level the dam can release in a given period.
* @param L_max	The maximum level below flood stage the dam level could be lowered in a given period.
* @param L_0	The initial/current level below flood stage the dam's water level is in.
* @param c_x	The cost per millimetre of dam level lowered.
* @param q_y	The cost per millimetre of water flood.
* @param q_w	The cost per millimetre of water imported.
* @param xi_volum	The matrix with possible precipitation volumes for all months.
* @param xi_probs	The matrix with precipitation probabilities for all months.
* @param obj_type	The integer describing if minimizing the expected cost (1) or the probability of violating the limits (2). The defult is `1`.
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
double Dam_CPLEX_RP(
	int			T,
	int			S,
	double		R_max,
	double		L_max,
	double		L_0,
	double		c_x,
	double		q_y,
	double		q_w,
	matrix_flt	xi_volum,
	matrix_flt	xi_probs,
	int			obj_type,
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
	IloNum		lb_L = 0.0;						// x variables lower bound
	IloNum		ub_L = IloInfinity;				// x variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	IloNum		rhs_lower_ct = R_max;			// Right-hand-side of lowering capacity constraints
	IloNum		rhs_level_ct = L_max;			// Right-hand-side of maximum water level decision variable
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	IloNum		big_M = 10000000.0;				// Big-M constraint
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// The CPLEX status of the solution found for the problem


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	int							O = num_scenarios(S, T);	// Total number of scenarios
	IloNumVarArray				x(env, T);					// Variables for lowering decision at period `t` (first-stage)
	IloArray<IloNumVarArray>	L(env, O);					// Decision variables for the dam water level at period `t` and scenario `s`
	IloArray<IloNumVarArray>	y(env, O);					// Decision variables for flooded water at period `t` and scenario `s`
	IloArray<IloNumVarArray>	w(env, O);					// Decision variables for imported water at period `t` and scenario `s`
	IloBoolVarArray				u(env, O);					// Decision variables for probability flood at period `t` and scenario `s`
	IloBoolVarArray				v(env, O);					// Decision variables for probability flood imported water at period `t` and scenario `s`
	vector_flt					P_s;						// Probability of scenario `s`
	matrix_int					Omega;						// The set of all possible scenarios
	matrix_flt					xi_omega;					// The list with all weather scenario configuration vectors
	
	Omega = get_scenarios_set(S, T);
	P_s = get_scenarios_probs(T, Omega, xi_probs);
	xi_omega = create_scenarios(T, xi_volum, Omega);

	for (int t = 0; t < T; t++) {
		x[t] = IloNumVar(env, lb_x, ub_x, ILOFLOAT);		// Defining x variables
	}

	for (int s = 0; s < O; s++) {
		L[s] = IloNumVarArray(env, T);						// Assigning a array with the number of periods
		y[s] = IloNumVarArray(env, T);						// Assigning a array with the number of periods
		w[s] = IloNumVarArray(env, T);						// Assigning a array with the number of periods

		u[s] = IloBoolVar(env, 0, 1);						// Defining u variables
		v[s] = IloBoolVar(env, 0, 1);						// Defining v variables
		
		for (int t = 0; t < T; t++) {
			L[s][t] = IloNumVar(env, lb_L, ub_L, ILOFLOAT);		// Defining L variables
			y[s][t] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);		// Defining y variables
			w[s][t] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);		// Defining w variables
		}
	}


	/* FIRST-STAGE CONSTRAINTS */
	for (int t = 0; t < T; t++) {
		model.add(x[t] >= rhs_nonneg_ct);	// Adding 'x_t >= 0.0' constraint to model
		model.add(x[t] <= rhs_lower_ct);	// Adding 'x_t <= R_max' constraint to model
	}


	/* LINKING LEVEL OF THE DAM TO y AND w VARIABLES */
	IloExpr	flood_prob(env);
	IloExpr	import_prob(env);

	for (int s = 0; s < O; s++) {

		for (int t = 0; t < T; t++) {
			int tau = (t - 1);

			if (t == 0) {																	// Calculating the level below flood stage at time `t`
				model.add(L[s][t] == (L_0 - x[t] + xi_omega[s][t] - y[s][t] + w[s][t]));	// Period `t = 1`
			}
			else {
				model.add(L[s][t] == (L[s][tau] - x[t] + xi_omega[s][t] - y[s][t] + w[s][t]));	// Periods `t = 2, ..., |H|`
			}

			model.add(L[s][t] >= rhs_nonneg_ct);		// Adding 'L_st >= 0.0' constraint to model
			model.add(L[s][t] <= rhs_level_ct);			// Adding 'L_st <= L_max' constraint to model

			model.add(y[s][t] >= rhs_nonneg_ct);		// Adding 'y_st >= 0.0' constraint to model
			model.add(w[s][t] >= rhs_nonneg_ct);		// Adding 'w_st >= 0.0' constraint to model

			model.add(y[s][t] <= (big_M * u[s]));		// Activating the binary variable if flood limit is violated
			model.add(w[s][t] <= (big_M * v[s]));		// Activating the binary variable if import limit is violated
		}
		flood_prob += (P_s[s] * u[s]);					// Flood probability
		import_prob += (P_s[s] * v[s]);					// Flood probability
	}

	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);										// Expression for the objective function flood part.
	IloExpr	obj_fnc_2(env);										// Expression for the objective function import part.

	if (obj_type == 1) {

		for (int t = 0; t < T; t++) {
			obj_fnc_1 += (c_x * x[t]);									// Water volume released
		}

		for (int s = 0; s < O; s++) {
			IloExpr	temp_expr(env);

			for (int t = 0; t < T; t++) {
				//temp_expr += ((q_y * y[s][t]) + (q_w * w[s][t]));
				temp_expr += ((u[s] * q_y * y[s][t]) + (v[s] * q_w * w[s][t]));
			}
			obj_fnc_2 += (P_s[s] * temp_expr);							// Expected flood and import water level costs over periods under scenario `s`
			temp_expr.end();
		}
		// obj_fnc_2 += (0.5 * (flood_prob + import_prob));
	}
	else {
		obj_fnc_1 += 0.0;
		obj_fnc_2 += (0.5 * (flood_prob + import_prob));
	}

	IloObjective objective_function = IloMinimize(		// Defining an objective function to be minimized
		env,
		(obj_fnc_1 + obj_fnc_2),
		"Objective_Function");
	model.add(objective_function);						// Adding objective function expression to the model
	obj_fnc_1.end();									// Ending objective fnc expression 1 to prevent memory leaks.
	obj_fnc_2.end();									// Ending objective fnc expression 2 to prevent memory leaks.
	flood_prob.end();
	import_prob.end();


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
	vector_flt	x_solution(T, 0.0);
	matrix_flt	L_solution(O, vector_flt(T, 0.0));
	matrix_flt	y_solution(O, vector_flt(T, 0.0));
	matrix_flt	w_solution(O, vector_flt(T, 0.0));
	vector_int	u_solution(O, 0);
	vector_int	v_solution(O, 0);

	if (solution_status != 3) {
		for (int t = 0; t < T; t++) {
			x_solution[t] = (double)mdl_solution.getValue(x[t]);

			for (int s = 0; s < O; s++) {
				L_solution[s][t] = (double)mdl_solution.getValue(L[s][t]);
				y_solution[s][t] = (double)mdl_solution.getValue(y[s][t]);
				w_solution[s][t] = (double)mdl_solution.getValue(w[s][t]);
				u_solution[s] = (int)mdl_solution.getValue(u[s]);
				v_solution[s] = (int)mdl_solution.getValue(v[s]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_vector(x_solution, false, true, "solution_x_t_CPLEX_RP.txt", "Outputs", " ");
			print_matrix(L_solution, false, true, "solution_L_st_CPLEX_RP.txt", "Outputs", " ");
			print_matrix(y_solution, false, true, "solution_y_st_CPLEX_RP.txt", "Outputs", " ");
			print_matrix(w_solution, false, true, "solution_w_st_CPLEX_RP.txt", "Outputs", " ");
			print_vector(u_solution, false, true, "solution_u_s_CPLEX_RP.txt", "Outputs", " ");
			print_vector(v_solution, false, true, "solution_v_s_CPLEX_RP.txt", "Outputs", " ");
			print_vector(P_s, false, true, "solution_P_s_CPLEX_RP.txt", "Outputs", " ");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	x.end();
	L.end();
	y.end();
	w.end();
	u.end();
	v.end();
	mdl_solution.end();
	env.end();
	//model.end();
	x_solution.clear();
	x_solution.shrink_to_fit();
	L_solution.clear();
	L_solution.shrink_to_fit();
	y_solution.clear();
	y_solution.shrink_to_fit();
	w_solution.clear();
	w_solution.shrink_to_fit();
	u_solution.clear();
	u_solution.shrink_to_fit();
	v_solution.clear();
	v_solution.shrink_to_fit();


	return objective_fnc_value;
}


/**
* Function solves the Clear Lake Dam expected value solution problem (Exercise 3, pp. 50-51, Birge & Louveaux, 2011) using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param T	The integer for the time horizon.
* @param S	The integer for the number of realizations in a period.
* @param R_max	The maximum possible level the dam can release in a given period.
* @param L_max	The maximum level below flood stage the dam level could be lowered in a given period.
* @param L_0	The initial/current level below flood stage the dam's water level is in.
* @param c_x	The cost per millimetre of dam level lowered.
* @param q_y	The cost per millimetre of water flood.
* @param q_w	The cost per millimetre of water imported.
* @param xi_volum	The matrix with possible precipitation volumes for all months.
* @param xi_probs	The matrix with precipitation probabilities for all months.
* @param obj_type	The integer describing if minimizing the expected cost (1) or the probability of violating the limits (2). The defult is `1`.
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
tuple<double, double> Dam_CPLEX_expected(
	int			T,
	int			S,
	double		R_max,
	double		L_max,
	double		L_0,
	double		c_x,
	double		q_y,
	double		q_w,
	matrix_flt	xi_volum,
	matrix_flt	xi_probs,
	int			obj_type,
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
	IloModel	model(env, "RP_expected");		// CPLEX model to store problem and naming it "RP_expected"
	IloNum		lb_L = 0.0;						// y variables lower bound
	IloNum		ub_L = IloInfinity;				// y variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	IloNum		rhs_level_ct = L_max;			// Right-hand-side of maximum water level decision variable
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	IloNum		big_M = 10000000.0;				// Big-M constraint
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	double		objective_fnc_deter = -999.0;	// Variable to store the objective function value of the deterministic model
	int			solution_status;				// The status of the solution for the problem


	/* CREATING VECTOR OF SIZE T WITH MEAN PRECIPITATION/EVAPORATION VOLUMES PER MONTH */
	vector_flt xi_mean;
	for (int t = 0; t < T; t++) {
		double sum = 0.0;

		for (int s = 0; s < S; s++) {
			sum += (xi_volum[t][s] * xi_probs[t][s]);
		}
		xi_mean.push_back(sum);
	}


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	int							O = num_scenarios(S, T);	// Total number of scenarios in the set Omega
	vector_flt					x;							// The release decisions for mean realiations 
	IloArray<IloNumVarArray>	L(env, O);					// Variables for water level at period t and scenario s
	IloArray<IloNumVarArray>	y(env, O);					// Variables for flooded water at period t and scenario s
	IloArray<IloNumVarArray>	w(env, O);					// Variables for imported water at period t and scenario s
	IloBoolVarArray				u(env, O);					// Decision variables for probability flood at period `t` and scenario `s`
	IloBoolVarArray				v(env, O);					// Decision variables for probability flood imported water at period `t` and scenario `s`
	matrix_int					Omega;						// The set of all possible scenarios
	vector_flt					P_s;						// Probability of scenario s
	matrix_flt					xi_omega;					// The list with all weather scenario configuration vectors

	Omega = get_scenarios_set(S, T);
	P_s = get_scenarios_probs(T, Omega, xi_probs);
	xi_omega = create_scenarios(T, xi_volum, Omega);

	for (int s = 0; s < O; s++) {
		L[s] = IloNumVarArray(env, T);							// Assigning a vector with the number of scenarios at period `t`
		y[s] = IloNumVarArray(env, T);							// Assigning a vector with the number of scenarios at period `t`
		w[s] = IloNumVarArray(env, T);							// Assigning a vector with the number of scenarios at period `t`
		
		u[s] = IloBoolVar(env, 0, 1);							// Defining u variables
		v[s] = IloBoolVar(env, 0, 1);							// Defining v variables

		for (int t = 0; t < T; t++) {
			L[s][t] = IloNumVar(env, lb_L, ub_L, ILOFLOAT);		// Defining L variables
			y[s][t] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);		// Defining y variables
			w[s][t] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);		// Defining w variables
		}
	}


	/* GETTING THE x VARIABLE VALUES FOR THE PROBLEM USING MEAN PRECIPITATION/EVAPORATION VALUES */
	tie(x, objective_fnc_deter) = Dam_CPLEX_deterministic(T, R_max, L_max, L_0, c_x, q_y, q_w, xi_mean, false);		// Saves the deterministic result

	/* LINKING LEVEL OF THE DAM TO y AND w VARIABLES */
	IloExpr	flood_prob(env);
	IloExpr	import_prob(env);

	for (int s = 0; s < O; s++) {

		for (int t = 0; t < T; t++) {
			int tau = (t - 1);

			if (t == 0) {							// Level linking constraint `t`
				model.add(L[s][t] == (L_0 - x[t] + xi_omega[s][t] - y[s][t] + w[s][t]));		// Period `t = 1`
			}
			else {
				model.add(L[s][t] == (L[s][tau] - x[t] + xi_omega[s][t] - y[s][t] + w[s][t]));	// Periods `t = 2, ..., |H|`
			}

			model.add(L[s][t] >= rhs_nonneg_ct);	// Adding 'L_st >= 0.0' constraint to model
			model.add(L[s][t] <= rhs_level_ct);		// Adding 'L_st <= L_max' constraint to model

			model.add(y[s][t] >= rhs_nonneg_ct);	// Adding 'y_st >= 0.0' constraint to model
			model.add(w[s][t] >= rhs_nonneg_ct);	// Adding 'w_st >= 0.0' constraint to model

			model.add(y[s][t] <= (big_M * u[s]));	// Activating the binary variable if flood limit is violated
			model.add(w[s][t] <= (big_M * v[s]));	// Activating the binary variable if import limit is violated
		}
		flood_prob += (P_s[s] * u[s]);				// Flood probability
		import_prob += (P_s[s] * v[s]);				// Flood probability
	}

	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);							// Expression for the objective function water level lowering part (first-stage).
	IloExpr	obj_fnc_2(env);							// Expression for the objective function second-stage part.

	if (obj_type == 1) {

		for (int t = 0; t < T; t++) {
			obj_fnc_1 += (c_x * x[t]);					// Expected water release levels cost
		}

		for (int s = 0; s < O; s++) {
			IloExpr	temp_expr(env);

			for (int t = 0; t < T; t++) {
				temp_expr += ((u[s] * q_y * y[s][t]) + (v[s] * q_w * w[s][t]));
			}
			obj_fnc_2 += (P_s[s] * temp_expr);					// Expected flood volume over periods
			temp_expr.end();
		}
		// obj_fnc_2 += (0.5 * (flood_prob + import_prob));
	}
	else {
		obj_fnc_1 += 0.0;
		obj_fnc_2 += (0.5 * (flood_prob + import_prob));
	}

	IloObjective objective_function = IloMinimize(env,
		(obj_fnc_1 + obj_fnc_2),
		"Objective_Function");								// Defining an objective function to be minimized
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
	matrix_flt	L_solution(O, vector_flt(T, 0.0));
	matrix_flt	y_solution(O, vector_flt(T, 0.0));
	matrix_flt	w_solution(O, vector_flt(T, 0.0));
	vector_int	u_solution(O, 0);
	vector_int	v_solution(O, 0);

	if (solution_status != 3) {
		for (int s = 0; s < O; s++) {
			u_solution[s] = (int)mdl_solution.getValue(u[s]);
			v_solution[s] = (int)mdl_solution.getValue(v[s]);

			for (int t = 0; t < T; t++) {
				L_solution[s][t] = (double)mdl_solution.getValue(L[s][t]);
				y_solution[s][t] = (double)mdl_solution.getValue(y[s][t]);
				w_solution[s][t] = (double)mdl_solution.getValue(w[s][t]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_vector(x, false, true, "solution_x_t_CPLEX_expected.txt", "Outputs", " ");
			print_matrix(y_solution, false, true, "solution_y_st_CPLEX_expected.txt", "Outputs", " ");
			print_matrix(w_solution, false, true, "solution_w_st_CPLEX_expected.txt", "Outputs", " ");
			print_matrix(L_solution, false, true, "solution_L_st_CPLEX_expected.txt", "Outputs", " ");
			print_vector(P_s, false, true, "solution_P_s_CPLEX_expected.txt", "Outputs", " ");
			print_vector(u_solution, false, true, "solution_u_s_CPLEX_expected.txt", "Outputs", " ");
			print_vector(v_solution, false, true, "solution_v_s_CPLEX_expected.txt", "Outputs", " ");
		}
	}


	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	y.end();
	w.end();
	L.end();
	u.end();
	v.end();
	mdl_solution.end();
	env.end();
	//model.end();
	x.clear();
	x.shrink_to_fit();
	L_solution.clear();
	L_solution.shrink_to_fit();
	y_solution.clear();
	y_solution.shrink_to_fit();
	w_solution.clear();
	w_solution.shrink_to_fit();
	u_solution.clear();
	u_solution.shrink_to_fit();
	v_solution.clear();
	v_solution.shrink_to_fit();

	return make_tuple(objective_fnc_value, objective_fnc_deter);
}


/**
* Function solves the Clear Lake Dam mean precipitation problem (Exercise 3, pp. 50-51, Birge & Louveaux, 2011) using CPLEX library to get solution
* for calculating the deterministic problem.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param T	The integer for the time horizon.
* @param R_max	The maximum possible level the dam can release in a given period.
* @param L_max	The maximum level below flood stage the dam level could be lowered in a given period.
* @param L_0	The initial/current level below flood stage the dam's water level is in.
* @param c_x	The cost per millimetre of dam level lowered.
* @param q_y	The cost per millimetre of water flood.
* @param q_w	The cost per millimetre of water imported.
* @param xi_mean	The vector with the mean evaporation/precipitation volumes for all months.
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
tuple<vector_flt, double> Dam_CPLEX_deterministic(
	int			T,
	double		R_max,
	double		L_max,
	double		L_0,
	double		c_x,
	double		q_y,
	double		q_w,
	vector_flt	xi_mean,
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
	IloModel	model(env, "SP_deterministic");	// CPLEX model to store problem and naming it "SP_deterministic"
	IloNum		lb_x = 0.0;						// x variables lower bound
	IloNum		ub_x = IloInfinity;				// x variables upper bound
	IloNum		lb_L = 0.0;						// L variables lower bound
	IloNum		ub_L = IloInfinity;				// L variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	IloNum		rhs_level_ct = L_max;			// Right-hand-side of maximum water level decision variable
	IloNum		rhs_lower_ct = R_max;			// Right-hand-side of water level lowering constraints
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// The status of the solution for the problem


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	IloNumVarArray			x(env, T);				// First-stage decision variables - Water level release per month based on mean weather
	IloNumVarArray			L(env, T);				// Second-stage varaibles - Level of water above flood stage at period t for the mean weather
	IloNumVarArray			y(env, T);				// Second-stage varaibles - Level of water above flood stage at period t for the mean weather
	IloNumVarArray			w(env, T);				// Second-stage varaibles - Level of water below 250mm below flood stage for the mean weather

	for (int t = 0; t < T; t++) {

		x[t] = IloNumVar(env, lb_x, ub_x, ILOFLOAT);				// Defining x variables
		L[t] = IloNumVar(env, lb_L, ub_L, ILOFLOAT);				// Defining L variables
		y[t] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);				// Defining y variables
		w[t] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);				// Defining w variables
	}


	/* RELEASE CAPACITY CONSTRAINTS */
	for (int t = 0; t < T; t++) {

		model.add(x[t] >= rhs_nonneg_ct);				// Adding 'x_t >= 0.0' constraint to model
		model.add(x[t] <= rhs_lower_ct);				// Adding 'x_t <= R_max' constraint to model
	}


	/* LINKING LEVEL OF THE DAM TO y AND w VARIABLES */
	for (int t = 0; t < T; t++) {
		int tau = (t - 1);

		if (t == 0) {
			model.add(L[t] == (L_0 - x[t] + xi_mean[t] - y[t] + w[t]));		// Adding linking constraints
		}
		else {
			model.add(L[t] == (L[tau] - x[t] + xi_mean[t] - y[t] + w[t]));
		}

		model.add(L[t] >= rhs_nonneg_ct);			// Adding 'L_t >= 0.0' constraint to model
		model.add(L[t] <= rhs_level_ct);			// Adding 'L_t <= L_max' constraint to model
		model.add(y[t] >= rhs_nonneg_ct);			// Adding 'y_t >= 0.0' constraint to model
		model.add(w[t] >= rhs_nonneg_ct);			// Adding 'w_t >= 0.0' constraint to model
	}


	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);									// Expression for the objective function release part.
	IloExpr	obj_fnc_2(env);									// Expression for the objective function flood/import part.

	for (int t = 0; t < T; t++) {

		obj_fnc_1 += (c_x * x[t]);							// Cost for releasing water over periods
		obj_fnc_2 += ((q_y * y[t]) + (q_w * w[t]));			// Cost for flood/import levels over periods
	}
	IloObjective objective_function = IloMinimize(				// Defining an objective function to be minimized
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
	vector_flt	x_solution(T, -999.0);
	vector_flt	L_solution(T, -999.0);
	vector_flt	y_solution(T, -999.0);
	vector_flt	w_solution(T, -999.0);

	if (solution_status != 3) {											// If problem is not infeasible
		for (int t = 0; t < T; t++) {

			x_solution[t] = (double)mdl_solution.getValue(x[t]);
			L_solution[t] = (double)mdl_solution.getValue(L[t]);
			y_solution[t] = (double)mdl_solution.getValue(y[t]);
			w_solution[t] = (double)mdl_solution.getValue(w[t]);
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_vector(xi_mean, false, true, "solution_mean_xi_t_CPLEX_deterministic.txt", "Outputs", " ");
			print_vector(x_solution, false, true, "solution_x_t_CPLEX_deterministic.txt", "Outputs", " ");
			print_vector(L_solution, false, true, "solution_L_t_CPLEX_deterministic.txt", "Outputs", " ");
			print_vector(y_solution, false, true, "solution_y_t_CPLEX_deterministic.txt", "Outputs", " ");
			print_vector(w_solution, false, true, "solution_w_t_CPLEX_deterministic.txt", "Outputs", " ");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	x.end();
	L.end();
	y.end();
	w.end();
	mdl_solution.end();
	env.end();
	L_solution.clear();
	L_solution.shrink_to_fit();
	y_solution.clear();
	y_solution.shrink_to_fit();
	w_solution.clear();
	w_solution.shrink_to_fit();
	//model.end();

	return make_tuple(x_solution, objective_fnc_value);
}


/**
* Function solves the Clear Lake Dam multi-stage stochastic problem (Exercise 3, pp. 50-51, Birge & Louveaux, 2011) using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param T	The integer for the time horizon.
* @param S	The integer for the number of realizations in a period.
* @param R_max	The maximum possible level the dam can release in a given period.
* @param L_max	The maximum level below flood stage the dam level could be lowered in a given period.
* @param L_0	The initial/current level below flood stage the dam's water level is in.
* @param c_x	The cost per millimetre of dam level lowered.
* @param q_y	The cost per millimetre of water flood.
* @param q_w	The cost per millimetre of water imported.
* @param xi_volum	The matrix with possible precipitation volumes for all months.
* @param xi_probs	The matrix with precipitation probabilities for all months.
* @param obj_type	The integer describing if minimizing the expected cost (1) or the probability of violating the limits (2). The defult is `1`.
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
double Dam_CPLEX_WS(
	int			T,
	int			S,
	double		R_max,
	double		L_max,
	double		L_0,
	double		c_x,
	double		q_y,
	double		q_w,
	matrix_flt	xi_volum,
	matrix_flt	xi_probs,
	int			obj_type,
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
	IloModel	model(env, "WS");				// CPLEX model to store problem and naming it "WS"
	IloNum		lb_x = 0.0;						// x variables lower bound
	IloNum		ub_x = IloInfinity;				// x variables upper bound
	IloNum		lb_L = 0.0;						// x variables lower bound
	IloNum		ub_L = IloInfinity;				// x variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	IloNum		rhs_lower_ct = R_max;			// Right-hand-side of lowering capacity constraints
	IloNum		rhs_level_ct = L_max;			// Right-hand-side of maximum water level decision variable
	IloNum		rhs_nonneg_ct = 0.0;			// Right-hand-side of non-negative constraints
	IloNum		big_M = 10000000.0;				// Big-M constraint
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// The CPLEX status of the solution found for the problem


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	int							O = num_scenarios(S, T);	// Total number of scenarios
	IloArray<IloNumVarArray>	x(env, O);					// Variables for lowering decision at period `t` (first-stage)
	IloArray<IloNumVarArray>	L(env, O);					// Decision variables for the dam water level at period `t` and scenario `s`
	IloArray<IloNumVarArray>	y(env, O);					// Decision variables for flooded water at period `t` and scenario `s`
	IloArray<IloNumVarArray>	w(env, O);					// Decision variables for imported water at period `t` and scenario `s`
	IloBoolVarArray				u(env, O);					// Decision variables for probability flood at period `t` and scenario `s`
	IloBoolVarArray				v(env, O);					// Decision variables for probability flood imported water at period `t` and scenario `s`
	vector_num					P_s;						// Probability of scenario `s`
	matrix_int					Omega;						// The set of all possible scenarios
	matrix_flt					xi_omega;					// The list with all weather scenario configuration vectors

	Omega = get_scenarios_set(S, T);
	P_s = get_scenarios_probs(T, Omega, xi_probs);
	xi_omega = create_scenarios(T, xi_volum, Omega);

	for (int t = 0; t < T; t++) {
		
	}

	for (int s = 0; s < O; s++) {
		x[s] = IloNumVarArray(env, T);						// Assigning a array with the number of periods
		L[s] = IloNumVarArray(env, T);						// Assigning a array with the number of periods
		y[s] = IloNumVarArray(env, T);						// Assigning a array with the number of periods
		w[s] = IloNumVarArray(env, T);						// Assigning a array with the number of periods

		u[s] = IloBoolVar(env, 0, 1);							// Defining u variables
		v[s] = IloBoolVar(env, 0, 1);							// Defining v variables

		for (int t = 0; t < T; t++) {
			x[s][t] = IloNumVar(env, lb_x, ub_x, ILOFLOAT);		// Defining x variables
			L[s][t] = IloNumVar(env, lb_L, ub_L, ILOFLOAT);		// Defining L variables
			y[s][t] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);		// Defining y variables
			w[s][t] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);		// Defining w variables
		}
	}


	/* LINKING LEVEL OF THE DAM TO y AND w VARIABLES */
	IloExpr	flood_prob(env);
	IloExpr	import_prob(env);

	for (int s = 0; s < O; s++) {

		for (int t = 0; t < T; t++) {
			int tau = (t - 1);

			model.add(x[s][t] >= rhs_nonneg_ct);	// Adding 'x_t >= 0.0' constraint to model
			model.add(x[s][t] <= rhs_lower_ct);	// Adding 'x_t <= R_max' constraint to model

			if (t == 0) {																	// Calculating the level below flood stage at time `t`
				model.add(L[s][t] == (L_0 - x[s][t] + xi_omega[s][t] - y[s][t] + w[s][t]));	// Period `t = 1`
			}
			else {
				model.add(L[s][t] == (L[s][tau] - x[s][t] + xi_omega[s][t] - y[s][t] + w[s][t]));	// Periods `t = 2, ..., |H|`
			}

			model.add(L[s][t] >= rhs_nonneg_ct);	// Adding 'L_st >= 0.0' constraint to model
			model.add(L[s][t] <= rhs_level_ct);		// Adding 'L_st <= L_max' constraint to model

			model.add(y[s][t] >= rhs_nonneg_ct);	// Adding 'y_st >= 0.0' constraint to model
			model.add(w[s][t] >= rhs_nonneg_ct);	// Adding 'w_st >= 0.0' constraint to model

			model.add(y[s][t] <= (big_M * u[s]));	// Activating the binary variable if flood limit is violated
			model.add(w[s][t] <= (big_M * v[s]));	// Activating the binary variable if import limit is violated
		}
		flood_prob += (P_s[s] * u[s]);				// Flood probability
		import_prob += (P_s[s] * v[s]);				// Flood probability
	}

	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc(env);										// Expression for the objective function flood part.

	if (obj_type == 1) {

		for (int s = 0; s < O; s++) {
			IloExpr	temp_expr(env);

			for (int t = 0; t < T; t++) {
				temp_expr += ((c_x * x[s][t]) + (u[s] * q_y * y[s][t]) + (v[s] * q_w * w[s][t]));
			}
			obj_fnc += (P_s[s] * temp_expr);							// Expected flood and import water level costs over periods under scenario `s`
			temp_expr.end();
		}
		// obj_fnc += (0.5 * (flood_prob + import_prob));
	}
	else {
		obj_fnc += (0.5 * (flood_prob + import_prob));
	}

	IloObjective objective_function = IloMinimize(		// Defining an objective function to be minimized
		env,
		obj_fnc,
		"Objective_Function");
	model.add(objective_function);						// Adding objective function expression to the model
	obj_fnc.end();										// Ending objective fnc expression 2 to prevent memory leaks.


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
	matrix_flt	x_solution(O, vector_flt(T, 0.0));
	matrix_flt	L_solution(O, vector_flt(T, 0.0));
	matrix_flt	y_solution(O, vector_flt(T, 0.0));
	matrix_flt	w_solution(O, vector_flt(T, 0.0));
	vector_int	u_solution(O, 0);
	vector_int	v_solution(O, 0);

	if (solution_status != 3) {
		for (int t = 0; t < T; t++) {

			for (int s = 0; s < O; s++) {
				x_solution[s][t] = (double)mdl_solution.getValue(x[s][t]);
				L_solution[s][t] = (double)mdl_solution.getValue(L[s][t]);
				y_solution[s][t] = (double)mdl_solution.getValue(y[s][t]);
				w_solution[s][t] = (double)mdl_solution.getValue(w[s][t]);

				u_solution[s] = (int)mdl_solution.getValue(u[s]);
				v_solution[s] = (int)mdl_solution.getValue(v[s]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_matrix(x_solution, false, true, "solution_x_st_CPLEX_EVPI.txt", "Outputs", " ");
			print_matrix(L_solution, false, true, "solution_L_st_CPLEX_EVPI.txt", "Outputs", " ");
			print_matrix(y_solution, false, true, "solution_y_st_CPLEX_EVPI.txt", "Outputs", " ");
			print_matrix(w_solution, false, true, "solution_w_st_CPLEX_EVPI.txt", "Outputs", " ");
			print_vector(u_solution, false, true, "solution_u_s_CPLEX_EVPI.txt", "Outputs", " ");
			print_vector(v_solution, false, true, "solution_v_s_CPLEX_EVPI.txt", "Outputs", " ");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	objective_function.end();
	x.end();
	L.end();
	y.end();
	w.end();
	u.end();
	v.end();
	mdl_solution.end();
	env.end();
	//model.end();
	x_solution.clear();
	x_solution.shrink_to_fit();
	L_solution.clear();
	L_solution.shrink_to_fit();
	y_solution.clear();
	y_solution.shrink_to_fit();
	w_solution.clear();
	w_solution.shrink_to_fit();
	u_solution.clear();
	u_solution.shrink_to_fit();
	v_solution.clear();
	v_solution.shrink_to_fit();

	return objective_fnc_value;
}


/**
* Function solves tomato problem using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
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
*/
double CPLEX_SP_tomato(
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
	IloModel	model(env, "tomato");			// CPLEX model to store problem and naming it "tomato"
	IloNum		lb_x = 0.0;						// x variables lower bound
	IloNum		ub_x = IloInfinity;				// x variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// CPLEX solution status

	int			T = 3;							// Time Horizon
	int			S = 2;							// Number of possible realizations per period
	int			P = 3;							// Number of products
	int			R = 4;							// Number of resources
	vector_num	R_max = {						// Maximum amount of resources
		200.0,
		250.0,
		300.0,
		100.0
	};
	vector_num	y_0 = {							// Initial inventory
		0.0,
		0.0,
		0.0
	};
	vector_num	c_production = {				// Cost of production per box of product
		1.0,
		1.5,
		2.5
	};
	vector_num	c_storage = {					// Cost of storage per box of product
		0.50,
		0.25,
		0.20
	};
	vector_num	c_oportunity = {				// Cost of oportunity for unmet demand
		2.0,
		3.0,
		6.0
	};
	matrix_num	A = {							// Amount of resources necessary per box of product
		{0.50 , 1.00 , 0.00 , 0.25},
		{0.80 , 0.50 , 0.50 , 1.00},
		{1.00 , 0.50 , 1.00 , 3.00}
	};
	vector_num	c_extra = {						// Cost per extra unit of resource
		2.0,
		0.5,
		1.0,
		1.0
	};


	/* CREATING MATRIX OF SIZE TxS WITH PRECIPITATION/EVAPORATION VOLUMES PER MONTH */
	vector<matrix_flt> xi_demand = {
		{
			{200.0 , 100.0},
			{40.0 , 30.0},
			{20.0 , 5.0}
		},
		{
			{200.0 , 100.0},
			{40.0 , 30.0},
			{20.0 , 5.0}
		},
		{
			{200.0 , 100.0},
			{40.0 , 30.0},
			{20.0 , 5.0}
		},
	};

	/* CREATING MATRIX OF SIZE TxS WITH PRECIPITATION/EVAPORATION VOLUME PROBABILITIES PER MONTH */
	vector<matrix_flt> xi_probs = {
		{
			{0.50 , 0.50},
			{0.50 , 0.50},
			{0.50 , 0.50}
		},
		{
			{0.50 , 0.50},
			{0.50 , 0.50},
			{0.50 , 0.50}
		},
		{
			{0.50 , 0.50},
			{0.50 , 0.50},
			{0.50 , 0.50}
		},
	};


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	int									O = (int)pow((double)S, (double)T);	// Total number of scenarios
	IloArray<IloNumVarArray>			x(env, T);							// Products production
	IloArray<IloNumVarArray>			u(env, T);							// Extra resource used
	IloArray<IloArray<IloNumVarArray>>	y(env, O);							// Surplus 
	IloArray<IloArray<IloNumVarArray>>	w(env, O);							// Shortage
	vector_num							P_s;								// Probability of scenarios
	matrix_int							Omega;								// The set of all possible scenarios

	Omega = get_scenarios_set(S, T);
	P_s = {
		0.1250 , 0.1250 , 0.1250 , 0.1250 , 0.1250 , 0.1250 , 0.1250 , 0.1250
	};
	vector<matrix_flt> xi;
	for (int s = 0; s < 8; s++) {
		matrix_flt mtrx;

		for (int t = 0; t < 3; t++) {
			vector_flt vec;

			for (int p = 0; p < 3; p++) {
				vec.push_back(xi_demand[t][p][Omega[s][t]]);
			}
			mtrx.push_back(vec);
			vec.clear();
		}
		xi.push_back(mtrx);
		mtrx.clear();
	}

	for (int t = 0; t < T; t++) {
		x[t] = IloNumVarArray(env, P);

		for (int p = 0; p < P; p++) {
			x[t][p] = IloNumVar(env, lb_x, ub_x, ILOFLOAT);							// Defining x variables
		}
	}

	for (int t = 0; t < T; t++) {
		u[t] = IloNumVarArray(env, R);

		for (int r = 0; r < R; r++) {
			u[t][r] = IloNumVar(env, lb_x, ub_x, ILOFLOAT);							// Defining u varaibles
		}
	}

	for (int s = 0; s < O; s++) {
		y[s] = IloArray<IloNumVarArray>(env, T);
		w[s] = IloArray<IloNumVarArray>(env, T);

		for (int t = 0; t < T; t++) {
			y[s][t] = IloNumVarArray(env, P);
			w[s][t] = IloNumVarArray(env, P);

			for (int p = 0; p < P; p++) {
				y[s][t][p] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);					// Defining y variables
				w[s][t][p] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);					// Defining w variables
			}
		}
	}


	/* FIRST-STAGE CONSTRAINTS */
	for (int t = 0; t < T; t++) {

		for (int r = 0; r < R; r++) {
			IloExpr	temp_p(env);

			for (int p = 0; p < P; p++) {
				temp_p += (A[p][r] * x[t][p]);
			}
			model.add((temp_p - u[t][r]) <= R_max[r]);
			temp_p.end();
		}
	}


	/* LINKING LEVEL OF THE DAM TO y AND w VARIABLES */
	for (int s = 0; s < O; s++) {

		for (int t = 0; t < T; t++) {
			int tau = (t - 1);

			for (int p = 0; p < P; p++) {

				if (t == 0) {																		// Defining inventory level
					model.add((x[t][p] - xi[s][t][p] + y_0[p] - y[s][t][p] + w[s][t][p]) == 0.0);
				}
				else {
					model.add((x[t][p] - xi[s][t][p] + y[s][tau][p] - y[s][t][p] + w[s][t][p]) == 0.0);
				}

				model.add(y[s][t][p] >= lb_y);														// Adding 'y_ts >= 0.0' constraint to model
				model.add(w[s][t][p] >= lb_w);														// Adding 'w_ts >= 0.0' constraint to model
			}
		}
	}


	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);									// Expression for the objective function flood part.
	IloExpr	obj_fnc_2(env);									// Expression for the objective function import part.

	for (int t = 0; t < T; t++) {

		for (int p = 0; p < P; p++) {

			obj_fnc_1 += (c_production[p] * x[t][p]);
		}

		for (int r = 0; r < R; r++) {

			obj_fnc_1 += (c_extra[r] * u[t][r]);
		}
	}

	for (int s = 0; s < O; s++) {
		IloExpr	temp_expr(env);

		for (int t = 0; t < T; t++) {

			for (int p = 0; p < P; p++) {

				temp_expr += ((c_storage[p] * y[s][t][p]) + (c_oportunity[p] * w[s][t][p]));	// Expected shortage and surplus cost
			}
		}
		obj_fnc_2 += (P_s[s] * temp_expr);
		temp_expr.end();
	}

	IloObjective objective_function = IloMinimize(env,
		(obj_fnc_1 + obj_fnc_2),
		"Objective_Function");									// Defining an objective function to be minimized
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
		if (solution_status != 3) {
			env.out() << "SP objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}


	/* RETRIEVING SOLUTION */
	matrix_flt	x_solution(T, vector_flt(P, 0.0));

	if (solution_status != 3) {
		for (int t = 0; t < T; t++) {

			for (int p = 0; p < P; p++) {
				x_solution[t][p] = (double)mdl_solution.getValue(x[t][p]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_matrix(x_solution, false, true, "solution_x_pt_CPLEX_tomatos.txt", "Outputs", " ");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	mdl_solution.end();
	objective_function.end();
	x.end();
	y.end();
	w.end();
	env.end();
	//model.end();
	x_solution.clear();
	x_solution.shrink_to_fit();

	return objective_fnc_value;
}


/**
* Function solves tomato problem using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
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
*/
double CPLEX_SP_tomato_vss(
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
	IloModel	model(env, "expect");			// CPLEX model to store problem and naming it "SP"
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;				// CPLEX solution status

	int			T = 3;							// Time Horizon
	int			S = 2;							// Number of possible realizations per period
	int			P = 3;							// Number of products
	int			R = 4;							// Number of resources
	vector_num	R_max = {						// Maximum amount of resources
		200.0,
		250.0,
		300.0,
		100.0
	};
	vector_num	y_0 = {							// Initial inventory
		0.0,
		0.0,
		0.0
	};
	vector_num	c_production = {				// Cost of production per box of product
		1.0,
		1.5,
		2.5
	};
	vector_num	c_storage = {					// Cost of storage per box of product
		0.50,
		0.25,
		0.20
	};
	vector_num	c_oportunity = {				// Cost of oportunity for unmet demand
		2.0,
		3.0,
		6.0
	};
	matrix_num	A = {							// Amount of resources necessary per box of product
		{0.50 , 1.00 , 0.00 , 0.25},
		{0.80 , 0.50 , 0.50 , 1.00},
		{1.00 , 0.50 , 1.00 , 3.00}
	};
	vector_num	c_extra = {						// Cost per extra unit of resource
		2.0,
		0.5,
		1.0,
		1.0
	};


	/* CREATING MATRIX OF SIZE TxS WITH PRECIPITATION/EVAPORATION VOLUMES PER MONTH */
	vector<matrix_flt> xi_demand = {
		{
			{200.0 , 100.0},
			{40.0 , 30.0},
			{20.0 , 5.0}
		},
		{
			{200.0 , 100.0},
			{40.0 , 30.0},
			{20.0 , 5.0}
		},
		{
			{200.0 , 100.0},
			{40.0 , 30.0},
			{20.0 , 5.0}
		},
	};

	/* CREATING MATRIX OF SIZE TxS WITH PRECIPITATION/EVAPORATION VOLUME PROBABILITIES PER MONTH */
	vector<matrix_flt> xi_probs = {
		{
			{0.50 , 0.50},
			{0.50 , 0.50},
			{0.50 , 0.50}
		},
		{
			{0.50 , 0.50},
			{0.50 , 0.50},
			{0.50 , 0.50}
		},
		{
			{0.50 , 0.50},
			{0.50 , 0.50},
			{0.50 , 0.50}
		},
	};

	//print_matrix(xi_volume[0], false, false, "", "", "\t");

	/* CREATING VECTOR WITH MEAN DEMANDS PER PERIOD */
	matrix_flt xi_mean = {
		{150.0, 35.0, 12.5},
		{150.0, 35.0, 12.5},
		{150.0, 35.0, 12.5}
	};


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	// VARIABLES ARE CREATED ALREADY DEFINING THE CONSTRAINTS THAT x_t, y_ts AND w_ts MUST BE GREATER OR EQUAL TO ZERO FOR ALL t,s
	int									O = (int)pow((double)S, (double)T);
	IloArray<IloNumVarArray>			u(env, T);				// 
	IloArray<IloArray<IloNumVarArray>>	y(env, O);				// Surplus
	IloArray<IloArray<IloNumVarArray>>	w(env, O);				// Shortage
	vector_num							P_s;					// Probability of scenario
	matrix_int							Omega;					// The set of all possible scenarios

	matrix_flt x = CPLEX_SP_tomato_determ(xi_mean, false, false);
	//matrix_flt x = {
	//	{200.0, 40.0, 20.0},
	//	{100.0, 30.0, 15.0},
	//	{100.0, 30.0, 5.0}
	//};

	Omega = get_scenarios_set(S, T);
	P_s = {
		0.1250 , 0.1250 , 0.1250 , 0.1250 , 0.1250 , 0.1250 , 0.1250 , 0.1250
	};
	vector<matrix_flt> xi;
	for (int s = 0; s < 8; s++) {
		matrix_flt mtrx;

		for (int t = 0; t < 3; t++) {
			vector_flt vec;

			for (int p = 0; p < 3; p++) {
				vec.push_back(xi_demand[t][p][Omega[s][t]]);
			}
			mtrx.push_back(vec);
			vec.clear();
		}
		xi.push_back(mtrx);
		mtrx.clear();
	}

	for (int t = 0; t < T; t++) {
		u[t] = IloNumVarArray(env, R);

		for (int r = 0; r < R; r++) {
			u[t][r] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
		}
	}

	for (int s = 0; s < O; s++) {
		y[s] = IloArray<IloNumVarArray>(env, T);
		w[s] = IloArray<IloNumVarArray>(env, T);

		for (int t = 0; t < T; t++) {
			y[s][t] = IloNumVarArray(env, P);
			w[s][t] = IloNumVarArray(env, P);

			for (int p = 0; p < P; p++) {
				y[s][t][p] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);	// Defining surplus variables
				w[s][t][p] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);	// Defining shortage variables
			}
		}
	}


	/* FIRST-STAGE CONSTRAINTS */
	for (int t = 0; t < T; t++) {

		for (int r = 0; r < R; r++) {
			IloExpr	temp_p(env);

			for (int p = 0; p < P; p++) {
				temp_p += (A[p][r] * x[t][p]);
			}
			model.add((temp_p - u[t][r]) <= R_max[r]);
			temp_p.end();
		}
	}


	/* LINKING LEVEL OF THE DAM TO y AND w VARIABLES */
	for (int s = 0; s < O; s++) {

		for (int t = 0; t < T; t++) {
			int tau = (t - 1);

			for (int p = 0; p < P; p++) {

				if (t == 0) {										// Defining inventory level
					model.add((x[t][p] - xi[s][t][p] + y_0[p] - y[s][t][p] + w[s][t][p]) == 0.0);
				}
				else {
					model.add((x[t][p] - xi[s][t][p] + y[s][tau][p] - y[s][t][p] + w[s][t][p]) == 0.0);
				}

				model.add(y[s][t][p] >= lb_y);						// Adding 'y_ts >= 0.0' constraint to model
				model.add(w[s][t][p] >= lb_w);						// Adding 'w_ts >= 0.0' constraint to model
			}
		}
	}


	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);									// Expression for the objective function flood part.
	IloExpr	obj_fnc_2(env);									// Expression for the objective function import part.

	for (int t = 0; t < T; t++) {

		for (int p = 0; p < P; p++) {

			obj_fnc_1 += (c_production[p] * x[t][p]);
		}

		for (int r = 0; r < R; r++) {

			obj_fnc_1 += (c_extra[r] * u[t][r]);
		}
	}

	for (int s = 0; s < O; s++) {
		IloExpr	temp_expr(env);

		for (int t = 0; t < T; t++) {

			for (int p = 0; p < P; p++) {

				temp_expr += ((c_storage[p] * y[s][t][p]) + (c_oportunity[p] * w[s][t][p]));		// Expected shortage and surplus
			}
		}
		obj_fnc_2 += (P_s[s] * temp_expr);
		temp_expr.end();
	}

	IloObjective objective_function = IloMinimize(env,
		(obj_fnc_1 + obj_fnc_2),
		"Objective_Function");									// Defining an objective function to be minimized
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
		if (solution_status != 3) {
			env.out() << "SP objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}


	/* RETRIEVING SOLUTION */
	vector<matrix_flt>	y_solution(O, matrix_flt(T, vector_flt(P, 0.0)));
	vector<matrix_flt>	w_solution(O, matrix_flt(T, vector_flt(P, 0.0)));

	if (solution_status != 3) {
		for (int s = 0; s < O; s++) {

			for (int t = 0; t < T; t++) {

				for (int p = 0; p < P; p++) {
					y_solution[s][t][p] = (double)mdl_solution.getValue(y[s][t][p]);
					w_solution[s][t][p] = (double)mdl_solution.getValue(w[s][t][p]);
				}
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			//print_matrix(y_solution, false, true, "solution_x_pt_CPLEX_tomatos.txt", "Outputs", " ");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	mdl_solution.end();
	objective_function.end();
	y.end();
	w.end();
	env.end();
	//model.end();
	y_solution.clear();
	y_solution.shrink_to_fit();
	w_solution.clear();
	w_solution.shrink_to_fit();

	return objective_fnc_value;
}


/**
* Function solves tomato problem using CPLEX library.
*
* CPLEX parameters -> https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-list-parameters
*
* @param volumes	The mean demands.
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
*/
matrix_flt CPLEX_SP_tomato_determ(
	matrix_flt demands,
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
	IloModel	model(env, "det");				// CPLEX model to store problem and naming it "det"
	IloNum		lb_x = 0.0;						// x variables lower bound
	IloNum		ub_x = IloInfinity;				// x variables upper bound
	IloNum		lb_y = 0.0;						// y variables lower bound
	IloNum		ub_y = IloInfinity;				// y variables upper bound
	IloNum		lb_w = 0.0;						// w variables lower bound
	IloNum		ub_w = IloInfinity;				// w variables upper bound
	double		objective_fnc_value = -999.0;	// Variable to store the objective function value of the model
	int			solution_status;

	int			T = 3;							// Time Horizon
	int			P = 3;							// Number of products
	int			R = 4;							// Number of resources
	vector_num	R_max = {						// Maximum amount of resources
		200.0,
		250.0,
		300.0,
		100.0
	};
	vector_num	y_0 = {							// Initial inventory
		0.0,
		0.0,
		0.0
	};
	vector_num	c_production = {				// Cost of production per box of product
		1.0,
		1.5,
		2.5
	};
	vector_num	c_storage = {					// Cost of storage per box of product
		0.50,
		0.25,
		0.20
	};
	vector_num	c_oportunity = {				// Cost of oportunity for unmet demand
		2.0,
		3.0,
		6.0
	};
	matrix_num	A = {							// Amount of resources necessary per box of product
		{0.50 , 1.00 , 0.00 , 0.25},
		{0.80 , 0.50 , 0.50 , 1.00},
		{1.00 , 0.50 , 1.00 , 3.00}
	};
	vector_num	c_extra = {						// Cost per extra unit of resource
		2.0,
		0.5,
		1.0,
		1.0
	};


	/* DEFINING AND CREATING OPTIMIZATION DECISION VARIABLES */
	IloArray<IloNumVarArray>	x(env, T);
	IloArray<IloNumVarArray>	u(env, T);
	IloArray<IloNumVarArray>	y(env, T);
	IloArray<IloNumVarArray>	w(env, T);

	for (int t = 0; t < T; t++) {
		x[t] = IloNumVarArray(env, P);
		y[t] = IloNumVarArray(env, P);
		w[t] = IloNumVarArray(env, P);

		for (int p = 0; p < P; p++) {
			x[t][p] = IloNumVar(env, lb_x, ub_x, ILOFLOAT);					// Defining x production variables
			y[t][p] = IloNumVar(env, lb_y, ub_y, ILOFLOAT);					// Defining y surplus variables
			w[t][p] = IloNumVar(env, lb_w, ub_w, ILOFLOAT);					// Defining w shortage variables
		}

		u[t] = IloNumVarArray(env, R);										

		for (int r = 0; r < R; r++) {
			u[t][r] = IloNumVar(env, lb_x, ub_x, ILOFLOAT);					// Defining u extra resource variables
		}
	}


	/* RESOURCE CONSTRAINTS */
	for (int t = 0; t < T; t++) {

		for (int r = 0; r < R; r++) {
			IloExpr	temp_p(env);

			for (int p = 0; p < P; p++) {
				temp_p += (A[p][r] * x[t][p]);
			}
			model.add((temp_p - u[t][r]) <= R_max[r]);
			temp_p.end();
		}
	}


	/* LINKING CONSTRAINT - INVENTORY, SURPLUS AND SHORTAGE VARIABLES */
	for (int t = 0; t < T; t++) {
		int tau = (t - 1);

		for (int p = 0; p < P; p++) {
			if (t == 0) {															// Defining inventory level
				model.add((x[t][p] - demands[t][p] + y_0[p] - y[t][p] + w[t][p]) == 0.0);
			}
			else {
				model.add((x[t][p] - demands[t][p] + y[tau][p] - y[t][p] + w[t][p]) == 0.0);
			}
			model.add(y[t][p] >= lb_y);												// Adding 'y_ts >= 0.0' constraint to model
			model.add(w[t][p] >= lb_w);												// Adding 'w_ts >= 0.0' constraint to model
		}
	}


	/* DEFINING AND ADDING OBJECTIVE FUNCTION TO THE MODEL */
	IloExpr	obj_fnc_1(env);									// Expression for the objective function surplus part.
	IloExpr	obj_fnc_2(env);									// Expression for the objective function shortage part.

	for (int t = 0; t < T; t++) {

		for (int p = 0; p < P; p++) {

			obj_fnc_1 += (c_production[p] * x[t][p]);
		}

		for (int r = 0; r < R; r++) {

			obj_fnc_1 += (c_extra[r] * u[t][r]);
		}
	}

	for (int t = 0; t < T; t++) {

		for (int p = 0; p < P; p++) {

			obj_fnc_2 += ((c_storage[p] * y[t][p]) + (c_oportunity[p] * w[t][p]));	// Expected flood volume over periods
		}
	}

	IloObjective objective_function = IloMinimize(env,
		(obj_fnc_1 + obj_fnc_2),
		"Objective_Function");									// Defining an objective function to be minimized
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
		if (solution_status != 3) {
			env.out() << "SP objective value = " << (double)mdl_solution.getObjValue() << endl;
		}
		env.out() << "###########################################################\n";
	}


	/* RETRIEVING SOLUTION */
	matrix_flt	x_solution(T, vector_flt(P, 0.0));

	if (solution_status != 3) {
		for (int t = 0; t < T; t++) {

			for (int p = 0; p < P; p++) {
				x_solution[t][p] = (double)mdl_solution.getValue(x[t][p]);
			}
		}

		/* SAVING SOLUTION FOUND TO TEXT FILES */
		if (save_results) {
			print_matrix(x_solution, false, true, "solution_x_tp_CPLEX_tomatos_determ.txt", "Outputs", " ");
		}
	}

	/* DESTROY THE ENVIRONMENT OBJECT AFTER THE PROBLEM HAS ENDED AND INFORMATION FROM IT (VARIABLES, VALUES ETC.) ARE NO LONGER NEEDED */
	mdl_solution.end();
	objective_function.end();
	x.end();
	env.end();
	//model.end();

	return x_solution;
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
