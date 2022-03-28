#ifndef A0_DEFINIT_H
#define A0_DEFINIT_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <ilcplex/ilocplex.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <tuple>
#include <vector>


using namespace std;

/* TYPEDEF'N MOST USED CONTAINERS TO CLEAN CODE */
typedef	vector<int>					vector_int;
typedef vector<vector<int>>			matrix_int;
//typedef	vector<float>				vector_flt;
//typedef vector<vector<float>>		matrix_flt;
typedef	vector<double>				vector_flt;
typedef vector<vector<double>>		matrix_flt;
typedef	vector<IloNum>				vector_num;
typedef vector<vector<IloNum>>		matrix_num;



/**
* Function for solving the Clear Lake Dam recourse problem with CPLEX
*/
double Dam_CPLEX_RP(
	int			T,						// The total number of periods
	int			S,						// The number of scenarios in each period
	double		R_max,					// The maximum possible level the dam can release in a given period
	double		L_max,					// The maximum level below flood stage the dam level could be lowered in a given period
	double		L_0,					// The initial/current level below flood stage the dam's water level is in
	double		c_x,					// The cost per millimetre of dam level lowered
	double		q_y,					// The cost per millimetre of water flood
	double		q_w,					// The cost per millimetre of water imported
	matrix_flt	xi_volum,				// The matrix with possible precipitation volumes for all the months in the time horizon
	matrix_flt	xi_probs,				// The matrix with with precipitation probabilities for all months in the time horizon
	int			obj_type = 1,			// If minimizing the expected cost or the probability of violating the limits
	bool show_cplex_display = true,		// Whether or not show the display
	bool save_results = true,			// Whether or not to save the results to files
	int display_intensity = 4,			// Information display intensity
	double gap_tolerance = -1.0,		// Percentage tolerance gap
	int time_tolerance = -1,			// Time limit in seconds
	double cuts_factor = -1.0,			// Cuts factor
	int gomory_cuts_tolerance = 0,		// Gomory cuts tolerance
	double tree_size_limit = -1.0,		// Memory size limit for Branch and Bound tree
	int mip_emphasis = 0,				// MIP emphasis
	int thread_limit = 0				// Limit count on the threads used
);
/**
* Function for solving the Clear Lake Dam deterministic problem with CPLEX
*/
tuple<vector_flt, double> Dam_CPLEX_deterministic(
	int			T,						// The total number of periods
	double		R_max,					// The maximum possible level the dam can release in a given period
	double		L_max,					// The maximum level below flood stage the dam level could be lowered in a given period
	double		L_0,					// The initial/current level below flood stage the dam's water level is in
	double		c_x,					// The cost per millimetre of dam level lowered
	double		q_y,					// The cost per millimetre of water flood
	double		q_w,					// The cost per millimetre of water imported
	vector_flt	xi_mean,				// The vector with the mean evaporation/precipitation volumes for all months
	bool show_cplex_display = true,		// Whether or not show the display
	bool save_results = true,			// Whether or not to save the results to files
	int display_intensity = 4,			// Information display intensity
	double gap_tolerance = -1.0,		// Percentage tolerance gap
	int time_tolerance = -1,			// Time limit in seconds
	double cuts_factor = -1.0,			// Cuts factor
	int gomory_cuts_tolerance = 0,		// Gomory cuts tolerance
	double tree_size_limit = -1.0,		// Memory size limit for Branch and Bound tree
	int mip_emphasis = 0,				// MIP emphasis
	int thread_limit = 0				// Limit count on the threads used
);
/**
* Function for solving the Clear Lake Dam expected value problem with CPLEX
*/
tuple<double, double> Dam_CPLEX_expected(
	int			T,						// The total number of periods
	int			S,						// The number of scenarios in each period
	double		R_max,					// The maximum possible level the dam can release in a given period
	double		L_max,					// The maximum level below flood stage the dam level could be lowered in a given period
	double		L_0,					// The initial/current level below flood stage the dam's water level is in
	double		c_x,					// The cost per millimetre of dam level lowered
	double		q_y,					// The cost per millimetre of water flood
	double		q_w,					// The cost per millimetre of water imported
	matrix_flt	xi_volum,				// The matrix with possible precipitation volumes for all the months in the time horizon
	matrix_flt	xi_probs,				// The matrix with with precipitation probabilities for all months in the time horizon
	int			obj_type = 1,			// If minimizing the expected cost or the probability of violating the limits
	bool show_cplex_display = true,		// Whether or not show the display
	bool save_results = true,			// Whether or not to save the results to files
	int display_intensity = 4,			// Information display intensity
	double gap_tolerance = -1.0,		// Percentage tolerance gap
	int time_tolerance = -1,			// Time limit in seconds
	double cuts_factor = -1.0,			// Cuts factor
	int gomory_cuts_tolerance = 0,		// Gomory cuts tolerance
	double tree_size_limit = -1.0,		// Memory size limit for Branch and Bound tree
	int mip_emphasis = 0,				// MIP emphasis
	int thread_limit = 0				// Limit count on the threads used
);
double Dam_CPLEX_WS(
	int			T,						// The total number of periods
	int			S,						// The number of scenarios in each period
	double		R_max,					// The maximum possible level the dam can release in a given period
	double		L_max,					// The maximum level below flood stage the dam level could be lowered in a given period
	double		L_0,					// The initial/current level below flood stage the dam's water level is in
	double		c_x,					// The cost per millimetre of dam level lowered
	double		q_y,					// The cost per millimetre of water flood
	double		q_w,					// The cost per millimetre of water imported
	matrix_flt	xi_volum,				// The matrix with possible precipitation volumes for all the months in the time horizon
	matrix_flt	xi_probs,				// The matrix with with precipitation probabilities for all months in the time horizon
	int			obj_type = 1,			// If minimizing the expected cost or the probability of violating the limits
	bool show_cplex_display = true,		// Whether or not show the display
	bool save_results = true,			// Whether or not to save the results to files
	int display_intensity = 4,			// Information display intensity
	double gap_tolerance = -1.0,		// Percentage tolerance gap
	int time_tolerance = -1,			// Time limit in seconds
	double cuts_factor = -1.0,			// Cuts factor
	int gomory_cuts_tolerance = 0,		// Gomory cuts tolerance
	double tree_size_limit = -1.0,		// Memory size limit for Branch and Bound tree
	int mip_emphasis = 0,				// MIP emphasis
	int thread_limit = 0				// Limit count on the threads used
);
IloExpr cumulate(IloNumVarArray arr, int tau, IloEnv& env);
IloNum cumulate(vector_flt vec, int tau);
int num_scenarios(int R, int H);
matrix_int get_scenarios_set(int R, int H);
vector_flt get_scenarios_probs(int H, matrix_int Omega, matrix_flt xi_probs);
matrix_flt create_scenarios(int H, matrix_flt xi, matrix_int Omega);



double CPLEX_SP_tomato(
	bool show_cplex_display = true,		// Whether or not show the display
	bool save_results = true,			// Whether or not to save the results to files
	int display_intensity = 4,			// Information display intensity
	double gap_tolerance = -1.0,		// Percentage tolerance gap
	int time_tolerance = -1,			// Time limit in seconds
	double cuts_factor = -1.0,			// Cuts factor
	int gomory_cuts_tolerance = 0,		// Gomory cuts tolerance
	double tree_size_limit = -1.0,		// Memory size limit for Branch and Bound tree
	int mip_emphasis = 0,				// MIP emphasis
	int thread_limit = 0				// Limit count on the threads used
);
double CPLEX_SP_tomato_vss(
	bool show_cplex_display = true,		// Whether or not show the display
	bool save_results = true,			// Whether or not to save the results to files
	int display_intensity = 4,			// Information display intensity
	double gap_tolerance = -1.0,		// Percentage tolerance gap
	int time_tolerance = -1,			// Time limit in seconds
	double cuts_factor = -1.0,			// Cuts factor
	int gomory_cuts_tolerance = 0,		// Gomory cuts tolerance
	double tree_size_limit = -1.0,		// Memory size limit for Branch and Bound tree
	int mip_emphasis = 0,				// MIP emphasis
	int thread_limit = 0				// Limit count on the threads used
);
matrix_flt CPLEX_SP_tomato_determ(
	matrix_flt demands = { {150.0, 35.0, 12.5}, {150.0, 35.0, 12.5}, {150.0, 35.0, 12.5} },
	bool show_cplex_display = true,		// Whether or not show the display
	bool save_results = true,			// Whether or not to save the results to files
	int display_intensity = 4,			// Information display intensity
	double gap_tolerance = -1.0,		// Percentage tolerance gap
	int time_tolerance = -1,			// Time limit in seconds
	double cuts_factor = -1.0,			// Cuts factor
	int gomory_cuts_tolerance = 0,		// Gomory cuts tolerance
	double tree_size_limit = -1.0,		// Memory size limit for Branch and Bound tree
	int mip_emphasis = 0,				// MIP emphasis
	int thread_limit = 0				// Limit count on the threads used
);


/**
* Functions for printing and saving data.
*/
void print_vector(vector_flt vec, bool transposed, bool save_to_file, string save_file_name, string folder_name, string sep);
void print_vector(vector_int vec, bool transposed, bool save_to_file, string save_file_name, string folder_name, string sep);
void print_matrix(matrix_flt matrix, bool transposed, bool save_to_file, string save_file_name, string folder_name, string sep);
void print_matrix(matrix_int matrix, bool transposed, bool save_to_file, string save_file_name, string folder_name, string sep);
void check_directory(string dir_name, bool create_dir);
string get_date_time();
double elapsed_cpu_time(clock_t start_time);



#endif