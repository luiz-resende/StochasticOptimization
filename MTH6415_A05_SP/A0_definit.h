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
* FUNCTION SOLVES THE RP ARILINE PROBLEM
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
);
/**
* FUNCTION SOLVES THE RP ARILINE PROBLEM
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
	bool show_cplex_display = true,
	bool save_results = false,
	int display_intensity = 2,
	double gap_tolerance = -1.0,
	int time_tolerance = -1.0,
	double cuts_factor = -1.0,
	int gomory_cuts_tolerance = 0,
	double tree_size_limit = -1.0,
	int mip_emphasis = 0,
	int thread_limit = 0
);
/**
* FUNCTION SOLVES THE EXPECTED ARILINE PROBLEM
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
	bool show_cplex_display = true,
	bool save_results = false,
	int display_intensity = 2,
	double gap_tolerance = -1.0,
	int time_tolerance = -1.0,
	double cuts_factor = -1.0,
	int gomory_cuts_tolerance = 0,
	double tree_size_limit = -1.0,
	int mip_emphasis = 0,
	int thread_limit = 0
);
/**
* FUNCTION SOLVES THE DETERMINISTIC ARILINE PROBLEM
*/
tuple<vector_flt, double> Airline_CPLEX_deterministic(
	int			S,
	int			I,
	double		Max_seats,
	vector_flt	A,
	vector_flt	c_x,
	vector_flt	r_y,
	vector_flt	xi_mean,
	bool show_cplex_display = true,
	bool save_results = false,
	int display_intensity = 2,
	double gap_tolerance = -1.0,
	int time_tolerance = -1.0,
	double cuts_factor = -1.0,
	int gomory_cuts_tolerance = 0,
	double tree_size_limit = -1.0,
	int mip_emphasis = 0,
	int thread_limit = 0
);
/**
* FUNCTION SOLVES THE WS ARILINE PROBLEM
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
	bool show_cplex_display = true,
	bool save_results = false,
	int display_intensity = 2,
	double gap_tolerance = -1.0,
	int time_tolerance = -1.0,
	double cuts_factor = -1.0,
	int gomory_cuts_tolerance = 0,
	double tree_size_limit = -1.0,
	int mip_emphasis = 0,
	int thread_limit = 0
);
vector_flt round_values(vector_flt vec);
IloExpr cumulate(IloNumVarArray arr, int tau, IloEnv& env);
IloNum cumulate(vector_flt vec, int tau);
int num_scenarios(int R, int H);
matrix_int get_scenarios_set(int R, int H);
vector_flt get_scenarios_probs(int H, matrix_int Omega, matrix_flt xi_probs);
matrix_flt create_scenarios(int H, matrix_flt xi, matrix_int Omega);
void solve_problem_1(double	ticket_economy_price);


/**
* FUNCTION SOLVES THE THREE-STAGE ARILINE PROBLEM
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
	bool show_cplex_display = true,
	bool save_results = false,
	int display_intensity = 2,
	double gap_tolerance = -1.0,
	int time_tolerance = -1.0,
	double cuts_factor = -1.0,
	int gomory_cuts_tolerance = 0,
	double tree_size_limit = -1.0,
	int mip_emphasis = 0,
	int thread_limit = 0
);
/**
* FUNCTION SOLVES THE THREE-STAGE ARILINE PROBLEM WITH FIXED PARTITION
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
	bool show_cplex_display = true,
	bool save_results = false,
	int display_intensity = 2,
	double gap_tolerance = -1.0,
	int time_tolerance = -1.0,
	double cuts_factor = -1.0,
	int gomory_cuts_tolerance = 0,
	double tree_size_limit = -1.0,
	int mip_emphasis = 0,
	int thread_limit = 0
);
void solve_problem_2();


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