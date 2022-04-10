#ifndef A1_PRINT_DATA_CPP
#define A1_PRINT_DATA_CPP

#include "A0_definit.h"


/**
* Prints vector of double values.
*
* @param vec The std::vector constainer with information to print.
* @param transposed Whether or not the vector is to be transposed.
* @param save_to_file Whether to save to a text file or output to screen.
* @param save_file_name String with name for saved file.
* @param folder_name Folder where to save files if save_to_file==true.
* @param sep The value separator string.
*/
void print_vector(vector_flt vec, bool transposed, bool save_to_file, string save_file_name, string folder_name, string sep) {
	string vector_print;

	if (transposed) {
		vector_print += "";

		for (unsigned int i = 0; i < vec.size(); i++) {
			vector_print += (to_string(vec[i]) + sep);
		}
		vector_print = (vector_print.substr(0, vector_print.size() - 1) + "\n");
	}
	else {

		for (unsigned int i = 0; i < vec.size(); i++) {
			vector_print += (to_string(vec[i]) + "\n");
		}
	}

	if (save_to_file) {
		string path_name = (get_date_time() + save_file_name);
		if (folder_name != "") {
			check_directory(folder_name, true);
			path_name = folder_name + "/" + path_name;
		}
		ofstream out_file(path_name);
		out_file << vector_print;
		out_file.close();
	}
	else {
		printf("%s", vector_print.c_str());
	}
}


/**
* Prints vector of int values.
*
* @param vec The std::vector constainer with information to print.
* @param transposed Whether or not the vector is to be transposed.
* @param save_to_file Whether to save to a text file or output to screen.
* @param save_file_name String with name for saved file.
* @param folder_name Folder where to save files if save_to_file==true.
* @param sep The value separator string.
*/
void print_vector(vector_int vec, bool transposed, bool save_to_file, string save_file_name, string folder_name, string sep) {
	string vector_print;

	if (transposed) {
		vector_print += "";

		for (unsigned int i = 0; i < vec.size(); i++) {
			vector_print += (to_string(vec[i]) + sep);
		}
		vector_print = (vector_print.substr(0, vector_print.size() - 1) + "\n");
	}
	else {

		for (unsigned int i = 0; i < vec.size(); i++) {
			vector_print += (to_string(vec[i]) + "\n");
		}
	}

	if (save_to_file) {
		string path_name = (get_date_time() + save_file_name);
		if (folder_name != "") {
			check_directory(folder_name, true);
			path_name = folder_name + "/" + path_name;
		}
		ofstream out_file(path_name);
		out_file << vector_print;
		out_file.close();
	}
	else {
		printf("%s", vector_print.c_str());
	}
}


/**
* Prints matrix of double values.
*
* @param matrix The std::vector<std::vector> constainer with information to print.
* @param transposed Whether or not the matrix is to be transposed.
* @param save_to_file Whether to save to a text file or output to screen.
* @param save_file_name String with name for saved file.
* @param folder_name Folder where to save files if save_to_file==true.
* @param sep The value separator string.
*/
void print_matrix(matrix_flt matrix, bool transposed, bool save_to_file, string save_file_name, string folder_name, string sep) {
	string matrix_print;

	if (!transposed) {
		int n_rows = (int)matrix.size();
		int n_cols = (int)matrix[0].size();

		for (int i = 0; i < n_rows; i++) {
			string temp_string;

			for (int j = 0; j < n_cols; j++) {
				if (j < n_cols - 1) {
					temp_string += (to_string(matrix[i][j]) + sep);
				}
				else {
					temp_string += (to_string(matrix[i][j]) + "\n");
				}
			}
			matrix_print += temp_string;
		}
	}
	else {
		int n_rows = (int)matrix.size();
		int n_cols = (int)matrix[0].size();

		vector_flt temp_vec(n_rows, (double)0.0);
		matrix_flt transpose_matrix;
		for (int col = 0; col < n_cols; col++) {
			transpose_matrix.push_back(temp_vec);
		}

		for (int i = 0; i < n_rows; ++i) {
			for (int j = 0; j < n_cols; ++j) {
				transpose_matrix[j][i] = matrix[i][j];
			}
		}
			
		for (int i = 0; i < n_cols; i++) {
			string temp_string;

			for (int j = 0; j < n_rows; j++) {
				if (j < n_rows - 1) {
					temp_string += (to_string(transpose_matrix[i][j]) + sep);
				}
				else {
					temp_string += (to_string(transpose_matrix[i][j]) + "\n");
				}
			}
			matrix_print += temp_string;
		}
	}

	if (save_to_file) {
		string path_name = (get_date_time() + save_file_name);
		if (folder_name != "") {
			check_directory(folder_name, true);
			path_name = folder_name + "/" + path_name;
		}
		ofstream out_file(path_name);
		out_file << matrix_print;
		out_file.close();
	}
	else {
		printf("%s", matrix_print.c_str());
	}
}


/**
* Prints matrix of int values.
*
* @param matrix The std::vector<std::vector> constainer with information to print.
* @param transposed Whether or not the matrix is to be transposed.
* @param save_to_file Whether to save to a text file or output to screen.
* @param save_file_name String with name for saved file.
* @param folder_name Folder where to save files if save_to_file==true.
* @param sep The value separator string.
*/
void print_matrix(matrix_int matrix, bool transposed, bool save_to_file, string save_file_name, string folder_name, string sep) {
	string matrix_print;

	if (!transposed) {
		int n_rows = (int)matrix.size();
		int n_cols = (int)matrix[0].size();

		for (int i = 0; i < n_rows; i++) {
			string temp_string;

			for (int j = 0; j < n_cols; j++) {
				if (j < n_cols - 1) {
					temp_string += (to_string(matrix[i][j]) + sep);
				}
				else {
					temp_string += (to_string(matrix[i][j]) + "\n");
				}
			}
			matrix_print += temp_string;
		}
	}
	else {
		int n_rows = (int)matrix.size();
		int n_cols = (int)matrix[0].size();

		vector_int temp_vec(n_rows, (int)0);
		matrix_int transpose_matrix;
		for (int col = 0; col < n_cols; col++) {
			transpose_matrix.push_back(temp_vec);
		}

		for (int i = 0; i < n_rows; ++i) {
			for (int j = 0; j < n_cols; ++j) {
				transpose_matrix[j][i] = matrix[i][j];
			}
		}

		for (int i = 0; i < n_cols; i++) {
			string temp_string;

			for (int j = 0; j < n_rows; j++) {
				if (j < n_rows - 1) {
					temp_string += (to_string(transpose_matrix[i][j]) + sep);
				}
				else {
					temp_string += (to_string(transpose_matrix[i][j]) + "\n");
				}
			}
			matrix_print += temp_string;
		}
	}

	if (save_to_file) {
		string path_name = (get_date_time() + save_file_name);
		if (folder_name != "") {
			check_directory(folder_name, true);
			path_name = folder_name + "/" + path_name;
		}
		ofstream out_file(path_name);
		out_file << matrix_print;
		out_file.close();
	}
	else {
		printf("%s", matrix_print.c_str());
	}
}


/**
* Checks and creates directory/folder.
*
* @param dir_name The std::string with the name of the folder to be verified/created.
* @param create_dir Whether or not to create the folder if it does not exist already.
*/
void check_directory(string dir_name, bool create_dir=false) {
	if (filesystem::is_directory(dir_name)) {
		string str = "Directory/folder '" + dir_name + "' exist!\n";
	}
	else {
		string str = "Directory/folder '" + dir_name + "' does not exist...\n";
		if (create_dir) {
			filesystem::create_directory(dir_name);
			string str = "Directory '" + dir_name + "' created!\n";
			printf("%s", str.c_str());
		}
	}
}


/**
* Date and time string.
*
* Function gets current date and time and returns a string in the format
* 'YYYY-MM-DD_hh-mm-ss_' for naming purposes.
* 
* @return date String with date and time in the above format.
*/
string get_date_time() {
	time_t rawtime;
	struct tm timeinfo;
	time(&rawtime);
	localtime_s(&timeinfo, &rawtime);

	vector_int date_int = {
		(timeinfo.tm_year + 1900), (timeinfo.tm_mon + 1), timeinfo.tm_mday, timeinfo.tm_hour, timeinfo.tm_min, timeinfo.tm_sec
	};

	vector<string> date_str;
	for (unsigned int i = 0; i < date_int.size(); i++) {

		if (date_int[i] < 10) {
			date_str.push_back("0" + to_string(date_int[i]));
		}
		else {
			date_str.push_back(to_string(date_int[i]));
		}
	}

	string date = (
		date_str[0] + "-"
		+ date_str[1] + "-"
		+ date_str[2] + "_"
		+ date_str[3] + "h"
		+ date_str[4] + "m"
		+ date_str[5] + "s_"
		);

	return date;
}


/**
* Computes elapsed CPU time.
* 
* Function is called immediately after some computation to calculate the time it required.
* 
* @param start_time The start time of the computation (typename clock_t).
* @return The time in seconds the operation took.
*/
double elapsed_cpu_time(clock_t start_time) {
	clock_t end_time = clock();

	return (((double)end_time - (double)start_time) / (double)CLOCKS_PER_SEC);
}


#endif