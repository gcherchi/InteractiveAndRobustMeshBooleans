#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

//#define USE_FLOATS

#ifdef USE_FLOATS
typedef float fpnumber;
#define FPN_EPSILON	FLT_EPSILON
#else
typedef double fpnumber;
#define FPN_EPSILON	DBL_EPSILON
#endif

void error(const char *s, int line)
{
	if (line) { cerr << "Line " << line << ": " << s << endl; exit(0); }
	else cerr << s << endl; exit(0);
}

void tokenize(string& line, vector<string>& tokens, char separator)
{
	tokens.clear();
	stringstream check1(line);
	string intermediate;
	while (getline(check1, intermediate, separator)) tokens.push_back(intermediate); 
}

void createHeadingComment(string& s)
{
	s += "/****************************************************************************\n";
	s += "* Indirect predicates for geometric constructions					        *\n";
	s += "*                                                                           *\n";
	s += "* Consiglio Nazionale delle Ricerche                                        *\n";
	s += "* Istituto di Matematica Applicata e Tecnologie Informatiche                *\n";
	s += "* Sezione di Genova                                                         * \n";
	s += "* IMATI-GE / CNR                                                            * \n";
	s += "*                                                                           *\n";
	s += "* Authors: Marco Attene                                                     * \n";
	s += "* Copyright(C) 2019: IMATI-GE / CNR                                         * \n";
	s += "* All rights reserved.                                                      * \n";
	s += "*                                                                           *\n";
	s += "* This program is free software; you can redistribute it and/or modify      * \n";
	s += "* it under the terms of the GNU Lesser General Public License as published  * \n";
	s += "* by the Free Software Foundation; either version 3 of the License, or (at  * \n";
	s += "* your option) any later version.                                           * \n";
	s += "*                                                                           *\n";
	s += "* This program is distributed in the hope that it will be useful, but       * \n";
	s += "* WITHOUT ANY WARRANTY; without even the implied warranty of                * \n";
	s += "* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  * \n";
	s += "* General Public License for more details.                                  * \n";
	s += "*                                                                           *\n";
	s += "* You should have received a copy of the GNU Lesser General Public License  * \n";
	s += "* along with this program.  If not, see http://www.gnu.org/licenses/.       *\n";
	s += "*                                                                           *\n";
	s += "****************************************************************************/ \n";

	s += "\n/* This code was generated automatically. Do not edit unless you exactly   */\n";
	s += "/* know what you are doing!                                                */\n\n";
}

#define MAX_STACK_SIZE	16384 // In bytes, includes all: input parameters, doubles, ints, ...
#define MAX_STATIC_SIZE	128 // This is the maximum in any case, but might be reduced to avoid stack overflows
#define MAX_WORTH_DEGREE 20 // If the polynomial degree is higher than this, avoid FP arithmetic at all

#define UNDERFLOW_GUARDING // If defined, exact functions check for underflows and try to fix

class lambda_variable
{
public:
	string name;		// Name of the implicit point type (e.g. implict2D_SSI)
	string ipoint_name;	// Name of the implicit point variable (e.g. p1)
	vector<int> output_pars; // Indices in all_vars of the results of getLambda() (e.s. l1x, l1y, d)

	lambda_variable(string& n) : name(n) {}

	void print_filtered(ofstream& file);
	void print_interval(ofstream& file);
	void print_bigfloat(ofstream& file);
	void print_exact(ofstream& file);
};

class variable
{
public:
	string name;			// Variable name
	int op1, op2;			// Operands (indexes in all_vars. -1 if input val)
	char op;				// Operation (+, -, *)
	int size;				// Variable size (expansion max length)
	string actual_length;	// Variable length (expansion actual length)

	fpnumber value_bound;		// Forward error analysis: bound on magnitude
	fpnumber error_bound;		// Forward error analysis: bound on error
	int error_degree;		// Forward error analysis: degree
	bool is_a_max;			// TRUE if magnitude is relevant for error bound
	bool is_lambda_out;		// TRUE if this var is the output of a lambda calculation
	bool is_used;			// TRUE if this var is the operand of another var
	bool non_zero;			// TRUE if this var is the operand of another var

	static bool append;		// TRUE to append code to an existing file
	static vector<variable> all_vars;	// List of all the variables
	static vector<variable*> sign_vars;	// List of variables that determine the sign (denominators with odd exponent)
	static vector<lambda_variable> all_lambda_vars;	// List of all the lambda variables
	static string func_name; // Name of the function
	static bool is_indirect; // True if function is an indirect predicate
	static bool is_lambda;	 // True if function defines a lambda
	static bool is_fma;		 // TRUE if at least an FMA instruction is in expression
	static bool boolean_predicate;	 // True if function defines a boolean predicate

	static fpnumber ulp(fpnumber d)
	{
		int exp;
#ifdef USE_FLOATS
		frexpf(d, &exp);
		fpnumber u = ldexpf(0.5, exp - 23);
#else
		frexp(d, &exp);
		fpnumber u = ldexp(0.5, exp - 52);
#endif
		//u += (u / (1 << 11)); // Comment this out to ignore INTEL's x387 extended precision
		return u;
	}

	static void init(bool _append) 
	{
		append = _append;
		is_indirect = is_lambda = is_fma = boolean_predicate = false;
		all_vars.push_back(variable(string("2"))); // "2" is a special variable representing the constant 2
		all_vars.push_back(variable(string("1"))); // "1" is a special variable representing the constant 2
	}

	// Plain declaration
	variable(string& s) : is_a_max(false), is_lambda_out(false), is_used(false), non_zero(false)
	{
		name = s;

		size = 1;
		value_bound = 1;
		error_bound = 0;
		error_degree = 1;

		op1 = op2 = -1; 
		actual_length = name + "_len";

		if (name == "2") value_bound = 2;
	}

	// Plain lambda declaration
	variable(bool is_lambda_v, string& s) : is_a_max(false), is_lambda_out(true), is_used(false), non_zero(false)
	{
		name = s;

		size = 128;
		value_bound = 1;
		error_bound = 0;
		error_degree = 1;

		op1 = op2 = -1;
		actual_length = name + "_len";

		if (name == "2") value_bound = 2;
		is_indirect = true;
	}

	// Plain assignment
	variable(string& s, int o1) : name(s), op1(o1), op2(-1), is_a_max(false), is_lambda_out(false), is_used(false), non_zero(false)
	{
		variable& ov1 = all_vars[o1];
		ov1.is_used = true;

		size = ov1.size;
		value_bound = ov1.value_bound;
		error_bound = ov1.error_bound;
		error_degree = ov1.error_degree;
		if (ov1.error_bound == 0) ov1.is_a_max = true;

		if (size > 2 * MAX_STATIC_SIZE) size = 2 * MAX_STATIC_SIZE;
		actual_length = to_string(size);
	}

	// Binary operation
	variable(string& s, int o1, int o2, char o) : name(s), op1(o1), op2(o2), op(o), is_a_max(false), is_lambda_out(false), is_used(false), non_zero(false)
	{
		variable& ov1 = all_vars[o1];
		variable& ov2 = all_vars[o2];
		ov1.is_used = ov2.is_used = true;

		if (o == '+' || o == '-')
		{
			size = ov1.size + ov2.size;
			if (o == '-' && ov1.error_bound == 0 && ov2.error_bound == 0) // Translation filter
			{
				// This is the same formula used in FPG. I am not sure it is correct.
				// Why is the value bound 1 and not 2? (1 - (-1)) = 2.
				value_bound = 1;
				error_bound = value_bound*(0.5*FPN_EPSILON);
				error_degree = 1;
				is_a_max = true;
			}
			else // Regular sum or subtraction
			{
				value_bound = ov1.value_bound + ov2.value_bound;
				fpnumber u = 0.5*ulp(value_bound);
				value_bound += u;
				error_bound = ov1.error_bound + ov2.error_bound + u;
				error_degree = std::max(ov1.error_degree, ov2.error_degree);

				if (ov1.error_bound == 0) ov1.is_a_max = true;
				if (ov2.error_bound == 0) ov2.is_a_max = true;
			}
		}
		else if (o == '*')
		{
			if (ov1.name == string("2")) {
				size = ov2.size;
				value_bound = ov1.value_bound * ov2.value_bound;
				error_bound = ov2.error_bound;
				error_degree = ov2.error_degree;
			}
			else {
				size = 2 * ov1.size * ov2.size;
				value_bound = ov1.value_bound * ov2.value_bound;
				fpnumber u = 0.5 * ulp(value_bound);
				value_bound += u;
				// The formula herebelow is slightly tighter than FPG's. Original as in FPG is commented below.
				error_bound = ov1.error_bound * ov2.error_bound +
					ov1.error_bound * (ov2.value_bound - ov2.error_bound) +
					ov2.error_bound * (ov1.value_bound - ov1.error_bound) +
					u;
				//error_bound = ov1.error_bound * ov2.error_bound +
				//			  ov1.error_bound * ov2.value_bound +
				//			  ov2.error_bound *ov1.value_bound +
				//			  u;
				error_degree = ov1.error_degree + ov2.error_degree;
			}

			if (ov1.error_bound == 0) ov1.is_a_max = true;
			if (ov2.error_bound == 0) ov2.is_a_max = true;
		}
		else error("Unknown operand", 0);
		if (size > 2 * MAX_STATIC_SIZE) size = 2 * MAX_STATIC_SIZE;
		actual_length = to_string(size);
	}


	bool isInput() const { return (op1 == -1 && op2 == -1); }

	static int getVarByName(string name)
	{
		for (int i = 0; i < (int)all_vars.size(); i++) if (all_vars[i].name == name) return i;
		return -1;
	}


	static void produceAllCode(string& func_name)
	{
		string filtered_funcname = func_name + "_filtered";
		string interval_funcname = func_name + "_interval";
		string exact_funcname = func_name + "_exact";
		string bigfloat_funcname = func_name + "_bigfloat";
		string heading_comment;

		createHeadingComment(heading_comment);

		bool prevHeader = false;
		ifstream checkopen("indirect_predicates.h");
		if (checkopen.is_open())
		{
			prevHeader = true;
			checkopen.close();
		}
		ofstream header("indirect_predicates.h", ios_base::app);

		if (!prevHeader)
		{
			header << heading_comment;
			header << "#include \"implicit_point.h\"\n\n";
		}

		if (is_lambda)
		{
			//filteredPrototype(filtered_funcname, header);
			intervalPrototype(interval_funcname, header);
			exactPrototype(exact_funcname, header);
			bigfloatPrototype(bigfloat_funcname, header);
		}
		else
		{
			multistagePrototype(func_name, header);
		}

		ofstream file;
		if (append)
		{
			bool prevFile = false;

			ifstream checkopen("indirect_predicates.hpp");
			if (checkopen.is_open())
			{
				prevFile = true;
				checkopen.close();
			}

			file.open("indirect_predicates.hpp", ios_base::app);

			if (!prevFile)
			{
				file << heading_comment;
				file << "#include \"implicit_point.h\"\n\n";
				file << "#pragma intrinsic(fabs)\n\n";
#ifdef UNDERFLOW_GUARDING
				file << "// Uncomment the following to activate overflow/underflow checks\n";
				file << "#define CHECK_FOR_XYZERFLOWS\n\n";
#endif
			}
		}
		else
		{
			file.open(func_name + ".hpp");
			file << "#include \"implicit_point.h\"\n\n";
			file << "#pragma intrinsic(fabs)\n\n";
#ifdef UNDERFLOW_GUARDING
			file << "// Uncomment the following to activate overflow/underflow checks\n";
			file << "#define CHECK_FOR_XYZERFLOWS\n\n";
#endif
		}

		if (!is_lambda && !is_indirect) produceFilteredCode(filtered_funcname, file);
		produceIntervalCode(interval_funcname, file);
		produceBigfloatCode(bigfloat_funcname, file);
		produceExactCode(exact_funcname, file);
		
		if (!is_lambda) produceMultiStageCode(func_name, filtered_funcname, interval_funcname, exact_funcname, file);
	}

	bool isParameter() const {
		return isInput() && name != "1" && name != "2" && !is_lambda_out;
	}

	static string createParameterProtoList(string mytype)
	{
		bool first = true;
		stringstream file;
		for (lambda_variable& l : all_lambda_vars) { file << ((first) ? ("const ") : (", const ")) << l.name << "& " << l.ipoint_name; first = false; }
		for (variable& v : all_vars) if (v.isParameter()) { file << ((first) ? (mytype + " ") : (", " + mytype + " ")) << v.name; first = false; }
		return file.str();
	}

	static string createParameterValueList()
	{
		bool first = true;
		stringstream file;
		for (lambda_variable& l : all_lambda_vars) { file << ((first) ? ("") : (", ")) << l.ipoint_name; first = false; }
		for (variable& v : all_vars) if (v.isParameter()) { file << ((first) ? ("") : (", ")) << v.name; first = false; }
		return file.str();
	}

	static void produceMultiStageCode(string& func_name, string& filtered_funcname, string& interval_funcname, string& exact_funcname, ofstream& file)
	{
		bool worth_a_semistatic_filter = ((*(all_vars.end() - 1)).error_degree <= MAX_WORTH_DEGREE);

		string all_pars = createParameterValueList();
		string comm = (worth_a_semistatic_filter) ? ("") : ("//");
		char* return_type = "inline int ";
		file << return_type << func_name << "(" << createParameterProtoList("double") << ")\n{\n";

		file << "   int ret;\n";
		if (!is_lambda && !is_indirect) {
			file << comm << "   ret = " << filtered_funcname << "(" << all_pars << ");\n";
			file << comm << "   if (ret != Filtered_Sign::UNCERTAIN) return ret;\n";
		}
		file << "   ret = " << interval_funcname << "(" << all_pars << ");\n";
		file << "   if (ret != Filtered_Sign::UNCERTAIN) return ret;\n";
		file << "   return " << exact_funcname << "(" << all_pars << ");\n";

		file << "}\n\n";
	}


	static void produceSemiStaticFilter(fpnumber epsilon, int degree, string& threshold_name, ofstream& file)
	{
		file << "   double " << threshold_name << " = max_var;\n";
		int td = 1;
		int l2d = std::floor(std::log2(degree));
		for (int i = 0; i < l2d; i++) { file << "   " << threshold_name << " *= " << threshold_name << ";\n"; td += td; }
		for (int i = 0; i < (degree-td); i++) file << "   " << threshold_name << " *= max_var;\n";
		file << "   " << threshold_name << " *= " << std::setprecision(std::numeric_limits<fpnumber>::digits10 + 1) << epsilon << ";\n";
	}

	static void produceFilteredCode(string& funcname, ofstream& file)
	{
		const string double_type = (is_fma) ? ("ALG16_Double") : ("double");
		char* return_type = (is_lambda || boolean_predicate) ? ("inline bool ") : ("inline int ");

		bool maxes = false;

		if (is_lambda)
		{
			file << return_type << funcname << "(";
			bool first = true;
			for (int i = 1; i < (int)all_vars.size(); i++)
			{
				variable& v = all_vars[i];
				if (v.name.substr(0, 6) == "lambda") file << ", double& " << v.name;
				if (v.isInput() && v.is_used && v.name != "1" && v.name != "2") { file << ((first) ? ("") : (", ")) << "double " << v.name; first = false; }
			}

			file << ", double& max_var)\n{\n";
		}
		else file << return_type << funcname << "(" << createParameterProtoList("double") << ")\n{\n";

		for (int i = 1; i < (int)all_vars.size(); i++) if (all_vars[i].is_a_max) maxes = true;

		// Calculate lambdas

		if (is_indirect)
		{
			bool first;
			file << "   " << double_type;
			first = true; for (variable& v : all_vars) if (v.isInput() && v.is_lambda_out) { file << ((first) ? (" ") : (", ")) << v.name; first = false; }
			file << ", max_var = 0;\n";

			file << "    if (\n";
			first = true; for (lambda_variable& l : all_lambda_vars) { file << ((first) ? ("       !") : ("       || !")); l.print_filtered(file); file << "\n"; first = false; }
			file << "   ) return Filtered_Sign::UNCERTAIN;\n\n";
		}

		// Function body - operations

		int i;
		variable* v;
		for (i = 1; i < (int)all_vars.size(); i++)
		{
			v = &all_vars[i];
			if (v->isInput()) continue;
			string& o1 = all_vars[v->op1].name;

			if (is_lambda && v->name.substr(0, 6) == "lambda") file << "   " << v->name;
			else file << "   const " << double_type << " " << v->name;

			if (v->op2 == -1) {
				file << " = " << o1 << ";\n";
			}
			else {
				string& o2 = all_vars[v->op2].name;

				file << " = " << o1 << " " << v->op << " " << o2 << ";\n";
			}
		}


		// Function body - filter calculation

		if (is_lambda) file << "\n   double _tmp_fabs;\n";
		else
		{
			if (maxes) file << "\n   double _tmp_fabs;\n";
			if (!is_indirect) file << "\n   double max_var = 0.0;\n";
		}

		for (i = 1; i < (int)all_vars.size(); i++) if (all_vars[i].is_a_max)
			file << "   if ((_tmp_fabs = fabs(" << all_vars[i].name << ")) > max_var) max_var = _tmp_fabs;\n";

		fpnumber epsilon = 0;
		int degree = 0;
		string eps_name = "epsilon";
		string outvar_name;
		for (i = 1; i < (int)all_vars.size(); i++) {
			epsilon = all_vars[i].error_bound;
			degree = all_vars[i].error_degree;
			outvar_name = all_vars[i].name;
			if (all_vars[i].name.substr(0, 8) == "lambda_d" && (!is_lambda || all_vars[i].op2 != -1)) {
				eps_name = all_vars[i].name + "_eps"; break;
			}
		}

		if (is_lambda && eps_name == "epsilon") // This should happen only when this is a lambda function and lambda_d = 1
		{
			file << "\n   return true;\n}\n\n";
			return;
		}

		produceSemiStaticFilter(epsilon, degree, eps_name, file);

		// Function body - compare with filter and return

		if (is_lambda) file << "\n   return ( " << "(" << outvar_name << " > " << outvar_name << "_eps || " << outvar_name << " < -" << outvar_name << "_eps) );\n";
		else if (boolean_predicate)
		{
			file << "\n   return ( ";
			bool first = true;
			for (i = 1; i < (int)all_vars.size(); i++) {
				v = &all_vars[i];
				if (v->non_zero)
				{
					if (!first) file << "|| ";
					file << "(" << v->name << " > epsilon || -" << v->name << " < epsilon) ";
					first = false;
				}
			}
			file << ");\n";
		}
		else
		{
			if (!sign_vars.empty())
			{
				file << "   if ((";
				bool first = true; for (variable* v : sign_vars) { file << ((first) ? ("") : (" + ")) << "(" << v->name << " < 0)"; first = false; }
				file << ((sign_vars.size() > 1) ? ") & 1) " : ")) ") << outvar_name << " = -" << outvar_name << ";\n";
			}

			file << "   if (" << outvar_name << " > epsilon) return IP_Sign::POSITIVE;\n";
			file << "   if (-" << outvar_name << " > epsilon) return IP_Sign::NEGATIVE;\n";
			file << "   return Filtered_Sign::UNCERTAIN;\n";
		}

		// Function end
		file << "}\n\n";
	}

	static void produceIntervalCode(string& funcname, ofstream& file)
	{
		char* return_type = (is_lambda || boolean_predicate) ? ("inline bool ") : ("inline int ");

		if (is_lambda)
		{
			file << return_type << funcname << "(";
			bool first = true;
			for (int i = 1; i < (int)all_vars.size(); i++)
			{
				variable& v = all_vars[i];
				if (v.name.substr(0, 6) == "lambda") file << ", interval_number& " << v.name;
				if (v.isInput() && v.is_used && v.name != "1" && v.name != "2") { file << ((first) ? ("") : (", ")) << "interval_number " << v.name; first = false; }
			}

			file << ")\n{\n";
		}
		else file << return_type << funcname << "(" << createParameterProtoList("interval_number") << ")\n{\n";

		// Calculate lambdas

		if (is_indirect)
		{
			file << "   interval_number";
			bool first = true; for (variable& v : all_vars) if (v.isInput() && v.is_lambda_out) { file << ((first) ? (" ") : (", ")) << v.name; first = false; }
			file << ";\n";
			file << "   if (\n";
			first = true; for (lambda_variable& l : all_lambda_vars) { file << ((first) ? ("   !") : ("   || !")); l.print_interval(file); file << "\n"; first = false; }
			file << "   ) return Filtered_Sign::UNCERTAIN;\n\n";
		}


		// Function body - operations

		file << "   setFPUModeToRoundUP();\n";

		variable* v;
		for (int i = 1; i < (int)all_vars.size(); i++)
		{
			v = &all_vars[i];
			if (v->isInput()) continue;
			string expr;
			string& o1 = all_vars[v->op1].name;

			if (v->op2 == -1) {
				expr = o1;
			}
			else {
				string& o2 = all_vars[v->op2].name;
				expr = ((o1 == "2") ? (o2 + " " + v->op + " " + o1) : (o1 + " " + v->op + " " + o2));
			}
			if (is_lambda && v->name.substr(0, 6) == "lambda") file << "   " << v->name << " = " << expr << ";\n";
			else file << "   const interval_number " << v->name << "(" << expr << ");\n";
		}

		file << "   setFPUModeToRoundNEAR();\n\n";

		// Function body - filter

		if (is_lambda)
		{
			file << "   return ";
			for (int i = 1; i < (int)all_vars.size(); i++)
			{
				variable& v = all_vars[i];
				if (v.name.substr(0, 8) == "lambda_d")
				{
					if (v.op2 == -1) file << "true;\n";
					else file << v.name << ".signIsReliable();\n";
					break;
				}
			}

			//file << "\n   return (\n";
			//bool first = true;
			//for (int i = 1; i < (int)all_vars.size(); i++)
			//{
			//	variable& v = all_vars[i];
			//	if (v.name.substr(0, 8) == "lambda_d")
			//	{
			//		if (!first) file << "      && "; else file << "      ";
			//		file << v.name << ".signIsReliable()\n";
			//		first = false;
			//	}
			//}
			//file << "   );\n";
		}
		else if (boolean_predicate)
		{
			file << "   return (\n";
			bool first = true;
			for (int i = 1; i < (int)all_vars.size(); i++)
			{
				variable& v = all_vars[i];
				if (v.non_zero)
				{
					if (!first) file << "      || "; else file << "      ";
					file << v.name << ".signIsReliable()\n";
					first = false;
				}
			}
			file << "   );\n";
		}
		else
		{
			file << "   if (!" << v->name << ".signIsReliable()) return Filtered_Sign::UNCERTAIN;\n";

			if (!sign_vars.empty())
			{
				file << "   if ((";
				bool first = true; for (variable* v : sign_vars) { file << ((first) ? ("") : (" + ")) << "(" << v->name << " < 0)"; first = false; }
				file << ((sign_vars.size() > 1) ? ") & 1) " : ")) ") << "return -" << v->name << ".sign();\n   else";
			}

			file << "   return " << v->name << ".sign();\n";
		}
		file << "}\n\n";
	}


	static void produceBigfloatCode(string& funcname, ofstream& file)
	{
		char* return_type;
		if (is_lambda) return_type = "inline void ";
		else if (boolean_predicate) return_type = "inline bool ";
		else return_type = "inline int ";

		if (is_lambda)
		{
			file << return_type << funcname << "(";
			bool first = true;
			for (int i = 1; i < (int)all_vars.size(); i++)
			{
				variable& v = all_vars[i];
				if (v.name.substr(0, 6) == "lambda") file << ", bigfloat& " << v.name;
				if (v.isInput() && v.is_used && v.name != "1" && v.name != "2") { file << ((first) ? ("") : (", ")) << "bigfloat " << v.name; first = false; }
			}

			file << ")\n{\n";
		}
		else file << return_type << funcname << "(" << createParameterProtoList("bigfloat") << ")\n{\n";

		// Calculate lambdas
		
		if (is_indirect)
		{
			file << "   bigfloat";
			bool first = true; for (variable& v : all_vars) if (v.isInput() && v.is_lambda_out) { file << ((first) ? (" ") : (", ")) << v.name; first = false; }
			file << ";\n";
			for (lambda_variable& l : all_lambda_vars) {
				file << "   ";
				l.print_bigfloat(file); 
				file << ";\n";
			}
		}


		// Function body - operations

		variable* v;
		for (int i = 1; i < (int)all_vars.size(); i++)
		{
			v = &all_vars[i];
			if (v->isInput()) continue;
			string expr;
			string& o1 = all_vars[v->op1].name;

			if (v->op2 == -1) {
				expr = o1;
			}
			else {
				string& o2 = all_vars[v->op2].name;
				expr = ((o1 == "2") ? (o2 + " " + v->op + " " + o1) : (o1 + " " + v->op + " " + o2));
			}

			if (is_lambda && v->name.substr(0, 6) == "lambda") file << "   " << v->name << " = " << expr << ";\n";
			else file << "   const bigfloat " << v->name << "(" << expr << ");\n";
		}

		// Function body - filter

		if (boolean_predicate)
		{
			file << "\n   return true;\n";
		}
		else if (!is_lambda)
		{
			if (!sign_vars.empty())
			{
				file << "   if ((";
				bool first = true; for (variable* v : sign_vars) { file << ((first) ? ("") : (" + ")) << "(" << v->name << " < 0)"; first = false; }
				file << ((sign_vars.size() > 1) ? ") & 1) " : ")) ") << "return -sgn(" << v->name << ");\n   else";
			}

			file << "   return sgn(" << v->name << ");\n";
		}
		file << "}\n";
		file << "\n";
	}


	static void produceExactCode(string& funcname, ofstream& file)
	{
		// Calculate stack size
		uint32_t fixed_stacksize = 152; // Accoount for expansionObject + return_value
		for (variable& v : all_vars) if (v.isParameter()) fixed_stacksize += 8; // Size of a pointer on 64bit systems
		if (is_lambda) for (variable& v : all_vars) if (v.name.substr(0, 6) == "lambda") fixed_stacksize += 16; // two pointers
		else for (variable& v : all_vars) if (v.isInput() && v.is_lambda_out) fixed_stacksize += 16; // two pointers

		if (fixed_stacksize >= MAX_STACK_SIZE) error("Too many parameters - stack overflow unavoidable.\n", 0);

		uint32_t local_max_size = MAX_STATIC_SIZE * 2;
		uint32_t variable_stacksize;

		do
		{
			local_max_size /= 2;
			variable_stacksize = 0;
			for (variable& v : all_vars)
				if ((!is_lambda || v.name.substr(0, 6) != "lambda") && v.name != "2" && v.name != "1")
				{
					if (v.size > local_max_size) variable_stacksize += (local_max_size + 1) * 8;
					else variable_stacksize += (v.size * 8);
				}
				else variable_stacksize += 4;
		} while ((fixed_stacksize + variable_stacksize) >= MAX_STACK_SIZE);
		printf("STACK --- %d\n", fixed_stacksize + variable_stacksize);

		// Return type and function name
		char* return_type = (is_lambda) ? ("inline void ") : ("inline int ");
		if (is_lambda)
		{
			file << return_type << funcname << "(";
			bool first = true;
			for (variable& v : all_vars) if (v.isInput() && v.name != "1" && v.name != "2" && v.is_used && !v.is_lambda_out) { file << ((first) ? ("double ") : (", double ")) << v.name; first = false; }

			for (variable& v : all_vars) { if (v.name.substr(0, 6) == "lambda") file << ((first) ? ("double **") : (", double **")) << v.name << ", int& " << v.name << "_len"; first = false; }

			file << ")\n{\n";
		}
		else file << return_type << funcname << "(" << createParameterProtoList("double") << ")\n{\n";


		if (is_indirect)
		{
			if (boolean_predicate) file << " int return_value = IP_Sign::UNDEFINED;\n";
			else file << " double return_value = NAN;\n";

#ifdef UNDERFLOW_GUARDING
			if (is_indirect && !boolean_predicate)
			{
				file << "#ifdef CHECK_FOR_XYZERFLOWS\n";
				file << "   feclearexcept(FE_ALL_EXCEPT);\n";
				file << "#endif\n";
			}
#endif
			bool first;
			first = true; for (variable& v : all_vars) if (v.isInput() && v.is_lambda_out)
			{
//				if (v.size > local_max_size) 
					file << ((first) ? (" double ") : (", ")) << v.name << "_p[" << local_max_size << "], *" << v.name << " = " << v.name << "_p";
//				else file << ((first) ? (" double ") : (", ")) << v.name << "[" << v.size << "]";
				first = false;
			}
			file << ";\n";
			first = true; for (variable& v : all_vars) if (v.isInput() && v.is_lambda_out)
			{
				file << ((first) ? (" int ") : (", ")) << v.name << "_len = " << local_max_size; first = false;
			}
			file << ";\n";

			for (lambda_variable& l : all_lambda_vars)
			{
				file << " "; l.print_exact(file); file << ";\n";
			}

			file << " if (";
			first = true; for (lambda_variable& l : all_lambda_vars)
			{
				file << ((first) ? ("(") : (" && (")) << all_vars[l.output_pars.back()].name << "[" << all_vars[l.output_pars.back()].name << "_len - 1] != 0)";
				first = false;
			}
			file << ")\n {\n";
		}

		int i;

		for (i = 1; i < (int)all_vars.size(); i++) if (!all_vars[i].isInput()) break;

		// Function body
		//file << "   expansionObject o;\n";

		variable* v;
		int first_op = i;

		for (i = 1; i < (int)all_vars.size(); i++)
		{
			v = &all_vars[i];
			if (v->isInput()) continue;
			variable& op1 = all_vars[v->op1];
			string o1 = ((op1.name.substr(0, 6) == "lambda") ? ("*") : ("")) + op1.name;
			int s1 = op1.size;
			string& al1 = op1.actual_length;
			string lendec = "   int ";
			if (!is_lambda || v->name.substr(0, 6) != "lambda")
			{
				if (v->size > local_max_size) file << "   double " << v->name << "_p[" << local_max_size << "], *" << v->name << " = " << v->name << "_p;\n";
				else file << "   double " << v->name << "[" << v->size << "];\n";
			}
			else lendec = "   ";

			if (v->op2 == -1) {
				if (o1 == string("1"))
				{
					file << "   (*" << v->name << ")[0] = 1;\n";
					file << "   " << v->name << "_len = 1;\n";
					v->actual_length = "1";
				}
				else error("Only plain assignment = 1 is supported!\n", -1);

				continue;
			}

			variable& op2 = all_vars[v->op2];
			string o2 = ((op2.name.substr(0, 6) == "lambda") ? ("*") : ("")) + op2.name;
			int s2 = op2.size;
			string& al2 = op2.actual_length;

			string alb, rlb, lms;
			if (v->name.substr(0, 6) == "lambda") {
				alb = "";
				rlb = "*";
				lms = v->name + "_len";
			}
			else {
				alb = "&";
				rlb = "";
				stringstream lls;
				lls << local_max_size;
				lms = lls.str();
			}

			// Casi noti
			// [k],n *
			if (o1 == string("2"))
			{
				if (v->size > local_max_size) file << lendec << v->name << "_len = expansionObject::Double_With_PreAlloc(" << al2 << ", " << o2 << ", " << alb << v->name << ", " << lms << ");\n";
				else file << "   expansionObject::Double(" << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
				v->actual_length = al2;// v->name + "_len";
			}
			// 1,1 +-*^
			else if (s1 == 1 && s2 == 1)
			{
				if (v->op == '+') file << "   expansionObject::two_Sum(" << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else if (v->op == '-') file << "   expansionObject::two_Diff(" << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else if (v->op == '*' && v->op1 != v->op2) file << "   expansionObject::Two_Prod(" << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else if (v->op == '*' && v->op1 == v->op2) file << "   expansionObject::Square(" << o1 << ", " << rlb << v->name << ");\n";
			}
			// 2,1 *
			else if (s1 == 2 && s2 == 1)
			{
				if (v->op == '*') file << "   expansionObject::Two_One_Prod(" << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else if (v->op == '-') file << "   expansionObject::two_One_Diff(" << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
			}
			// 1,2 *
			else if (s1 == 1 && s2 == 2)
			{
				if (v->op == '*') file << "   expansionObject::Two_One_Prod(" << o2 << ", " << o1 << ", " << rlb << v->name << ");\n";
			}
			// 2,2 +-*^
			else if (s1 == 2 && s2 == 2 && v->op != '*') // Add Two_Square
			{
				if (v->op == '+') file << "   expansionObject::Two_Two_Sum(" << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else if (v->op == '-') file << "   expansionObject::Two_Two_Diff(" << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else if (v->op == '*') file << "   expansionObject::Two_Two_Prod(" << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
			}
			// n,n +-*
			else
			{
				if (v->size > local_max_size)
				{
					//uint32_t lms = (!is_lambda || v->name.substr(0, 6) != "lambda") ? local_max_size : MAX_STATIC_SIZE;
					if (!is_lambda || v->name.substr(0, 6) != "lambda") {
						stringstream lls;
						lls << local_max_size;
						lms = lls.str();
					}
					if (v->op == '+')
					{
						if (s2 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Sum_With_PreAlloc(" << al1 << ", " << o1 << ", " << 1 << ", &" << o2 << ", " << alb << v->name << ", " << lms << ");\n";
						else if (s1 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Sum_With_PreAlloc(" << 1 << ", &" << o1 << ", " << al2 << ", " << o2 << ", " << alb << v->name << ", " << lms << ");\n";
						else file << lendec << v->name << "_len = expansionObject::Gen_Sum_With_PreAlloc(" << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", " << alb << v->name << ", " << lms << ");\n";
					}
					else if (v->op == '-')
					{
						if (s2 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Diff_With_PreAlloc(" << al1 << ", " << o1 << ", " << 1 << ", &" << o2 << ", " << alb << v->name << ", " << lms << ");\n";
						else if (s1 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Diff_With_PreAlloc(" << 1 << ", &" << o1 << ", " << al2 << ", " << o2 << ", " << alb << v->name << ", " << lms << ");\n";
						else file << lendec << v->name << "_len = expansionObject::Gen_Diff_With_PreAlloc(" << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", " << alb << v->name << ", " << lms << ");\n";
					}
					else if (v->op == '*')
					{
						if (s2 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Scale_With_PreAlloc(" << al1 << ", " << o1 << ", " << o2 << ", " << alb << v->name << ", " << lms << ");\n";
						else if (s1 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Scale_With_PreAlloc(" << al2 << ", " << o2 << ", " << o1 << ", " << alb << v->name << ", " << lms << ");\n";
						else file << lendec << v->name << "_len = expansionObject::Gen_Product_With_PreAlloc(" << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", " << alb << v->name << ", " << lms << ");\n";
					}
				}
				else
					if (v->op == '+')
					{
						if (s2 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Sum(" << al1 << ", " << o1 << ", " << 1 << ", &" << o2 << ", " << rlb << v->name << ");\n";
						else if (s1 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Sum(" << 1 << ", &" << o1 << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
						else file << lendec << v->name << "_len = expansionObject::Gen_Sum(" << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
					}
					else if (v->op == '-')
					{
						if (s2 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Diff(" << al1 << ", " << o1 << ", " << 1 << ", &" << o2 << ", " << rlb << v->name << ");\n";
						else if (s1 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Diff(" << 1 << ", &" << o1 << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
						else file << lendec << v->name << "_len = expansionObject::Gen_Diff(" << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
					}
					else if (v->op == '*')
					{
						if (s2 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Scale(" << al1 << ", " << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
						else if (s1 == 1) file << lendec << v->name << "_len = expansionObject::Gen_Scale(" << al2 << ", " << o2 << ", " << o1 << ", " << rlb << v->name << ");\n";
						else file << lendec << v->name << "_len = expansionObject::Gen_Product(" << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
					}
				v->actual_length = v->name + "_len";
			}
		}

		file << "\n";

		if (boolean_predicate)
		{
			file << ((is_indirect)?("   return_value = ( "):(" return ( "));
			bool first = true;
			for (i = 1; i < (int)all_vars.size(); i++) {
				v = &all_vars[i];
				if (v->non_zero)
				{
					if (!first) file << "|| ";
					file << "(" << v->name << "[" << v->actual_length << " - 1] != 0) ";
					first = false;
				}
			}
			file << ");\n";

			if (is_indirect)
			{
				file << " }\n\n";
				for (variable& v : all_vars) if (v.isInput() && v.is_lambda_out)
				{
					if (v.size > local_max_size) file << " if (" << v.name << "_p != " << v.name << ") FreeDoubles(" << v.name << ");\n";
				}
				file << " return return_value;\n";
			}
		}
		else
		{
			if (!is_lambda) file << ((is_indirect) ? ("   ") : ("   double ")) << "return_value = " << v->name << "[" << v->actual_length << " - 1];\n";

			for (i--; i >= first_op; i--)
			{
				v = &all_vars[i];
				if (v->size > local_max_size && (!is_lambda || v->name.substr(0, 6) != "lambda")) file << "   if (" << v->name << "_p != " << v->name << ") FreeDoubles(" << v->name << ");\n";
			}

			if (!is_lambda)
			{
				if (!sign_vars.empty())
				{
					file << "   if (( ";
					bool first = true; for (variable* v : sign_vars) { file << ((first) ? ("") : (" + ")) << "(" << v->name << "[" << v->name << "_len -1] < 0)"; first = false; }
					file << ((sign_vars.size() > 1) ? ") & 1) " : ")) ") << "return_value = -return_value;\n";
				}

				if (is_indirect)
				{
					file << " }\n\n";
#ifdef UNDERFLOW_GUARDING
					file << "#ifdef CHECK_FOR_XYZERFLOWS\n";
					file << "   if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) return_value = ";

					string bigfloatname = funcname.substr(0, funcname.size() - 6) + "_bigfloat";

					file << bigfloatname << "(" << createParameterValueList() << ");\n";
					file << "#endif\n\n";
#endif
					for (variable& v : all_vars) if (v.isInput() && v.is_lambda_out)
					{
						if (v.size > local_max_size) file << " if (" << v.name << "_p != " << v.name << ") FreeDoubles(" << v.name << ");\n";
					}
				}

				file << "\n if (return_value > 0) return IP_Sign::POSITIVE;\n";
				file << " if (return_value < 0) return IP_Sign::NEGATIVE;\n";
				file << " if (return_value == 0) return IP_Sign::ZERO;\n";
				file << " return IP_Sign::UNDEFINED;\n";
			}
		}

		file << "}\n\n";
	}


	static void multistagePrototype(string& func_name, ofstream& file)
	{
		bool first;
		char* return_type = "int ";

		if (is_indirect)
		{
			file << return_type << func_name << "(";
			first = true;
			for (lambda_variable& l : all_lambda_vars) { file << ((first) ? ("const ") : (", const ")) << l.name << "& " << l.ipoint_name; first = false; }
			for (variable& v : all_vars) if (v.isInput() && v.name != "1" && v.name != "2" && !v.is_lambda_out) { file << ((first) ? ("double ") : (", double ")) << v.name; first = false; }
			file << ");\n";
		}
		else
		{
			file << return_type << func_name << "(";

			int i;
			string cur_par, all_pars, mvs = ");\n";

			bool first = true;
			for (i = 1; i < (int)all_vars.size(); i++)
			{
				variable& v = all_vars[i];
				if (v.isInput() && v.name != "1" && v.name != "2")
				{
					if (first) cur_par = v.name;
					else cur_par = ", " + v.name;
					if (first) file << "double " << v.name;
					else file << ", double " << v.name;
					all_pars += cur_par;
					first = false;
				}
			}
			file << ");\n";
		}
	}

	static void filteredPrototype(string& funcname, ofstream& file)
	{
		file << ((is_lambda || boolean_predicate) ? ("bool ") : ("int ")) << funcname << "(";

		bool first = true;
		bool maxes = false;
		for (int i = 1; i < (int)all_vars.size(); i++)
		{
			variable& v = all_vars[i];
			if (v.is_a_max) maxes = true;
			if (is_lambda && v.name.substr(0, 6) == "lambda") file << ", double& " << v.name;
			if (v.isInput() && v.is_used && v.name != "1" && v.name != "2")
			{
				if (first) file << "double " << v.name;
				else file << ", double " << v.name;
				first = false;
			}
		}

		if (variable::is_indirect) file << ", double prev_maxval";

		if (is_lambda) file << ", double& max_var";

		file << ");\n";
	}

	static void intervalPrototype(string& funcname, ofstream& file)
	{
		file << ((is_lambda || boolean_predicate) ? ("bool ") : ("int ")) << funcname << "(";

		int i;
		bool first = true;
		for (i = 1; i < (int)all_vars.size(); i++)
		{
			variable& v = all_vars[i];
			if (is_lambda && v.name.substr(0, 6) == "lambda") file << ", interval_number& " << v.name;
			if (v.isInput() && v.is_used && v.name != "1" && v.name != "2")
			{
				if (first) file << "interval_number " << v.name;
				else file << ", interval_number " << v.name;
				first = false;
			}
		}
		file << ");\n";
	}

	static void bigfloatPrototype(string& funcname, ofstream& file)
	{
		if (is_lambda) file << "void " << funcname << "(";
		else if (boolean_predicate) file << "bool " << funcname << "(";
		else file << "int " << funcname << "(";

		int i;
		bool first = true;
		for (i = 1; i < (int)all_vars.size(); i++)
		{
			variable& v = all_vars[i];
			if (is_lambda && v.name.substr(0, 6) == "lambda") file << ", bigfloat& " << v.name;
			if (v.isInput() && v.is_used && v.name != "1" && v.name != "2")
			{
				if (first) file << "bigfloat " << v.name;
				else file << ", bigfloat " << v.name;
				first = false;
			}
		}
		file << ");\n";
	}

	static void exactPrototype(string& funcname, ofstream& file)
	{
		char* return_type = (is_lambda) ? ("void ") : ("int ");
		file << return_type << funcname << "(";
		bool first = true; for (variable& v : all_vars) if (v.isInput() && v.name != "1" && v.name != "2" && v.is_used && !v.is_lambda_out) { file << ((first) ? ("double ") : (", double ")) << v.name; first = false; }

		if (is_lambda) for (variable& v : all_vars) { if (v.name.substr(0, 6) == "lambda") file << ((first) ? ("double **") : (", double **")) << v.name << ", int& " << v.name << "_len"; first = false; }
		else for (variable& v : all_vars) { if (v.isInput() && v.is_lambda_out) { file << ((first) ? ("double *") : (", double *")) << v.name << ", int& " << v.name << "_len";  first = false; } }

		for (int i = 1; i < (int)all_vars.size(); i++) if (!all_vars[i].isInput()) break;

		file << ");\n";
	}

	static void printErrorBounds()
	{
		int i;
		for (i = 1; i < (int)all_vars.size(); i++)
		{
			variable& v = all_vars[i];
			cout << v.name << " " << v.error_degree << " " << v.size << " ";
			cout << std::setprecision(std::numeric_limits<fpnumber>::digits10 + 1) << v.error_bound << " ";
			cout << std::setprecision(std::numeric_limits<fpnumber>::digits10 + 1) << v.value_bound << "\n";
		}
		cout << "NAME ERR_DEGREE EXP_SIZE ERR_BOUND VAL_BOUND\n";
	}

};

void lambda_variable::print_filtered(ofstream& file)
{
	bool first;
	file << ipoint_name << ".getFilteredLambda(";
	first = true; for (int v : output_pars) { file << ((first) ? ("") : (", ")) << variable::all_vars[v].name; first = false; }
	file << ", max_var)";
}

void lambda_variable::print_interval(ofstream& file)
{
	bool first;
	file << ipoint_name << ".getIntervalLambda(";
	first = true; for (int v : output_pars) { file << ((first) ? ("") : (", ")) << variable::all_vars[v].name; first = false; }
	file << ")";
}

void lambda_variable::print_bigfloat(ofstream& file)
{
	bool first;
	file << ipoint_name << ".getBigfloatLambda(";
	first = true; for (int v : output_pars) { file << ((first) ? ("") : (", ")) << variable::all_vars[v].name; first = false; }
	file << ")";
}

void lambda_variable::print_exact(ofstream& file)
{
	bool first;
	file << ipoint_name << ".getExactLambda(";
	first = true; for (int v : output_pars) { file << ((first) ? ("&") : (", &")) << variable::all_vars[v].name << ", " << variable::all_vars[v].name << "_len"; first = false; }
	file << ")";
}

/*
void lambda_variable::print_filtered(ofstream& file)
{
	bool first;
	file << name << "_filtered(";
	first = true; for (int v : input_pars) { file << ((first) ? ("") : (", ")) << variable::all_vars[v].name; first = false; }
	for (int v : output_pars) file << ", " << variable::all_vars[v].name;
	file << ", max_var)";
}

void lambda_variable::print_filtered_proto(ofstream& file)
{
	bool first;
	file << "bool " << name << "_filtered(";
	first = true; for (int v : input_pars) { file << ((first) ? ("fpnumber ") : (", fpnumber ")) << variable::all_vars[v].name; first = false; }
	for (int v : output_pars) file << ", fpnumber& " << variable::all_vars[v].name;
	file << ", fpnumber& max_var);\n";
}

void lambda_variable::print_interval(ofstream& file)
{
	bool first;
	file << name << "_interval(";
	first = true; for (int v : input_pars) { file << ((first) ? ("") : (", ")) << variable::all_vars[v].name; first = false; }
	for (int v : output_pars) file << ", " << variable::all_vars[v].name;
	file << ")";
}

void lambda_variable::print_interval_proto(ofstream& file)
{
	bool first;
	file << "bool " << name << "_interval(";
	first = true; for (int v : input_pars) { file << ((first) ? ("interval_number ") : (", interval_number ")) << variable::all_vars[v].name; first = false; }
	for (int v : output_pars) file << ", interval_number& " << variable::all_vars[v].name;
	file << ");\n";
}

void lambda_variable::print_exact(ofstream& file)
{
	bool first;
	file << name << "_exact(";
	first = true; for (int v : input_pars) { file << ((first) ? ("") : (", ")) << variable::all_vars[v].name; first = false; }
	for (int v : output_pars) file << ", " << variable::all_vars[v].name << ", " << variable::all_vars[v].name << "_len";
	file << ")";
}

void lambda_variable::print_exact_proto(ofstream& file)
{
	bool first;
	file << "void " << name << "_exact(";
	first = true; for (int v : input_pars) { file << ((first) ? ("fpnumber ") : (", fpnumber ")) << variable::all_vars[v].name; first = false; }
	for (int v : output_pars) file << ", fpnumber *" << variable::all_vars[v].name << ", int& " << variable::all_vars[v].name << "_len";
	file << ");\n";
}
*/
vector<variable> variable::all_vars;
vector<variable *> variable::sign_vars;
vector<lambda_variable> variable::all_lambda_vars;
string variable::func_name;
bool variable::is_indirect;
bool variable::is_lambda;
bool variable::boolean_predicate;
bool variable::is_fma;
bool variable::append;

void parseLambdaVar(string& line, int ln)
{
	size_t pos_func = line.find('(');
	string func_name = line.substr(0, pos_func); // implicit2d or implicit3d
	string parameters = line.substr(pos_func + 1, line.size() - (pos_func + 2));
	pos_func = parameters.find(':');
	string input_pars = parameters.substr(0, pos_func); // p1
	string output_pars = parameters.substr(pos_func + 1, parameters.size() - (pos_func + 1)); // l1x,l1y,l1z,d1
	//cout << func_name << endl;
	//cout << input_pars << endl;
	//cout << output_pars << endl;

	variable::all_lambda_vars.push_back(lambda_variable(func_name));
	lambda_variable& l = variable::all_lambda_vars.back();

	vector<string> tokens;
	tokenize(input_pars, tokens, ',');

	l.ipoint_name = tokens[0];

	tokenize(output_pars, tokens, ',');
	int num_outvars = tokens.size();

	for (int i = 0; i < (int)tokens.size(); i ++)
	{
		variable::all_vars.push_back(variable(true, tokens[i]));
		l.output_pars.push_back((variable::all_vars.size()-1));
	}
}

bool parseLine(ifstream& file, int ln)
{
	string line;
	vector<string> all_toks, tokens;
	if (!getline(file, line)) return false;
	tokenize(line, all_toks, ' ');

	for (unsigned int i = 0; i < all_toks.size(); i++)
		if (all_toks[i].size() > 1 && all_toks[i][0] == '/' && all_toks[i][1] == '/') break;
		else tokens.push_back(all_toks[i]);

	if (tokens.size() == 0) return true; // Skip blank lines

	if (tokens.size() == 1) // Single parameter of lambda_var
	{
		size_t pos_func = tokens[0].find('(');
		if (pos_func == string::npos) // Single parameter
		{
			variable::all_vars.push_back(variable(tokens[0]));
		}
		else // lambda_var
		{
			parseLambdaVar(tokens[0], ln);
		}
	}
	else
	{
		if (tokens[1] == "=") // Assignment
		{
			if (tokens.size() == 5) // Binary op
			{
				string& newvarname = tokens[0];
				if (variable::getVarByName(newvarname) != -1) error("Error: variable already declared", ln);
				int op1 = variable::getVarByName(tokens[2]); if (op1 == -1) error("Error: unknown variable used", ln);
				int op2 = variable::getVarByName(tokens[4]); if (op2 == -1) error("Error: unknown variable used", ln);
				char op = tokens[3][0];
				variable::all_vars.push_back(variable(newvarname, op1, op2, op));
			}
			else if (tokens.size() == 3) // Plain assignment
			{
				string& newvarname = tokens[0];
				if (variable::getVarByName(newvarname) != -1) error("Error: variable already declared", ln);
				int op1 = variable::getVarByName(tokens[2]); if (op1 == -1) error("Error: unknown variable used", ln);
				variable::all_vars.push_back(variable(newvarname, op1));
			}
			else error("Error: unexpected number of tokens", ln);
		}
		else // Parameter list
		{
			for (unsigned int i = 0; i < tokens.size(); i++) variable::all_vars.push_back(variable(tokens[i]));
		}
	}

	
	return true;
}

int main(int argc, char *argv[])
{
	if (argc < 2) error("USAGE: converter filename.txt [-p]\nIf '-p', all the code is appended to 'indirect_predicates.cpp'\n", 0);

	_controlfp(_RC_UP, _MCW_RC);

	ifstream file(argv[1]);
	if (!file.is_open()) error("Error: Can't open file", 0);

	variable::init(argc>2);

	string fullname(argv[1]);
	int sls_idx = fullname.find_last_of("\\"); if (sls_idx == string::npos) sls_idx = -1;
	int ext_idx = fullname.find_last_of(".");
	variable::func_name = fullname.substr(sls_idx + 1, ext_idx - sls_idx - 1);

	if (variable::func_name.substr(0, 6) == "lambda") variable::is_lambda = true;

	int ln = 1;

	while (parseLine(file, ln)) ln++;

	// Print all
	//variable::printErrorBounds();
	variable::produceAllCode(variable::func_name);

	return 0;
}
