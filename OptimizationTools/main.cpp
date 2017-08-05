#ifdef _WIN32
#define _SCL_SECURE_NO_WARNINGS
#endif
#include "GA_model_solver.hpp"
#include "simplex.hpp"
#include "PSO_solver.hpp"
#include<string>
#include<iostream>
#include<fstream>

int main() {

	{
		//rosenbrock's (global minimum at 1,1)
		//std::string equation = "100*(y-x^2)^2+(1-x)^2";
		//rastrigin's (global minimum at 0,0)
		//std::string equation = "20+(x^2-10*cos(2*pi*x))+(y^2-10*cos(2*pi*y))";
		
		std::string equation = "-(1.5*x1+2.5*x2+3*x3+4.5*x4)";
		std::unordered_map<std::string, PSO_solver<>::Var_Traits> dec_vars;
		dec_vars["x1"] = PSO_solver<>::Var_Traits{
			0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x2"] = PSO_solver<>::Var_Traits{
			0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x3"] = PSO_solver<>::Var_Traits{
			0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x4"] = PSO_solver<>::Var_Traits{
			0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		std::unordered_map<std::string, PSO_solver<>::cont_type> params;
		std::vector<std::string> constraints;
		//constraints.push_back("x <= y");
		constraints.push_back("2*x1+4*x2+3*x3+7*x4 <= 100000");
		constraints.push_back("3*x1+2*x2+3*x3+4*x4 <= 50000");
		constraints.push_back("2*x1+3*x2+2*x3+5*x4 <= 60000");
		PSO_solver<> solver{ equation, dec_vars, params, constraints, 100, 500 };
		//PSO_solver<>::Hint hint;
		//solver.hint(hint);
		PSO_solver<>::Solution* soln = solver.solve();

		std::cout << soln << std::endl;

		delete soln;
	}

	{
		Simplex_solver<> solver;
		const char pathSeparator =
#if defined(WIN32) || defined(_WIN32)
			'\\';
#else
			'/';
#endif
		std::string file_path;
		file_path.append("lp_examples") += pathSeparator;
		file_path.append("lp3.txt");
		std::ifstream in{ file_path };
		in >> solver;
		Simplex_solver<>::Solution * soln = solver.solve();

		for (auto itr = soln->begin(), end = soln->end(); itr != end; ++itr) {
			std::cout << *itr << " ";
		}
		std::cout << std::endl;
		delete soln;
	}

	{
		std::string equation = "(1.5*x1+2.5*x2+3*x3+4.5*x4)";
		std::unordered_map<std::string, GA_model_solver<>::Gene_Traits> dec_vars;
		auto* traits = &dec_vars["x1"];
		traits->type = traits->CONT;
		traits->lower_bound.continuous = 0.0;
		traits->upper_bound.continuous = 10000.0;
		traits->step_size.continuous = 100.0;
		traits = &dec_vars["x2"];
		traits->type = traits->CONT;
		traits->lower_bound.continuous = 0.0;
		traits->upper_bound.continuous = 10000.0;
		traits->step_size.continuous = 100.0;
		traits = &dec_vars["x3"];
		traits->type = traits->CONT;
		traits->lower_bound.continuous = 0.0;
		traits->upper_bound.continuous = 10000.0;
		traits->step_size.continuous = 100.0;
		traits = &dec_vars["x4"];
		traits->type = traits->CONT;
		traits->lower_bound.continuous = 0.0;
		traits->upper_bound.continuous = 10000.0;
		traits->step_size.continuous = 100.0;
		std::unordered_map<std::string, GA_model_solver<>::cont_type> params;
		std::vector<std::string> constraints;
		//constraints.push_back("x <= y");
		constraints.push_back("2*x1+4*x2+3*x3+7*x4 <= 100000");
		constraints.push_back("3*x1+2*x2+3*x3+4*x4 <= 50000");
		constraints.push_back("2*x1+3*x2+2*x3+5*x4 <= 60000");
		GA_model_solver<> solver{ equation, dec_vars, params, constraints, 100, 1500 };
		GA_model_solver<>::Solution* soln = solver.solve();

		std::cout << soln << std::endl;

		delete soln;
	}

	return 0;
}