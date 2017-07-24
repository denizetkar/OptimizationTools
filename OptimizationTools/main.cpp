#ifdef _WIN32
#define _SCL_SECURE_NO_WARNINGS
#endif
//#include "GA_model_solver.h"
//#include "simplex.h"
#include "PSO_solver.h"
#include<string>
#include<iostream>
#include<fstream>

int main() {

	{
		std::string equation = "100*(y-x^2)^2+(1-x)^2";
		std::unordered_map<std::string, PSO_solver<>::Var_Traits> dec_vars;
		auto& traits = dec_vars["x"];
		traits.type = traits.CONT;
		traits.lower_bound = -5.0;
		traits.upper_bound = 5.0;
		auto& traits2 = dec_vars["y"];
		traits2.type = traits2.CONT;
		traits2.lower_bound = -5.0;
		traits2.upper_bound = 5.0;
		std::unordered_map<std::string, PSO_solver<>::cont_type> params;
		std::vector<std::string> constraints;
		//constraints.push_back("x <= y");
		PSO_solver<> solver{ equation, dec_vars, params, constraints };
		PSO_solver<>::Solution* soln = solver.solve();

		std::cout << soln << std::endl;

		delete soln;
	}

	/*{
		Simplex_solver<> solver;
		const char pathSeparator =
#if defined(WIN32) || defined(_WIN32)
			'\\';
#else
			'/';
#endif
		std::string file_path;
		file_path.append("lp_examples") += pathSeparator;
		file_path.append("lp.txt");
		std::ifstream in{ file_path };
		in >> solver;
		Simplex_solver<>::Solution * soln = solver.solve();

		for (auto itr = soln->begin(), end = soln->end(); itr != end; ++itr) {
			std::cout << *itr << " ";
		}

		delete soln;
	}*/

	/*{
		std::string equation = "-(100*(y-x^2)^2+(1-x)^2)";
		std::unordered_map<std::string, GA_model_solver<>::Gene_Traits> dec_vars;
		auto& traits = dec_vars["x"];
		traits.type = traits.CONT;
		traits.lower_bound.continuous = -5.0;
		traits.upper_bound.continuous = 5.0;
		traits.step_size.continuous = 1.0;
		auto& traits2 = dec_vars["y"];
		traits2.type = traits2.CONT;
		traits2.lower_bound.continuous = -5.0;
		traits2.upper_bound.continuous = 5.0;
		traits2.step_size.continuous = 1.0;
		std::unordered_map<std::string, GA_model_solver<>::cont_type> params;
		std::vector<std::string> constraints;
		constraints.push_back("x <= y");
		GA_model_solver<> solver{ equation, dec_vars, params, constraints };
		GA_model_solver<>::Solution* soln = solver.solve();

		std::cout << soln << std::endl;

		delete soln;
	}*/

	return 0;
}