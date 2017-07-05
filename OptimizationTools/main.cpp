#ifdef _WIN32
#define _SCL_SECURE_NO_WARNINGS
#endif
#include<iostream>
#include "GA_model_solver.h"

int main() {

	std::string equation = "-(100*(y-x^2)^2+(1-x)^2)";
	std::unordered_map<std::string, GA_model_solver::Gene_Traits> dec_vars;
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
	std::unordered_map<std::string, GA_model_solver::cont> params;
	GA_model_solver solver{ equation, dec_vars, params, 4, 10000 };
	solver.print_individual(solver.solve());
	
	return 0;
}