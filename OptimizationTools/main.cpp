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
		//rosenbrock's (global minimum at 1,1)
		//std::string equation = "100*(y-x^2)^2+(1-x)^2";
		//rastrigin's (global minimum at 0,0)
		//std::string equation = "20+(x^2-10*cos(2*pi*x))+(y^2-10*cos(2*pi*y))";
		
		std::string equation = "(0.05*(3*x11+2*x21+d1^2+d2^3+1-2*1) + 0.175*(3*x12+2*x22+d1^2+d2^3+3-2*1) + 0.025*(3*x13+2*x23+d1^2+d2^3+5-2*1) + 0.1*(3*x14+2*x24+d1^2+d2^3+1-2*3) + 0.35*(3*x15+2*x25+d1^2+d2^3+3-2*3) + 0.05*(3*x16+2*x26+d1^2+d2^3+5-2*3) + 0.05*(3*x17+2*x27+d1^2+d2^3+1-2*5) + 0.175*(3*x18+2*x28+d1^2+d2^3+3-2*5) + 0.025*(3*x19+2*x29+d1^2+d2^3+5-2*5))";
		std::unordered_map<std::string, PSO_solver<>::Var_Traits> dec_vars;
		dec_vars["x11"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x12"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x13"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x14"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x15"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x16"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x17"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x18"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x19"] = PSO_solver<>::Var_Traits{
			0.0, 10.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x21"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x22"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x23"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x24"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x25"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x26"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x27"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x28"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["x29"] = PSO_solver<>::Var_Traits{
			-1000.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["d1"] = PSO_solver<>::Var_Traits{
			0.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["d2"] = PSO_solver<>::Var_Traits{
			0.0, 1000.0, PSO_solver<>::Var_Traits::CONT };
		std::unordered_map<std::string, PSO_solver<>::cont_type> params;
		std::vector<std::string> constraints;
		//constraints.push_back("x <= y");
		constraints.push_back("d1 > 2*(x11+x21)-1");
		constraints.push_back("d1 > 2*(x12+x22)-1");
		constraints.push_back("d1 > 2*(x13+x23)-1");
		constraints.push_back("d1 > 2*(x14+x24)-3");
		constraints.push_back("d1 > 2*(x15+x25)-3");
		constraints.push_back("d1 > 2*(x16+x26)-3");
		constraints.push_back("d1 > 2*(x17+x27)-5");
		constraints.push_back("d1 > 2*(x18+x28)-5");
		constraints.push_back("d1 > 2*(x19+x29)-5");
		constraints.push_back("d2 > 1/x11+2*x21^2+2*1+0.5*1");
		constraints.push_back("d2 > 1/x12+2*x22^2+2*3+0.5*1");
		constraints.push_back("d2 > 1/x13+2*x23^2+2*5+0.5*1");
		constraints.push_back("d2 > 1/x14+2*x24^2+2*1+0.5*3");
		constraints.push_back("d2 > 1/x15+2*x25^2+2*3+0.5*3");
		constraints.push_back("d2 > 1/x16+2*x26^2+2*5+0.5*3");
		constraints.push_back("d2 > 1/x17+2*x27^2+2*1+0.5*5");
		constraints.push_back("d2 > 1/x18+2*x28^2+2*3+0.5*5");
		constraints.push_back("d2 > 1/x19+2*x29^2+2*5+0.5*5");
		constraints.push_back("x11^2+2*x21 < 30");
		constraints.push_back("x12^2+2*x22 < 30");
		constraints.push_back("x13^2+2*x23 < 30");
		constraints.push_back("x14^2+2*x24 < 30");
		constraints.push_back("x15^2+2*x25 < 30");
		constraints.push_back("x16^2+2*x26 < 30");
		constraints.push_back("x17^2+2*x27 < 30");
		constraints.push_back("x18^2+2*x28 < 30");
		constraints.push_back("x19^2+2*x29 < 30");
		PSO_solver<> solver{ equation, dec_vars, params, constraints, 
			0.5, 1.4, 1.4, 0.1, 1.0e-10, 0.7, 0.2, 100, 500 };
		PSO_solver<>::Hint hint;
		hint["x11"] = 0.299645;
		hint["x24"] = -1.75135;
		hint["x19"] = 4.89214;
		hint["x12"] = 0.359591;
		hint["x13"] = 0.59274;
		hint["x29"] = -0.0109772;
		hint["x14"] = 0.325731;
		hint["x28"] = -0.955697;
		hint["x15"] = 0.38244;
		hint["x16"] = 0.893842;
		hint["x17"] = 0.324952;
		hint["x18"] = 0.420533;
		hint["x25"] = -1.13795;
		hint["x21"] = -1.85302;
		hint["x22"] = -1.30838;
		hint["x23"] = -0.508464;
		hint["x26"] = -0.207214;
		hint["x27"] = -1.60109;
		hint["d1"] = 4.76333;
		hint["d2"] = 12.7057;
		//solver.hint(hint);
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