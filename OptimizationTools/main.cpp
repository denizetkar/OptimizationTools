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
		
		std::string equation = "100*(x1+x2+x3) + 0.01*250*(95*(0.04*y11+2*y12+2.88*(1-y11-y12)) + 157*(0.08*y21+1.36*y22+1.32*(1-y21-y22)) + 46*(0.36*y31+0.08*y32+1.04*(1-y31-y32)) + 234*(0.88*y41+0.1*y42+0.52*(1-y41-y42)) + 75*(1.52*y51+1.8*y52+0.12*(1-y51-y52)) + 192*(3.36*y61+2.28*y62+0.08*(1-y61-y62))) + x1*(10*n1+0.01*(1.3+(0.24*(250*(95*y11+157*y21+46*y31+234*y41+75*y51+192*y61))/n1))*n1+(0.06*(250*(95*y11+157*y21+46*y31+234*y41+75*y51+192*y61))/n1))+x2*(10*n2+0.01*(1+(0.2*(250*(95*y12+157*y22+46*y32+234*y42+75*y52+192*y62))/n2))*n2+(0.06*(250*(95*y12+157*y22+46*y32+234*y42+75*y52+192*y62))/n2))+x3*(10*n3+0.01*(1.4+(0.28*(250*(95*(1-y11-y12)+157*(1-y21-y22)+46*(1-y31-y32)+234*(1-y41-y42)+75*(1-y51-y52)+192*(1-y61-y62)))/n3))*n3+(0.06*(250*(95*(1-y11-y12)+157*(1-y21-y22)+46*(1-y31-y32)+234*(1-y41-y42)+75*(1-y51-y52)+192*(1-y61-y62)))/n3)) + 0.2352*((7*y11*900+7*y21*2500+7*y31*625+7*y41*6400+7*y51*625+7*y61*6400)^0.5+(7*y12*900+7*y22*2500+7*y32*625+7*y42*6400+7*y52*625+7*y62*6400)^0.5+(7*(1-y11-y12)*900+7*(1-y21-y22)*2500+7*(1-y31-y32)*625+7*(1-y41-y42)*6400+7*(1-y51-y52)*625+7*(1-y61-y62)*6400)^0.5)";
		std::unordered_map<std::string, PSO_solver<>::Var_Traits> dec_vars;
		dec_vars["x1"] = PSO_solver<>::Var_Traits{ 
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["x2"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["x3"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y11"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y12"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		//dec_vars["y13"] = PSO_solver<>::Var_Traits{
		//	0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y21"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y22"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		//dec_vars["y23"] = PSO_solver<>::Var_Traits{
		//	0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y31"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y32"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		//dec_vars["y33"] = PSO_solver<>::Var_Traits{
		//	0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y41"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y42"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		//dec_vars["y43"] = PSO_solver<>::Var_Traits{
		//	0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y51"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y52"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		//dec_vars["y53"] = PSO_solver<>::Var_Traits{
		//	0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y61"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["y62"] = PSO_solver<>::Var_Traits{
			0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		//dec_vars["y63"] = PSO_solver<>::Var_Traits{
		//	0.0, 1.0, PSO_solver<>::Var_Traits::DISC };
		dec_vars["n1"] = PSO_solver<>::Var_Traits{
			0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["n2"] = PSO_solver<>::Var_Traits{
			0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		dec_vars["n3"] = PSO_solver<>::Var_Traits{
			0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		//dec_vars["d1"] = PSO_solver<>::Var_Traits{
		//	0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		//dec_vars["d2"] = PSO_solver<>::Var_Traits{
		//	0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		//dec_vars["d3"] = PSO_solver<>::Var_Traits{
		//	0.0, 10000.0, PSO_solver<>::Var_Traits::CONT };
		std::unordered_map<std::string, PSO_solver<>::cont_type> params;
		std::vector<std::string> constraints;
		//constraints.push_back("x <= y");
		//constraints.push_back("y11+y12+y13 == 1");
		//constraints.push_back("y21+y22+y23 == 1");
		//constraints.push_back("y31+y32+y33 == 1");
		//constraints.push_back("y41+y42+y43 == 1");
		//constraints.push_back("y51+y52+y53 == 1");
		//constraints.push_back("y61+y62+y63 == 1");
		constraints.push_back("y11 < x1");
		constraints.push_back("y21 < x1");
		constraints.push_back("y31 < x1");
		constraints.push_back("y41 < x1");
		constraints.push_back("y51 < x1");
		constraints.push_back("y61 < x1");
		constraints.push_back("y12 < x2");
		constraints.push_back("y22 < x2");
		constraints.push_back("y32 < x2");
		constraints.push_back("y42 < x2");
		constraints.push_back("y52 < x2");
		constraints.push_back("y62 < x2");
		constraints.push_back("(1-y11-y12) < x3");
		constraints.push_back("(1-y21-y22) < x3");
		constraints.push_back("(1-y31-y32) < x3");
		constraints.push_back("(1-y41-y42) < x3");
		constraints.push_back("(1-y51-y52) < x3");
		constraints.push_back("(1-y61-y62) < x3");
		//constraints.push_back("d1 == 250*(95*y11+157*y21+46*y31+234*y41+75*y51+192*y61)");
		//constraints.push_back("d2 == 250*(95*y12+157*y22+46*y32+234*y42+75*y52+192*y62)");
		//constraints.push_back("d3 == 250*(95*(1-y11-y12)+157*(1-y21-y22)+46*(1-y31-y32)+234*(1-y41-y42)+75*(1-y51-y52)+192*(1-y61-y62))");
		
		//EQUIVALENTLY
		//y13 = 1-y11-y12
		//y23 = 1-y21-y22
		//y33 = 1-y31-y32
		//y43 = 1-y41-y42
		//y53 = 1-y51-y52
		//y63 = 1-y61-y62
		//d1 = 250*(95*y11+157*y21+46*y31+234*y41+75*y51+192*y61)
		//d2 = 250*(95*y12+157*y22+46*y32+234*y42+75*y52+192*y62)
		//d3 = 250*(95*(1-y11-y12)+157*(1-y21-y22)+46*(1-y31-y32)+234*(1-y41-y42)+75*(1-y51-y52)+192*(1-y61-y62))
		PSO_solver<> solver{ equation, dec_vars, params, constraints, 
		0.5, 1.4, 1.4, 0.7, 0.2, 100, 1000 };
		PSO_solver<>::Hint hint;
		hint["x1"] = 0;
		hint["x2"] = 0;
		hint["x3"] = 1;
		hint["y11"] = 0;
		hint["y12"] = 1;
		hint["y21"] = 0;
		hint["y22"] = 1;
		hint["y31"] = 0;
		hint["y32"] = 0;
		hint["y41"] = 1;
		hint["y42"] = 1;
		hint["y51"] = 1;
		hint["y52"] = 0;
		hint["y61"] = 0;
		hint["y62"] = 1;
		hint["n1"] = 873.699;
		hint["n2"] = 5108.58;
		hint["n3"] = 519.966;
		solver.hint(hint);
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