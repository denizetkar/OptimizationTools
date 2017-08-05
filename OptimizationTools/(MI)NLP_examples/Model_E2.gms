
******************************************************************
******************************************************************
*
* DOWNLOADED FROM MINLP CYBER-INFRASTRUCTURE 
* www.minlp.org
*
* PROBLEM : A Novel Bi-Level Optimization Method for the Identification of Critical Points in Flow Sheet Synthesis under Uncertainty
*
* AUTHOR(S) : Zorka Novak Pintaric, Zdravko Kravanja
*
* SUBMITTED BY : Zorka Novak Pintaric

******************************************************************
******************************************************************* Model E2
	std::string equation = "(0.05*(3*x11+2*x21+d1^2+d2^3+1-2*1) + 0.175*(3*x12+2*x22+d1^2+d2^3+3-2*1) + 0.025*(3*x13+2*x23+d1^2+d2^3+5-2*1) + 0.1*(3*x14+2*x24+d1^2+d2^3+1-2*3) + 0.35*(3*x15+2*x25+d1^2+d2^3+3-2*3) + 0.05*(3*x16+2*x26+d1^2+d2^3+5-2*3) + 0.05*(3*x17+2*x27+d1^2+d2^3+1-2*5) + 0.175*(3*x18+2*x28+d1^2+d2^3+3-2*5) + 0.025*(3*x19+2*x29+d1^2+d2^3+5-2*5))";
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
* x11: 0.299645 x24: -1.75135 x19: 4.89214 x12: 0.359591 x13: 0.59274 x29: -0.0109
* 772 x14: 0.325731 x28: -0.955697 x15: 0.38244 x16: 0.893842 x17: 0.324952 x18: 0
* .420533 x25: -1.13795 x21: -1.85302 x22: -1.30838 x23: -0.508464 x26: -0.207214
* x27: -1.60109 d1: 4.76333 d2: 12.7057
* obj_val: 2069.82