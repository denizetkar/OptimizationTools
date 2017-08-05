
******************************************************************
******************************************************************
*
* DOWNLOADED FROM MINLP CYBER-INFRASTRUCTURE 
* www.minlp.org
*
* PROBLEM : Mixed-Integer Nonlinear Programming Models and Algorithms for Supply Chain Design with Stochastic Inventory Management
*
* AUTHOR(S) : Fengqi You, Ignacio Grossmann
*
* SUBMITTED BY : Fengqi You

******************************************************************
******************************************************************* MODEL P1
*"Mixed-Integer Nonlinear Programming Models and Algorithms for Large-Scale
* Supply Chain Design with Stochastic Inventory Management"
* Fengqi You* and Ignacio E. Grossmann**, Department of Chemical Engineering,
* Carnegie Mellon University, Pittsburgh, PA 15213
* *e-mail: yfq@andrew.cmu.edu, **e-mail: grossmann@cmu.edu
* NOTE: nomenclature missing
	std::string equation = "100*(x1+x2+x3) + 0.01*250*(95*(0.04*y11+2*y12+2.88*(1-y11-y12)) + 157*(0.08*y21+1.36*y22+1.32*(1-y21-y22)) + 46*(0.36*y31+0.08*y32+1.04*(1-y31-y32)) + 234*(0.88*y41+0.1*y42+0.52*(1-y41-y42)) + 75*(1.52*y51+1.8*y52+0.12*(1-y51-y52)) + 192*(3.36*y61+2.28*y62+0.08*(1-y61-y62))) + x1*(10*n1+0.01*(1.3+(0.24*(250*(95*y11+157*y21+46*y31+234*y41+75*y51+192*y61))/n1))*n1+(0.06*(250*(95*y11+157*y21+46*y31+234*y41+75*y51+192*y61))/n1))+x2*(10*n2+0.01*(1+(0.2*(250*(95*y12+157*y22+46*y32+234*y42+75*y52+192*y62))/n2))*n2+(0.06*(250*(95*y12+157*y22+46*y32+234*y42+75*y52+192*y62))/n2))+x3*(10*n3+0.01*(1.4+(0.28*(250*(95*(1-y11-y12)+157*(1-y21-y22)+46*(1-y31-y32)+234*(1-y41-y42)+75*(1-y51-y52)+192*(1-y61-y62)))/n3))*n3+(0.06*(250*(95*(1-y11-y12)+157*(1-y21-y22)+46*(1-y31-y32)+234*(1-y41-y42)+75*(1-y51-y52)+192*(1-y61-y62)))/n3)) + 0.2352*((7*y11*900+7*y21*2500+7*y31*625+7*y41*6400+7*y51*625+7*y61*6400)^0.5+(7*y12*900+7*y22*2500+7*y32*625+7*y42*6400+7*y52*625+7*y62*6400)^0.5+(7*(1-y11-y12)*900+7*(1-y21-y22)*2500+7*(1-y31-y32)*625+7*(1-y41-y42)*6400+7*(1-y51-y52)*625+7*(1-y61-y62)*6400)^0.5)";
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
* x1: 1 y11: 1 x2: 1 y12: 0 y62: 0 x3: 1 y31: 0 y21: 1 y22: 0 y32: 1 y61: 0 y41: 0
*  y42: 1 y51: 0 y52: 0 n1: 19.4296 n2: 20.4837 n3: 19.9985
* obj_val: 2287.91