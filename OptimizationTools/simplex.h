#pragma once
#include "split.h"
#include<vector>
#include<unordered_map>
#include<unordered_set>
#include<iostream>
#include<utility>
#include<limits>
#include<string>
#include<sstream>

template <typename numeric_type = double>
class Simplex_solver {
	//using numeric_type = double;
	enum Ineq_type { EQ, LTE, GTE, URS };
	struct Equation {
		std::unordered_map<size_t, numeric_type> var_to_coef;
		size_t basic_var_index;	//ONLY FIELD TO OMIT WHILE INITIALIZATION
		numeric_type rhs;
		Ineq_type ineq_type;
	};
	struct Variable {
		numeric_type val;	//ONLY FIELD TO OMIT WHILE INITIALIZATION
		enum { DEC_VAR, SLACK_VAR, SURP_VAR, ARTF_VAR, URS_VAR } var_type;
		Ineq_type ineq_type;
	};
	struct Variable_seperator {
		bool operator()(int ch) {
			return (ch == ';');
		}
	};
	struct Value_seperator {
		bool operator()(int ch) {
			return (ch == ',');
		}
	};

	std::vector<Equation> equations;
	//DEC_VAR, ... , SLACK_VAR/SURP_VAR, ... , ARTF_VAR (stored in this order!)
	std::vector<Variable> variables;
	numeric_type TOLERANCE;
	bool isInit;

	inline bool zeroCheck(const numeric_type& num) const {
		return (num > -TOLERANCE && num < TOLERANCE);
	}
public:
	using Solution = std::vector<numeric_type>;

	bool init(std::istream& in) {
		equations.clear();
		variables.clear();
		isInit = false;

		std::vector<std::string> splits;
		std::string line;
		std::istringstream ss;
		numeric_type tmp;
		size_t vsz;

		std::getline(in, line);	//FIRST WE EXPECT EITHER "sparse" OR "dense"
		if (line == "sparse") {
			if (!std::getline(in, line)) {
				return false;
			}
			//ADD OBJECTIVE VALUE 'z'
			variables.push_back({ numeric_type{}, Variable::DEC_VAR, Ineq_type::URS });
			//WE EXPECT THE SIGN RESTRICTION DECLARATIONS FOR EACH VARIABLE
			std::split(line, splits);
			for (size_t i = 0, sz = splits.size(); i < sz; ++i) {
				if (splits[i] == "LTE") {
					variables.push_back({ numeric_type{}, Variable::DEC_VAR, Ineq_type::LTE });
				}
				else if (splits[i] == "GTE") {
					variables.push_back({ numeric_type{}, Variable::DEC_VAR, Ineq_type::GTE });
				}
				else if (splits[i] == "URS") {
					variables.push_back({ numeric_type{}, Variable::DEC_VAR, Ineq_type::URS });
				}
				else {
					return false;
				}
			}
			vsz = variables.size();
			//WE EXPECT THE OBJECTIVE FUNCTION COEFFICIENTS (which is to be MAXIMIZED)
			if (!std::getline(in, line)) {
				return false;
			}
			splits.clear();
			std::split(line, splits, Variable_seperator{});
			if (splits.size() > vsz) {
				return false;
			}
			equations.push_back(Equation{});
			std::vector<std::string> sub_splits;
			size_t index;
			equations[0].var_to_coef[0] = numeric_type{ 1.0 };
			for (size_t i = 0, sz = splits.size(); i < sz; ++i) {
				std::split(splits[i], sub_splits, Value_seperator{});
				if (sub_splits.size() != 2) {
					return false;
				}
				ss.str(sub_splits[0]);
				ss.seekg(0);
				ss >> index;
				if (index >= vsz) {
					return false;
				}
				ss.str(sub_splits[1]);
				ss.seekg(0);
				ss >> tmp;
				if (tmp != 0.0) {	//ADD IF AND ONLY IF TMP IS NOT ABSOLUTELY ZERO (instead of 'zeroCheck')
					equations[0].var_to_coef[index] = tmp;
				}
				sub_splits.clear();
			}
			ss.str(splits.back());
			ss.seekg(0);
			ss >> tmp;
			equations[0].rhs = tmp;
			//READ CONSTRAINTS FROM THIS POINT ONWARDS
			for (; std::getline(in, line); ) {
				if (line.empty()) {
					continue;
				}
				splits.clear();
				std::split(line, splits, Variable_seperator{});
				//INEQ_TYPE - VARIABLE COEFFICIENTS - RHS
				//VARIABLE COEFFICIENTS: "index,value"
				if (splits.size() > (vsz + 1)) {
					return false;
				}
				equations.push_back(Equation{});
				if (splits[0] == "EQ") {
					equations.back().ineq_type = Ineq_type::EQ;
				}
				else if (splits[0] == "LTE") {
					equations.back().ineq_type = Ineq_type::LTE;
				}
				else if (splits[0] == "GTE") {
					equations.back().ineq_type = Ineq_type::GTE;
				}
				else {
					return false;
				}
				for (size_t i = 0, sz = splits.size() - 2; i < sz; ++i) {
					std::split(splits[i + 1], sub_splits, Value_seperator{});
					if (sub_splits.size() != 2) {
						return false;
					}
					ss.str(sub_splits[0]);
					ss.seekg(0);
					ss >> index;
					if (index >= vsz) {
						return false;
					}
					ss.str(sub_splits[1]);
					ss.seekg(0);
					ss >> tmp;
					if (tmp != 0.0) {	//ADD IF AND ONLY IF TMP IS NOT ABSOLUTELY ZERO (instead of 'zeroCheck')
						equations.back().var_to_coef[index] = tmp;
					}
					sub_splits.clear();
				}
				ss.str(splits.back());
				ss.seekg(0);
				ss >> tmp;
				equations.back().rhs = tmp;
			}
		}
		else if (line == "dense") {
			if (!std::getline(in, line)) {
				return false;
			}
			//ADD OBJECTIVE VALUE 'z'
			variables.push_back({ numeric_type{}, Variable::DEC_VAR, Ineq_type::URS });
			//WE EXPECT THE SIGN RESTRICTION DECLARATIONS FOR EACH VARIABLE
			std::split(line, splits);
			for (size_t i = 0, sz = splits.size(); i < sz; ++i) {
				if (splits[i] == "LTE") {
					variables.push_back({ numeric_type{}, Variable::DEC_VAR, Ineq_type::LTE });
				}
				else if (splits[i] == "GTE") {
					variables.push_back({ numeric_type{}, Variable::DEC_VAR, Ineq_type::GTE });
				}
				else if (splits[i] == "URS") {
					variables.push_back({ numeric_type{}, Variable::DEC_VAR, Ineq_type::URS });
				}
				else {
					return false;
				}
			}
			vsz = variables.size();
			//WE EXPECT THE OBJECTIVE FUNCTION COEFFICIENTS (which is to be MAXIMIZED)
			if (!std::getline(in, line)) {
				return false;
			}
			splits.clear();
			std::split(line, splits);
			if (splits.size() != vsz) {
				return false;
			}
			equations.push_back(Equation{});
			equations[0].var_to_coef[0] = numeric_type{ 1.0 };
			for (size_t i = 0, sz = vsz - 1; i < sz; ++i) {
				ss.str(splits[i]);
				ss.seekg(0);
				ss >> tmp;
				if (tmp != 0.0) {	//ADD IF AND ONLY IF TMP IS NOT ABSOLUTELY ZERO (instead of 'zeroCheck')
					equations[0].var_to_coef[i + 1] = tmp;
				}
			}
			ss.str(splits.back());
			ss.seekg(0);
			ss >> tmp;
			equations[0].rhs = tmp;
			//READ CONSTRAINTS FROM THIS POINT ONWARDS
			for (; std::getline(in, line); ) {
				if (line.empty()) {
					continue;
				}
				splits.clear();
				std::split(line, splits);
				//INEQ_TYPE - VARIABLE COEFFICIENTS - RHS
				if (splits.size() != (vsz + 1)) {
					return false;
				}
				equations.push_back(Equation{});
				if (splits[0] == "EQ") {
					equations.back().ineq_type = Ineq_type::EQ;
				}
				else if (splits[0] == "LTE") {
					equations.back().ineq_type = Ineq_type::LTE;
				}
				else if (splits[0] == "GTE") {
					equations.back().ineq_type = Ineq_type::GTE;
				}
				else {
					return false;
				}
				for (size_t i = 0, sz = vsz - 1; i < sz; ++i) {
					ss.str(splits[i + 1]);
					ss.seekg(0);
					ss >> tmp;
					if (tmp != 0.0) {	//ADD IF AND ONLY IF TMP IS NOT ABSOLUTELY ZERO (instead of 'zeroCheck')
						equations.back().var_to_coef[i + 1] = tmp;
					}
				}
				ss.str(splits.back());
				ss.seekg(0);
				ss >> tmp;
				equations.back().rhs = tmp;
			}
		}
		else {
			return false;
		}
		
		isInit = true;
		return true;
	}

	friend std::istream& operator>> (std::istream& is, Simplex_solver& obj)
	{
		if (!obj.init(is)) {
			is.setstate(std::ios::failbit);
		}
		return is;
	}

	Solution * solve() {
		if (!isInit) {
			return nullptr;
		}
		isInit = false;
		Solution * soln = new Solution;
		//MAPS URS VARIABLE 'x' TO 'x1' (same as 'x') AND 'x2' WHERE 'x = x1 - x2'
		std::unordered_map<size_t, size_t> urs_map;
		Equation two_phase_z;
		two_phase_z.var_to_coef[0] = numeric_type{ 1.0 };
		two_phase_z.basic_var_index = 0;
		two_phase_z.rhs = numeric_type{ 0.0 };
		two_phase_z.ineq_type = Ineq_type::EQ;
		equations[0].basic_var_index = 0;
		//EQUATION SIZE WILL NOT CHANGE DURING SOLVING THE LP
		size_t esz = equations.size();
		bool easyBFS = true, isOptimal, isUnbounded;
		//___CONVERT THE LP INTO STANDARD FORM___
		//CONVERT ALL VARIABLES TO GTE (greater than or equal to)
		variables[0].val = numeric_type{ 0.0 };
		for (size_t i = 1, vsz = variables.size(), j = vsz; i < vsz; ++i) {
			variables[i].val = numeric_type{ 0.0 };
			switch (variables[i].ineq_type) {
			case Ineq_type::LTE:
				for (size_t k = 0; i < esz; ++k) {
					//UPDATE THE NEGATED VARIABLE COEFFICIENTS
					auto search = equations[k].var_to_coef.find(i);
					if (search != equations[k].var_to_coef.end()) {
						search->second = -search->second;
					}
				}
				//variables[i].ineq_type = Ineq_type::GTE;
			case Ineq_type::GTE:
				break;
			case Ineq_type::URS:
				for (size_t k = 0; k < esz; ++k) {
					//UPDATE THE DECOMPOSED VARIABLE COEFFICIENTS
					auto search = equations[k].var_to_coef.find(i);
					if (search != equations[k].var_to_coef.end()) {
						equations[k].var_to_coef[j] = -search->second;
					}
				}
				variables.push_back({ numeric_type{ 0.0 }, Variable::URS_VAR, Ineq_type::GTE });
				urs_map[i] = j;
				++j;
				//variables[i].ineq_type = Ineq_type::GTE;
				break;
			default:
				throw "ERROR: wtf???";
			}
		}
		//EQUATIONS HAS TO INCLUDE THE OBJECTIVE FUNCTION (size >= 1) !!!
		//MAKE SURE THE RIGHT HAND SIDES ARE NON-NEGATIVE AND ADD NECESSARY VARIABLES
		for (size_t i = 1; i < esz; ++i) {
			if (equations[i].rhs < 0.0) {
				equations[i].rhs = -equations[i].rhs;
				for (auto itr = equations[i].var_to_coef.begin(),
					end = equations[i].var_to_coef.end(); itr != end; ++itr) {
					itr->second = -itr->second;
				}
				switch (equations[i].ineq_type) {
				case Ineq_type::LTE:
					equations[i].ineq_type = Ineq_type::GTE;
					break;
				case Ineq_type::GTE:
					equations[i].ineq_type = Ineq_type::LTE;
					break;
				default: break;
				}
			}
			switch (equations[i].ineq_type) {
			case Ineq_type::LTE:
				//ADD SLACK VARIABLE
				equations[i].var_to_coef[variables.size()] = numeric_type{ 1.0 };
				equations[i].basic_var_index = variables.size();
				variables.push_back({ equations[i].rhs, Variable::SLACK_VAR, Ineq_type::GTE });
				break;
			case Ineq_type::GTE:
				//ADD SURPLUS VARIABLE
				equations[i].var_to_coef[variables.size()] = numeric_type{ -1.0 };
				variables.push_back({ numeric_type{ 0.0 }, Variable::SURP_VAR, Ineq_type::GTE });
			case Ineq_type::EQ:
				//ADD ARTIFICIAL VARIABLE
				easyBFS = false;
				//CONSTRUCT THE OBJECTIVE FUNCTION OF 2-PHASE METHOD
				for (auto itr = equations[i].var_to_coef.begin(),
					end = equations[i].var_to_coef.end(); itr != end; ++itr) {
					numeric_type& coef = two_phase_z.var_to_coef[itr->first];
					coef -= itr->second;
					if (zeroCheck(coef)) {
						two_phase_z.var_to_coef.erase(itr->first);
					}
				}
				two_phase_z.rhs -= equations[i].rhs;
				//two_phase_z.var_to_coef[variables.size()] = numeric_type{ 1.0 };
				equations[i].var_to_coef[variables.size()] = numeric_type{ 1.0 };
				equations[i].basic_var_index = variables.size();
				variables.push_back({ equations[i].rhs, Variable::ARTF_VAR, Ineq_type::GTE });
				break;
			default:
				throw "ERROR: invalid constraint";
			}
			equations[i].ineq_type = Ineq_type::EQ;
		}
		//CONVERT THE OBJECTIVE FUNCTION INTO ROW 0 FORM (assumes variables[0] is z)
		for (auto itr = equations[0].var_to_coef.begin(),
			end = equations[0].var_to_coef.end(); itr != end; ++itr) {
			//DO NOT NEGATE COEFFICIENT OF THE 'z' VARIABLE (objective value)
			if (itr->first != 0) {
				itr->second = -itr->second;
			}
		}
		//___OBTAIN AN INITIAL BFS___
		std::vector<numeric_type> coef_ratios;
		coef_ratios.resize(esz);
		numeric_type tmp, pivot_bland_coef, two_phase_coef_ratio;
		if (!easyBFS) {
			//___EXECUTE 1ST PHASE OF THE 2-PHASE METHOD___
			while (true) {
				size_t bland_index = std::numeric_limits<size_t>::max(),
					pivot_index = std::numeric_limits<size_t>::max();
				numeric_type min_ratio = std::numeric_limits<numeric_type>::max();
				isOptimal = true;
				//SEARCH FOR BLAND INDEX (least index from among the negative coefficients)
				for (auto itr = two_phase_z.var_to_coef.begin(),
					end = two_phase_z.var_to_coef.end(); itr != end; ++itr) {
					if (itr->first < bland_index && itr->second < 0.0) {
						bland_index = itr->first;
						isOptimal = false;
					}
				}
				//IF THERE IS NO NEGATIVE COEFFICIENT THEN WE REACHED OPTIMALITY
				if (isOptimal) {
					if (!zeroCheck(two_phase_z.rhs)) {	//CASE1
						//throw "ERROR: phase-1 optimality where rhs is non-zero";
						delete soln;
						return nullptr;
					}
					break;
				}
				//SEARCH FOR PIVOT INDEX (least index from among the least ratios)
				isUnbounded = true;
				for (size_t i = 1; i < esz; ++i) {
					auto i_bland_search = equations[i].var_to_coef.find(bland_index);
					if (i_bland_search != equations[i].var_to_coef.end()) {
						if (i_bland_search->second > 0.0) {
							if (min_ratio > (tmp = equations[i].rhs / i_bland_search->second)) {
								isUnbounded = false;
								min_ratio = tmp;
								pivot_index = i;
							}
							else if (min_ratio == tmp && pivot_index > i) {
								pivot_index = i;
							}
						}
					}
				}
				if (isUnbounded) {
					throw "ERROR: the LP is unbounded (impossible here)";
				}
				//CALCULATE COEFFICIENT RATIOS (-a/b)
				//DELETE BLAND VARIABLE FROM ALL ROWS (in advance) THEN RE-ADD IT TO PIVOT ROW
				auto obj_bland_search = two_phase_z.var_to_coef.find(bland_index);
				auto pivot_bland_search = equations[pivot_index].var_to_coef.find(bland_index);
				pivot_bland_coef = pivot_bland_search->second;
				two_phase_coef_ratio = -(obj_bland_search->second) / pivot_bland_coef;
				two_phase_z.var_to_coef.erase(obj_bland_search);
				for (size_t i = 0; i < esz; ++i) {
					auto i_bland_search = equations[i].var_to_coef.find(bland_index);
					if (i_bland_search != equations[i].var_to_coef.end()) {
						coef_ratios[i] = -(i_bland_search->second) / pivot_bland_coef;
						equations[i].var_to_coef.erase(i_bland_search);
						//NOT SAFE TO USE 'pivot_bland_search' AFTER THIS POINT !!!
					}
					else {
						coef_ratios[i] = numeric_type{ 0.0 };
					}
				}
				//VARIABLE AT 'equations[pivot_index].basic_var_index' BECAME NON-BASIC
				variables[equations[pivot_index].basic_var_index].val = numeric_type{ 0.0 };
				//VARIABLE AT 'bland_index' BECAME BASIC VARIABLE OF PIVOT ROW
				equations[pivot_index].basic_var_index = bland_index;
				//___PERFORM ERO's IN ORDER TO PIVOT THE BLAND INDEX___
				//MAKE CHANGES TO OTHER ROWS FOR EACH VARIABLE (with non-zero coefficient) IN PIVOT ROW
				for (auto itr = equations[pivot_index].var_to_coef.begin(),
					end = equations[pivot_index].var_to_coef.end(); itr != end; ++itr) {
					numeric_type& coef = two_phase_z.var_to_coef[itr->first];
					coef += two_phase_coef_ratio * (itr->second);
					if (zeroCheck(coef)) {
						two_phase_z.var_to_coef.erase(itr->first);
					}
					for (size_t i = 0; i < esz; ++i) {
						if (i != pivot_index) {
							numeric_type& _coef = equations[i].var_to_coef[itr->first];
							_coef += coef_ratios[i] * (itr->second);
							if (zeroCheck(_coef)) {
								equations[i].var_to_coef.erase(itr->first);
							}
						}
					}
					itr->second /= pivot_bland_coef;
				}
				two_phase_z.rhs += two_phase_coef_ratio * equations[pivot_index].rhs;
				if (zeroCheck(two_phase_z.rhs)) {
					two_phase_z.rhs = numeric_type{ 0.0 };
				}
				for (size_t i = 0; i < esz; ++i) {
					if (i != pivot_index) {
						equations[i].rhs += coef_ratios[i] * equations[pivot_index].rhs;
						if (zeroCheck(equations[i].rhs)) {
							equations[i].rhs = numeric_type{ 0.0 };
						}
					}
				}
				equations[pivot_index].rhs /= pivot_bland_coef;
				//ADD BLAND VARIABLE BACK TO PIVOT ROW
				equations[pivot_index].var_to_coef[bland_index] = numeric_type{ 1.0 };
			}
			//BASIC VARIABLES OF ROWS ARE EQUAL TO THE RHS OF THAT ROW
			for (size_t i = 0; i < esz; ++i) {
				variables[equations[i].basic_var_index].val = equations[i].rhs;
			}
			//EITHER CASE2 OR CASE3 (since it didn't return)
			bool isCase2 = true;
			for (size_t i = 0; i < esz; ++i) {
				if (variables[equations[i].basic_var_index].var_type == Variable::ARTF_VAR) {
					isCase2 = false;
					artf_var_indexes.erase(equations[i].basic_var_index);
				}
			}
			//IN EITHER CASE WE DROP THE NON-BASIC ARTIFICIAL VARIABLES FROM ALL ROWS
			for (size_t i = 0; i < esz; ++i) {
				for (auto itr = artf_var_indexes.begin(),
					end = artf_var_indexes.end(); itr != end; ++itr) {
					equations[i].var_to_coef.erase(*itr);
				}
			}
			if (!isCase2) {	//DROP DECISION VARIABLES WITH POSITIVE COEFFICIENT IN 'two_phase_z'
				for (auto itr = two_phase_z.var_to_coef.begin(),
					end = two_phase_z.var_to_coef.end(); itr != end; ++itr) {
					switch (variables[itr->first].var_type) {
					case Variable::DEC_VAR:
					case Variable::URS_VAR:
						if (itr->second > 0.0) {
							equations[0].var_to_coef.erase(itr->first);
							variables[itr->first].val = numeric_type{ 0.0 };
						}
					default:
						break;
					}
				}
			}
		}

		while (true) {
			size_t bland_index = std::numeric_limits<size_t>::max(),
				pivot_index = std::numeric_limits<size_t>::max();
			numeric_type min_ratio = std::numeric_limits<numeric_type>::max();
			isOptimal = true;
			//SEARCH FOR BLAND INDEX (least index from among the negative coefficients)
			for (auto itr = equations[0].var_to_coef.begin(),
				end = equations[0].var_to_coef.end(); itr != end; ++itr) {
				if (itr->first < bland_index && itr->second < 0.0) {
					bland_index = itr->first;
					isOptimal = false;
				}
			}
			//IF THERE IS NO NEGATIVE COEFFICIENT THEN WE REACHED OPTIMALITY
			if (isOptimal) {
				break;
			}
			//SEARCH FOR PIVOT INDEX (least index from among the least ratios)
			isUnbounded = true;
			for (size_t i = 1; i < esz; ++i) {
				auto i_bland_search = equations[i].var_to_coef.find(bland_index);
				if (i_bland_search != equations[i].var_to_coef.end()) {
					if (i_bland_search->second > 0.0) {
						if (min_ratio > (tmp = equations[i].rhs / i_bland_search->second)) {
							isUnbounded = false;
							min_ratio = tmp;
							pivot_index = i;
						}
						else if (min_ratio == tmp && pivot_index > i) {
							pivot_index = i;
						}
					}
				}
			}
			if (isUnbounded) {
				throw "ERROR: the LP is unbounded";
			}
			//CALCULATE COEFFICIENT RATIOS (-a/b)
			//DELETE BLAND VARIABLE FROM ALL ROWS (in advance) THEN RE-ADD IT TO PIVOT ROW
			auto pivot_bland_search = equations[pivot_index].var_to_coef.find(bland_index);
			pivot_bland_coef = pivot_bland_search->second;
			for (size_t i = 0; i < esz; ++i) {
				auto i_bland_search = equations[i].var_to_coef.find(bland_index);
				if (i_bland_search != equations[i].var_to_coef.end()) {
					coef_ratios[i] = -(i_bland_search->second) / pivot_bland_coef;
					equations[i].var_to_coef.erase(i_bland_search);
					//NOT SAFE TO USE 'pivot_bland_search' AFTER THIS POINT !!!
				}
				else {
					coef_ratios[i] = numeric_type{ 0.0 };
				}
			}
			//VARIABLE AT 'equations[pivot_index].basic_var_index' BECAME NON-BASIC
			variables[equations[pivot_index].basic_var_index].val = numeric_type{ 0.0 };
			//VARIABLE AT 'bland_index' BECAME BASIC VARIABLE OF PIVOT ROW
			equations[pivot_index].basic_var_index = bland_index;
			//___PERFORM ERO's IN ORDER TO PIVOT THE BLAND INDEX___
			//MAKE CHANGES TO OTHER ROWS FOR EACH VARIABLE (with non-zero coefficient) IN PIVOT ROW
			for (auto itr = equations[pivot_index].var_to_coef.begin(),
				end = equations[pivot_index].var_to_coef.end(); itr != end; ++itr) {
				for (size_t i = 0; i < esz; ++i) {
					if (i != pivot_index) {
						numeric_type& _coef = equations[i].var_to_coef[itr->first];
						_coef += coef_ratios[i] * (itr->second);
						if (zeroCheck(_coef)) {
							equations[i].var_to_coef.erase(itr->first);
						}
					}
				}
				itr->second /= pivot_bland_coef;
			}
			for (size_t i = 0; i < esz; ++i) {
				if (i != pivot_index) {
					equations[i].rhs += coef_ratios[i] * equations[pivot_index].rhs;
					if (zeroCheck(equations[i].rhs)) {
						equations[i].rhs = numeric_type{ 0.0 };
					}
				}
			}
			equations[pivot_index].rhs /= pivot_bland_coef;
			//ADD BLAND VARIABLE BACK TO PIVOT ROW
			equations[pivot_index].var_to_coef[bland_index] = numeric_type{ 1.0 };
		}
		//BASIC VARIABLES OF ROWS ARE EQUAL TO THE RHS OF THAT ROW
		for (size_t i = 0; i < esz; ++i) {
			variables[equations[i].basic_var_index].val = equations[i].rhs;
		}

		//COPY VALUES OF THE DECISION VARIABLES AND RETURN IT AS RESULT
		variables[0].ineq_type = Ineq_type::GTE;
		for (size_t i = 0, vsz = variables.size(); i < vsz; ++i) {
			if (variables[i].var_type == Variable::URS_VAR) {
				break;
			}
			else if (variables[i].var_type == Variable::DEC_VAR) {
				switch (variables[i].ineq_type) {
				case Ineq_type::URS:
					variables[i].val -= variables[urs_map[i]].val;
					break;
				case Ineq_type::LTE:
					variables[i].val = -variables[i].val;
					break;
				default:
					break;
				}
				soln->push_back(variables[i].val);
			}
		}
		return soln;
	}

	Simplex_solver(numeric_type tolerance = numeric_type{ 0.00001 }) : isInit{ false }, TOLERANCE{ tolerance } {
		if (TOLERANCE < 0.0) {
			TOLERANCE = -TOLERANCE;
		}
	}
};
