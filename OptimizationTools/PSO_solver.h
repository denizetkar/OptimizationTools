#pragma once
#include "exprtk.hpp"
#include "prime/primesieve.hpp"
#include<unordered_map>
#include<string>
#include<vector>
#include<random>
#include<iostream>
#include<limits>
#include<algorithm>
#include<cmath>

static std::random_device rd;
static std::mt19937_64 gen(rd());
static std::uniform_real_distribution<double> unif{ 0.0, 1.0 };
static std::normal_distribution<double> norm{};

//TRIES TO MINIMIZE THE OBJECTIVE FUNCTION AND SATISFY ALL CONSTRAINTS
template <typename numeric_type = double, typename discrete_type = long long>
class PSO_solver {
protected:
	//using numeric_type = double;
	//using discrete_type = long long;
	struct Particle {
		std::vector<numeric_type> variables;
		std::vector<numeric_type> prev_vars;
		std::vector<numeric_type> velocity;
		std::vector<numeric_type> pbest;
		numeric_type pbest_constr_viol;
		numeric_type pbest_obj_val;
		numeric_type constr_viol;
		numeric_type obj_val;
	};

public:
	struct Var_Traits {
		numeric_type lower_bound, upper_bound;
		enum Type { DISC, CONT } type;
	};

protected:
	struct Var_Traits_Internal : Var_Traits {
		std::string name;
		numeric_type max_velocity;
	};
	typedef exprtk::symbol_table<numeric_type> symbol_table_t;
	typedef exprtk::expression<numeric_type>     expression_t;
	typedef exprtk::parser<numeric_type>             parser_t;
	typedef exprtk::parser_error::type error_t;

	symbol_table_t symbol_table;
	std::vector<expression_t> constraint_expressions;
	expression_t expression;
	parser_t parser;
	std::vector<Particle*> particles;
	std::vector<numeric_type> dec_var_values;
	std::vector<Var_Traits_Internal> dec_var_traits;
	size_t population_coef, max_iteration;
	//CONSTRUCT ITERATION DEPENDENT var_min AND var_max (only for discrete variables(deprecated))
	std::vector<numeric_type> var_min, var_max;
	numeric_type inertial_weight, local_accelerator, global_accelerator,
		cont_div_coef_not, cont_div_coef_min, cont_div_metrix_exp,
		disc_div_coef_not, lambda, div_metric_const, l_param;
	//other parameters !!!!!!!
public:
	using cont_type = numeric_type;
	using disc_type = discrete_type;
	using Hint = std::unordered_map<std::string, numeric_type>;
	struct Solution {
		std::unordered_map<std::string, numeric_type> solution;
		numeric_type obj_val;
		friend std::ostream& operator<<(std::ostream& out, const Solution* soln) {
			if (soln) {
				for (auto itr = soln->solution.begin(),
					end = soln->solution.end(); itr != end; ++itr) {
					out << itr->first << ": " << itr->second << " ";
				}
				out << std::endl << "obj_val: " << soln->obj_val;
			}
			return out;
		}
	};

	void hint(const Hint& _hint) {
		Particle* res = new Particle;
		size_t vsz = dec_var_traits.size();
		bool hintFound = false;
		res->variables.resize(vsz);
		res->prev_vars.resize(vsz);
		res->velocity.resize(vsz);
		res->pbest.resize(vsz);
		res->pbest_constr_viol = std::numeric_limits<numeric_type>::max();
		res->pbest_obj_val = std::numeric_limits<numeric_type>::max();
		for (size_t i = 0; i < vsz; ++i) {
			numeric_type lb = dec_var_traits[i].lower_bound;
			numeric_type ub = dec_var_traits[i].upper_bound;
			res->prev_vars[i] = lb + (ub - lb) * static_cast<numeric_type>(unif(gen));
			res->velocity[i] = dec_var_traits[i].max_velocity *
				static_cast<numeric_type>(2.0 * unif(gen) - 1.0);
			auto var_search = _hint.find(dec_var_traits[i].name);
			if (var_search != _hint.end()) {
				hintFound = true;
				res->variables[i] = var_search->second;
			}
			else {
				res->variables[i] = lb + (ub - lb) * static_cast<numeric_type>(unif(gen));
			}
		}
		if (hintFound) {
			if (particles.empty()) {
				particles.push_back(res);
			}
			else {
				if (particles[0] != nullptr)
					delete particles[0];
				particles[0] = res;
			}
			return;
		}
		//IF REACHED HERE THEN NO HINT EXISTS
		delete res;
	}

	Particle* get_particle(bool isRandom = false) {
		Particle* res = new Particle;
		size_t vsz = dec_var_values.size();
		res->variables.resize(vsz);
		res->prev_vars.resize(vsz);
		res->velocity.resize(vsz);
		res->pbest.resize(vsz);
		res->pbest_constr_viol = std::numeric_limits<numeric_type>::max();
		res->pbest_obj_val = std::numeric_limits<numeric_type>::max();
		if (isRandom) {
			for (size_t i = 0; i < vsz; ++i) {
				numeric_type lb = dec_var_traits[i].lower_bound;
				numeric_type ub = dec_var_traits[i].upper_bound;
				res->variables[i] = lb + (ub - lb) * static_cast<numeric_type>(unif(gen));
				res->prev_vars[i] = lb + (ub - lb) * static_cast<numeric_type>(unif(gen));
				res->velocity[i] = dec_var_traits[i].max_velocity * static_cast<numeric_type>(2.0 * unif(gen) - 1.0);
			}
		}
		return res;
	}

	numeric_type get_phi(size_t prime, size_t k) {
		size_t prime_ = prime, k_ = k;
		numeric_type phi = numeric_type{ 0.0 };
		while (k_ > 0) {
			phi += static_cast<numeric_type>(k_ % prime) / 
				static_cast<numeric_type>(prime_);
			k_ = k_ / prime;
			prime_ *= prime;
		}
		return phi;
	}

	Solution* solve() {
		size_t vsz = dec_var_values.size();
		var_min.resize(vsz);
		var_max.resize(vsz);
		size_t psz = population_coef * vsz;
		particles.resize(psz, nullptr);
		size_t gbest = 0;
		std::vector<numeric_type> disc_div_coef;
		std::vector<size_t> primes;
		disc_div_coef.resize(vsz);
		numeric_type adjusted_var_min, adjusted_var_max, 
			cont_div_coef, _psz = static_cast<numeric_type>(psz);
		//INITIALIZE THE PARTICLES RANDOMLY
		primesieve::generate_n_primes(vsz - 1, &primes);
		for (size_t i = (particles[0]==nullptr?0:1); i < psz; ++i) {
			Particle* res = new Particle;
			res->variables.resize(vsz);
			res->prev_vars.resize(vsz);
			res->velocity.resize(vsz);
			res->pbest.resize(vsz);
			res->pbest_constr_viol = std::numeric_limits<numeric_type>::max();
			res->pbest_obj_val = std::numeric_limits<numeric_type>::max();
			numeric_type lb = dec_var_traits[0].lower_bound;
			numeric_type ub = dec_var_traits[0].upper_bound;
			res->variables[0] = lb + (ub - lb) * 
				static_cast<numeric_type>(i) / _psz;
			res->prev_vars[0] = lb + (ub - lb) * static_cast<numeric_type>(unif(gen));
			res->velocity[0] = dec_var_traits[0].max_velocity * 
				static_cast<numeric_type>(2.0 * unif(gen) - 1.0);
			for (size_t j = 1; j < vsz; ++j) {
				numeric_type lb = dec_var_traits[j].lower_bound;
				numeric_type ub = dec_var_traits[j].upper_bound;
				res->variables[j] = lb + (ub - lb) * get_phi(primes[j - 1], i);
				res->prev_vars[j] = lb + (ub - lb) * static_cast<numeric_type>(unif(gen));
				res->velocity[j] = dec_var_traits[j].max_velocity * 
					static_cast<numeric_type>(2.0 * unif(gen) - 1.0);
			}
			particles[i] = res;
		}
		//COMPUTE FITNESSES, DETERMINE pbest AND gbest, CONSTRAINT VIOLATIONS
		for (auto itr = var_min.begin(), end = var_min.end(); itr != end; ++itr) {
			*itr = std::numeric_limits<numeric_type>::max();
		}
		for (auto itr = var_max.begin(), end = var_max.end(); itr != end; ++itr) {
			*itr = -std::numeric_limits<numeric_type>::max();
		}
		for (size_t i = 0; i < psz; ++i) {
			for (size_t j = 0; j < vsz; ++j) {
				numeric_type& var_min_j = var_min[j];
				numeric_type& var_max_j = var_max[j];
				var_min_j = std::min(var_min_j, particles[i]->variables[j]);
				var_max_j = std::max(var_max_j, particles[i]->variables[j]);
				if (dec_var_traits[j].type == Var_Traits::DISC) {
					particles[i]->variables[j] = std::round(particles[i]->variables[j]);
				}
				dec_var_values[j] = particles[i]->variables[j];
			}
			if (std::isnan(particles[i]->obj_val = expression.value())) {
				particles[i]->obj_val = std::numeric_limits<numeric_type>::max();
			}
			//CALCULATE CONSTRAINT VIOLATION
			particles[i]->constr_viol = numeric_type{ 0.0 };
			for (size_t j = 0, consz = constraint_expressions.size(); j < consz; ++j) {
				numeric_type cons_val = constraint_expressions[j].value();
				if (std::isnan(cons_val)) {
					cons_val = std::numeric_limits<numeric_type>::max();
				}
				if (cons_val > 0.0) {	//CONSTRAINT VIOLATED!!!
					particles[i]->constr_viol += cons_val;
				}
			}
			//UPDATE pbest IF NECESSARY
			bool updatePbest = false;
			if (particles[i]->constr_viol > 0.0) {	//CURRENT PARTICLE IS INFEASIBLE
				if (particles[i]->pbest_constr_viol > 0.0) {	//pbest PARTICLE IS INFEASIBLE
					if (particles[i]->constr_viol < particles[i]->pbest_constr_viol) {
						updatePbest = true;
					}
					else if (particles[i]->constr_viol == particles[i]->pbest_constr_viol) {
						if (particles[i]->pbest_obj_val > particles[i]->obj_val) {
							updatePbest = true;
						}
					}
				}
			}
			else {	//CURRENT PARTICLE IS FEASIBLE
				if (particles[i]->pbest_constr_viol > 0.0) {	//pbest PARTICLE IS INFEASIBLE
					updatePbest = true;
				}
				else {	//pbest PARTICLE IS FEASIBLE
					if (particles[i]->pbest_obj_val > particles[i]->obj_val) {
						updatePbest = true;
					}
				}
			}
			if (updatePbest) {
				for (size_t j = 0; j < vsz; ++j) {
					particles[i]->pbest[j] = particles[i]->variables[j];
				}
				particles[i]->pbest_constr_viol = particles[i]->constr_viol;
				particles[i]->pbest_obj_val = particles[i]->obj_val;
			}
			//UPDATE gbest IF NECESSARY
			bool updateGbest = false;
			if (particles[i]->constr_viol > 0.0) {	//CURRENT PARTICLE IS INFEASIBLE
				if (particles[gbest]->pbest_constr_viol > 0.0) {	//gbest PARTICLE IS INFEASIBLE
					if (particles[i]->constr_viol < particles[gbest]->pbest_constr_viol) {
						updateGbest = true;
					}
					else if (particles[i]->constr_viol == particles[gbest]->pbest_constr_viol) {
						if (particles[gbest]->pbest_obj_val > particles[i]->obj_val) {
							updateGbest = true;
						}
					}
				}
			}
			else {	//CURRENT PARTICLE IS FEASIBLE
				if (particles[gbest]->pbest_constr_viol > 0.0) {	//gbest PARTICLE IS INFEASIBLE
					updateGbest = true;
				}
				else {	//gbest PARTICLE IS FEASIBLE
					if (particles[gbest]->pbest_obj_val > particles[i]->obj_val) {
						updateGbest = true;
					}
				}
			}
			if (updateGbest) {
				gbest = i;
			}
		}
		//___START THE MAIN OPTIMIZATION LOOP___
		for (size_t iteration = 0; iteration < max_iteration; ++iteration) {
			//CALCULATE DIVERSITY METRICS (only for discrete variables)
			cont_div_coef = numeric_type{ 1.0 };
			for (size_t i = 0; i < vsz; ++i) {
				numeric_type& var_min_i = var_min[i];
				numeric_type& var_max_i = var_max[i];
				adjusted_var_min =
					std::min(
						var_max_i - lambda * (var_max_i - var_min_i),
						std::max(
							particles[gbest]->pbest[i] - 0.5 * lambda * (var_max_i - var_min_i),
							var_min_i
						)
					);
				adjusted_var_max =
					std::max(
						var_min_i + lambda * (var_max_i - var_min_i),
						std::min(
							particles[gbest]->pbest[i] + 0.5 * lambda * (var_max_i - var_min_i),
							var_max_i
						)
					);
				switch (dec_var_traits[i].type) {
				case Var_Traits::DISC:
					disc_div_coef[i] =
						disc_div_coef_not *
						std::pow(
						(dec_var_traits[i].upper_bound - dec_var_traits[i].lower_bound + 1.0),
							-std::pow(
								div_metric_const * (adjusted_var_max - adjusted_var_min) /
								(dec_var_traits[i].upper_bound - dec_var_traits[i].lower_bound),
								2.0
							)
						);
					break;
				case Var_Traits::CONT:
					cont_div_coef *= std::pow((adjusted_var_max - adjusted_var_min) /
						(dec_var_traits[i].upper_bound - dec_var_traits[i].lower_bound),
						cont_div_metrix_exp);
					break;
				default: break;
				}
			}
			cont_div_coef = cont_div_coef_not * 
				std::pow(cont_div_coef_min,
					std::pow(cont_div_coef*div_metric_const, 2)
				);
			//UPDATE THE PARTICLE VELOCITIES
			for (size_t i = 0; i < psz; ++i) {
				for (size_t j = 0; j < vsz; ++j) {
					//BETTER RESULTS WITHOUT LcRiPSO
					numeric_type sigma_kt = std::abs(particles[i]->pbest[j] - particles[i]->variables[j]);
					if (sigma_kt == 0.0) {
						sigma_kt = std::abs(particles[i]->pbest[j] - particles[i]->prev_vars[j]);
					}
					sigma_kt *= l_param;
					numeric_type f_1t =
						static_cast<numeric_type>(norm(gen)) * sigma_kt + particles[i]->pbest[j];
					numeric_type f_2t =
						static_cast<numeric_type>(norm(gen)) * sigma_kt + particles[gbest]->pbest[j];
					particles[i]->velocity[j] = std::min(
						inertial_weight * particles[i]->velocity[j] +
						local_accelerator * static_cast<numeric_type>(unif(gen)) *
						(f_1t - particles[i]->variables[j]) +
						global_accelerator * static_cast<numeric_type>(unif(gen)) *
						(f_2t - particles[i]->variables[j]), 
						dec_var_traits[j].max_velocity) +
						cont_div_coef * static_cast<numeric_type>(unif(gen)) *
						particles[i]->velocity[j];
					//UPDATE PARTICLE POSITIONS IN THE MEANTIME
					particles[i]->prev_vars[j] = particles[i]->variables[j];
					if (particles[i]->velocity[j] > (dec_var_traits[j].upper_bound - particles[i]->variables[j])) {
						particles[i]->variables[j] = dec_var_traits[j].upper_bound;
					}
					else if (particles[i]->velocity[j] < (dec_var_traits[j].lower_bound - particles[i]->variables[j])) {
						particles[i]->variables[j] = dec_var_traits[j].lower_bound;
					}
					else {
						particles[i]->variables[j] += particles[i]->velocity[j];
					}
					if (dec_var_traits[j].type == Var_Traits::DISC) {
						if (static_cast<numeric_type>(unif(gen)) <= disc_div_coef[j]) {
							if (unif(gen) <= 0.5) {
								particles[i]->variables[j] = std::floor(particles[i]->variables[j]);
							}
							else {
								particles[i]->variables[j] = std::ceil(particles[i]->variables[j]);
							}
						}
						else {
							particles[i]->variables[j] = std::round(particles[i]->variables[j]);
						}
					}
				}
			}
			//COMPUTE FITNESSES, DETERMINE pbest AND gbest, CONSTRAINT VIOLATIONS
			for (auto itr = var_min.begin(), end = var_min.end(); itr != end; ++itr) {
				*itr = std::numeric_limits<numeric_type>::max();
			}
			for (auto itr = var_max.begin(), end = var_max.end(); itr != end; ++itr) {
				*itr = -std::numeric_limits<numeric_type>::max();
			}
			for (size_t i = 0; i < psz; ++i) {
				for (size_t j = 0; j < vsz; ++j) {
					dec_var_values[j] = particles[i]->variables[j];
					if (dec_var_traits[j].type == Var_Traits::DISC) {
						numeric_type& var_min_j = var_min[j];
						numeric_type& var_max_j = var_max[j];
						var_min_j = std::min(var_min_j, particles[i]->variables[j]);
						var_max_j = std::max(var_max_j, particles[i]->variables[j]);
					}
				}
				if (std::isnan(particles[i]->obj_val = expression.value())) {
					particles[i]->obj_val = std::numeric_limits<numeric_type>::max();
				}
				//CALCULATE CONSTRAINT VIOLATION
				particles[i]->constr_viol = numeric_type{ 0.0 };
				for (size_t j = 0, consz = constraint_expressions.size(); j < consz; ++j) {
					numeric_type cons_val = constraint_expressions[j].value();
					if (std::isnan(cons_val)) {
						cons_val = std::numeric_limits<numeric_type>::max();
					}
					if (cons_val > 0.0) {	//CONSTRAINT VIOLATED!!!
						particles[i]->constr_viol += cons_val;
					}
				}
				//UPDATE pbest IF NECESSARY
				bool updatePbest = false;
				if (particles[i]->constr_viol > 0.0) {	//CURRENT PARTICLE IS INFEASIBLE
					if (particles[i]->pbest_constr_viol > 0.0) {	//pbest PARTICLE IS INFEASIBLE
						if (particles[i]->constr_viol < particles[i]->pbest_constr_viol) {
							updatePbest = true;
						}
						else if (particles[i]->constr_viol == particles[i]->pbest_constr_viol) {
							if (particles[i]->pbest_obj_val > particles[i]->obj_val) {
								updatePbest = true;
							}
						}
					}
				}
				else {	//CURRENT PARTICLE IS FEASIBLE
					if (particles[i]->pbest_constr_viol > 0.0) {	//pbest PARTICLE IS INFEASIBLE
						updatePbest = true;
					}
					else {	//pbest PARTICLE IS FEASIBLE
						if (particles[i]->pbest_obj_val > particles[i]->obj_val) {
							updatePbest = true;
						}
					}
				}
				if (updatePbest) {
					for (size_t j = 0; j < vsz; ++j) {
						particles[i]->pbest[j] = particles[i]->variables[j];
					}
					particles[i]->pbest_constr_viol = particles[i]->constr_viol;
					particles[i]->pbest_obj_val = particles[i]->obj_val;
				}
				//UPDATE gbest IF NECESSARY
				bool updateGbest = false;
				if (particles[i]->constr_viol > 0.0) {	//CURRENT PARTICLE IS INFEASIBLE
					if (particles[gbest]->pbest_constr_viol> 0.0) {	//gbest PARTICLE IS INFEASIBLE
						if (particles[i]->constr_viol < particles[gbest]->pbest_constr_viol) {
							updateGbest = true;
						}
						else if (particles[i]->constr_viol == particles[gbest]->pbest_constr_viol) {
							if (particles[gbest]->pbest_obj_val > particles[i]->obj_val) {
								updateGbest = true;
							}
						}
					}
				}
				else {	//CURRENT PARTICLE IS FEASIBLE
					if (particles[gbest]->pbest_constr_viol > 0.0) {	//gbest PARTICLE IS INFEASIBLE
						updateGbest = true;
					}
					else {	//gbest PARTICLE IS FEASIBLE
						if (particles[gbest]->pbest_obj_val > particles[i]->obj_val) {
							updateGbest = true;
						}
					}
				}
				if (updateGbest) {
					gbest = i;
				}
			}
		}
		//CREATE NEW COPY OF THE FITTEST PARTICLE AND MAKE NECESSARY CLEANUPS
		Solution* res = new Solution;
		Particle* source = particles[gbest];
		for (size_t i = 0; i < vsz; ++i) {
			res->solution[dec_var_traits[i].name] = source->pbest[i];
		}
		res->obj_val = source->pbest_obj_val;
		for (size_t i = 0; i < psz; ++i) {
			delete particles[i];
		}
		particles.clear();
		return res;
	}

	PSO_solver(
		const std::string& objective_func, const std::unordered_map<std::string, Var_Traits>& dec_vars,
		const std::unordered_map<std::string, numeric_type>& params, const std::vector<std::string>& constraints,
		numeric_type inertial_weight = numeric_type{ 0.5 }, numeric_type local_accelerator = numeric_type{ 1.4 },
		numeric_type global_accelerator = numeric_type{ 1.4 }, numeric_type cont_div_coef_not = numeric_type{ 0.1 },
		numeric_type cont_div_coef_min = numeric_type{ 1.0e-10 }, numeric_type disc_div_coef_not = numeric_type{ 0.7 },
		numeric_type lambda = numeric_type{ 0.2 }, size_t population_coef = 10, size_t max_iteration = 500,
		numeric_type cons_tol = numeric_type{ 0.001 }, numeric_type space_intervals = numeric_type{ 10.0 }) {
		size_t vsz = dec_vars.size(), disc_count = 0;
		if (vsz == 0) {
			throw "ERROR: dec_vars == 0";
		}
		if (inertial_weight <= 0.0 || inertial_weight > 1.0) {
			throw "ERROR: invalid inertial weight";
		}
		this->inertial_weight = inertial_weight;
		if (local_accelerator <= 0.0 || global_accelerator <= 0.0) {
			throw "ERROR: invalid accelerators";
		}
		this->local_accelerator = local_accelerator;
		this->global_accelerator = global_accelerator;
		if (cont_div_coef_not <= 0.0 || cont_div_coef_min <= 0.0 || disc_div_coef_not <= 0.0 || lambda < 0.0 || lambda > 1.0) {
			throw "ERROR: invalid diversity coefficient parameters";
		}
		this->cont_div_coef_not = cont_div_coef_not;
		this->cont_div_coef_min = cont_div_coef_min;
		this->disc_div_coef_not = disc_div_coef_not;
		this->lambda = lambda;
		if (population_coef == 0) {
			throw "ERROR: invalid population_coef";
		}
		this->population_coef = population_coef;
		numeric_type psz = static_cast<numeric_type>(vsz * population_coef);
		this->div_metric_const = std::pow(
			(psz + 1.0) / (psz * lambda + 1.0),
			numeric_type{ 1.0 } / static_cast<numeric_type>(vsz));
		this->l_param = (0.91*0.51) /
			(std::pow(psz, 0.21) * std::pow(static_cast<numeric_type>(vsz), 0.58));
		if (max_iteration <= 10) {
			throw "ERROR: too few iterations";
		}
		this->max_iteration = max_iteration;
		numeric_type cons_tolerance = std::abs(cons_tol);
		symbol_table.add_constants();
		for (auto itr = params.begin(), end = params.end(); itr != end; ++itr) {
			if (symbol_table.add_constant(itr->first, itr->second) == false) {
				throw "ERROR: invalid constant";
			}
		}
		if (space_intervals < 1.0) {
			throw "ERROR: invalid space interval number";
		}
		size_t i = 0;
		dec_var_values.resize(vsz);
		for (auto itr = dec_vars.begin(), end = dec_vars.end(); itr != end; ++itr, ++i) {
			if (symbol_table.add_variable(itr->first, dec_var_values[i]) == false) {
				throw "ERROR: invalid variable";
			}
			dec_var_traits.push_back(Var_Traits_Internal{});
			Var_Traits_Internal& var_traits = dec_var_traits.back();
			var_traits.lower_bound = itr->second.lower_bound;
			var_traits.upper_bound = itr->second.upper_bound;
			var_traits.type = itr->second.type;
			var_traits.name = itr->first;
			if (var_traits.type == Var_Traits::DISC) {
				var_traits.lower_bound = std::round(var_traits.lower_bound);
				var_traits.upper_bound = std::round(var_traits.upper_bound);
				var_traits.max_velocity = std::max(
					(var_traits.upper_bound - var_traits.lower_bound) / space_intervals, 
					numeric_type{ 0.5 });
				++disc_count;
			}
			else {
				var_traits.max_velocity =
					(var_traits.upper_bound - var_traits.lower_bound) / space_intervals;
			}
			if (var_traits.lower_bound >= var_traits.upper_bound) {
				throw "ERROR: invalid lower/upper bounds";
			}
		}
		this->cont_div_metrix_exp = numeric_type{ 1.0 } / static_cast<numeric_type>(vsz - disc_count);
		expression.register_symbol_table(symbol_table);
		for (size_t i = 0, csz = constraints.size(); i < csz; ++i) {
			std::string normalized;
			size_t sz = constraints[i].size();
			bool isValid = false;
			if (sz) --sz;
			for (size_t j = 0; j < sz; ++j) {
				if (constraints[i][j] == '>') {
					isValid = true;
					if (constraints[i][j + 1] == '=') {
						normalized.append("-(").append(constraints[i].substr(0, j)).
							append(")+(").append(constraints[i].substr(j + 2)).append(")");
					}
					else {
						normalized.append("-(").append(constraints[i].substr(0, j - 1)).
							append(")+(").append(constraints[i].substr(j + 1)).append(")+").
							append(std::to_string(cons_tolerance));
					}
					break;
				}
				else if (constraints[i][j] == '<') {
					isValid = true;
					if (constraints[i][j + 1] == '=') {
						normalized.append(constraints[i].substr(0, j)).
							append("-(").append(constraints[i].substr(j + 2)).append(")");
					}
					else {
						normalized.append(constraints[i].substr(0, j)).
							append("-(").append(constraints[i].substr(j + 1)).append(")").append("+").
							append(std::to_string(cons_tolerance));
					}
					break;
				}
				else if (constraints[i][j] == '=') {
					if (constraints[i][j + 1] == '=') {
						isValid = true;
						//MAKE ADDITIONAL push_back to constraint_expressions
						//h(x)-E <= 0, -h(x)-E <= 0
						normalized.append(constraints[i].substr(0, j)).
							append("-(").append(constraints[i].substr(j + 2)).append("+").
							append(std::to_string(cons_tolerance)).append(")");
						if (parser.compile(normalized, expression) == false) {
							throw parser.error().c_str();
						}
						constraint_expressions.push_back(expression);
						normalized.clear();
						normalized.append("-(").append(constraints[i].substr(0, j)).
							append(")+(").append(constraints[i].substr(j + 2)).append(")-").
							append(std::to_string(cons_tolerance));
						break;
					}
				}
			}
			if (!isValid) {
				throw "ERROR: invalid constraint";
			}
			if (parser.compile(normalized, expression) == false) {
				throw parser.error().c_str();
			}
			constraint_expressions.push_back(expression);
		}
		if (parser.compile(objective_func, expression) == false) {
			throw parser.error().c_str();
		}
	}
	~PSO_solver() {
		for (auto itr = particles.begin(), end = particles.end(); itr != end; ++itr) {
			delete *itr;
		}
	}

};