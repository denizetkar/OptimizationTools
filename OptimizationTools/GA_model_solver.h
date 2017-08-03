#include "GA_solver.h"
#include "exprtk.hpp"
#include "binary_heap.h"
#include<vector>
#include<limits>
#include<cmath>
#include<random>
#include<string>
#include<sstream>
#include<unordered_map>
#include<unordered_set>
#include<algorithm>
#include<iostream>

static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<double> unif{ 0.0, 1.0 };
static std::normal_distribution<double> norm{};

//TRIES TO MAXIMIZE THE OBJECTIVE FUNCTION AND SATISFY ALL CONSTRAINTS
template <typename numeric_type = double, typename discrete_type = long long>
class GA_model_solver : public GA_solver<numeric_type> {

	class Anneal_impl : public Annealing {
		t_type base;
	public:
		Anneal_impl(t_type min, t_type max, double c_rate, double convergence_rate) : Annealing{ min, max, c_rate },
			base{ (max - min) * (1.0 - c_rate) * convergence_rate } {}
		void init() {
			temperature = temp_max;
		}
		double get_annealing() {
			return ((double)(temperature - temp_min)) / (temp_max - temp_min);
		}
		void cool_down() {
			temperature = temp_min + base + (t_type)(cool_rate * (temperature - temp_min));
		}
	};
	struct Gene_disc : Gene {
		discrete_type val, lower, upper, step_size;
		Gene_disc(discrete_type l, discrete_type u, discrete_type s_s) : lower{ l }, upper{ u }, step_size{ s_s } {
			if (l >= u) {
				throw "ERROR: l >= u";
			}
			if (s_s <= 0) {
				throw "ERROR: s_s <= 0";
			}
		}
		~Gene_disc() {}
		void init_rand() {
			val = lower + static_cast<discrete_type>(static_cast<numeric_type>(unif(gen))*static_cast<numeric_type>(upper - lower));
		}
		void mutate(double annealing) {
			discrete_type delta = static_cast<discrete_type>(std::round(static_cast<numeric_type>(step_size) * static_cast<numeric_type>(norm(gen) * annealing)));
			//discrete_type delta = static_cast<discrete_type>(std::llround(2.0 * step_size * (2.0 * unif(gen) - 1.0) * annealing));
			if (delta > upper - val) {
				delta = upper - val;
			}
			else if (delta < lower - val) {
				delta = lower - val;
			}
			val += delta;
		}
	};
	struct Gene_cont : Gene {
		numeric_type val, lower, upper, step_size;
		Gene_cont(numeric_type l, numeric_type u, numeric_type s_s) : lower{ l }, upper{ u }, step_size{ s_s } {
			if (l >= u) {
				throw "ERROR: l >= u";
			}
			if (s_s <= 0.0) {
				throw "ERROR: s_s <= 0";
			}
		}
		~Gene_cont() {}
		void init_rand() {
			val = lower + (upper - lower) * static_cast<numeric_type>(unif(gen));
		}
		void mutate(double annealing) {
			numeric_type delta = step_size * static_cast<numeric_type>(norm(gen) * annealing);
			//numeric_type delta = static_cast<numeric_type>(2.0 * step_size * (2.0 * unif(gen) - 1.0) * annealing);
			if (delta > upper - val) {
				delta = upper - val;
			}
			else if (delta < lower - val) {
				delta = lower - val;
			}
			val += delta;
		}
	};
	struct Individual_ext : Individual {
		prob cumul_s;
		numeric_type diversity;
		numeric_type cons_viol;
		size_t num_of_cons_viol;
		Individual_ext() {}
		~Individual_ext() {}
	};
public:
	struct Gene_Traits : Gene_Traits_Base {
		enum { DISC, CONT } type;
		union Lower {
			discrete_type discrete;
			numeric_type continuous;
		} lower_bound;
		union Upper {
			discrete_type discrete;
			numeric_type continuous;
		} upper_bound;
		union Stpsz {
			discrete_type discrete;
			numeric_type continuous;
		} step_size;
	};
	using cont_type = numeric_type;
	using disc_type = discrete_type;
private:
	using Ind_ext_t = std::pair<Individual*, size_t>;
	struct Comp_Ind_Ext_Ptr_Selection {
		bool operator()(const Individual* o1, const Individual* o2) {
			if (((const Individual_ext*)o1)->num_of_cons_viol == 0) {
				if (((const Individual_ext*)o2)->num_of_cons_viol == 0) {
					return o1->criteria < o2->criteria;
				}
				else {
					return true;
				}
			}
			else {
				if (((const Individual_ext*)o2)->num_of_cons_viol == 0) {
					return false;
				}
				else {
					if (((const Individual_ext*)o1)->num_of_cons_viol == ((const Individual_ext*)o2)->num_of_cons_viol) {
						return ((const Individual_ext*)o1)->cons_viol < ((const Individual_ext*)o2)->cons_viol;
					}
					else {
						return ((const Individual_ext*)o1)->num_of_cons_viol < ((const Individual_ext*)o2)->num_of_cons_viol;
					}
				}
			}
		}
	};
	struct Comp_Ind_Ext_t_Ptr_Elite_Selection {
		bool operator()(const Ind_ext_t& o1, const Ind_ext_t& o2) {
			if (((const Individual_ext*)o1.first)->num_of_cons_viol == ((const Individual_ext*)o2.first)->num_of_cons_viol) {
				return o1.first->fitness < o2.first->fitness;
			}
			else if (((const Individual_ext*)o1.first)->num_of_cons_viol < ((const Individual_ext*)o2.first)->num_of_cons_viol) {
				return false;
			}
			else {
				return true;
			}
		}
	};
	struct Comp_Ind_Ext_Ptr_Elite_Selection {
		bool operator()(const Individual* o1, const Individual* o2) {
			if (((const Individual_ext*)o1)->num_of_cons_viol == ((const Individual_ext*)o2)->num_of_cons_viol) {
				return o1->fitness < o2->fitness;
			}
			else if (((const Individual_ext*)o1)->num_of_cons_viol < ((const Individual_ext*)o2)->num_of_cons_viol) {
				return false;
			}
			else {
				return true;
			}
		}
	};
	typedef exprtk::symbol_table<numeric_type> symbol_table_t;
	typedef exprtk::expression<numeric_type>     expression_t;
	typedef exprtk::parser<numeric_type>             parser_t;
	typedef exprtk::parser_error::type error_t;

	symbol_table_t symbol_table;
	std::vector<expression_t> constraint_expressions;
	expression_t expression;
	parser_t parser;
	std::vector<numeric_type> dec_var_values;
	std::vector<std::string> dec_var_names;
	std::vector<Gene_Traits> dec_var_traits;
	prob p_s;
	size_t population_coef, elit_num, generation_max;
	numeric_type prev_fit, cur_fit;

	bool check_prob(prob p) {
		return (p >= 0.0 && p <= 1.0);
	}
	Individual* get_individual(bool isRandom = false) {
		Individual* res = new Individual_ext;
		for (auto itr = dec_var_traits.begin(), end = dec_var_traits.end(); itr != end; ++itr) {
			Gene* g;
			switch (itr->type) {
			case Gene_Traits::DISC:
				g = (Gene*) new Gene_disc{ itr->lower_bound.discrete, itr->upper_bound.discrete, itr->step_size.discrete };
				break;
			case Gene_Traits::CONT:
				g = (Gene*) new Gene_cont{ itr->lower_bound.continuous, itr->upper_bound.continuous, itr->step_size.continuous };
				break;
			default:
				throw "ERROR: unrecognized gene type";
				break;
			}
			if (isRandom) {
				g->init_rand();
			}
			res->genes.push_back(g);
		}
		return res;
	}
	f_type fitness(Individual const * ind) {
		numeric_type x;
		for (size_t i = 0, gsz = dec_var_traits.size(); i < gsz; ++i) {
			switch (dec_var_traits[i].type) {
			case Gene_Traits::DISC:
				x = static_cast<numeric_type>(((Gene_disc*)(ind->genes[i]))->val);
				break;
			case Gene_Traits::CONT:
				x = static_cast<numeric_type>(((Gene_cont*)(ind->genes[i]))->val);
				break;
			default:
				x = std::numeric_limits<numeric_type>::max();
				break;
			}
			dec_var_values[i] = x;
		}
		return static_cast<f_type>(expression.value());
	}
	void mutate(Individual* ind) {
		double ann = an->get_annealing();
		for (size_t i = 0, gsz = ind->genes.size(); i < gsz; ++i) {
			if (unif(gen) <= p_m) {
				ind->genes[i]->mutate(ann);
			}
		}
	}
	Individual* cross_over(Individual* ind1, Individual* ind2) {
		Individual* res = new Individual_ext;
		for (size_t i = 0, gsz = dec_var_traits.size(); i < gsz; ++i) {
			Gene* g;
			switch (dec_var_traits[i].type) {
			case Gene_Traits::DISC:
				g = (Gene*) new Gene_disc{ dec_var_traits[i].lower_bound.discrete, dec_var_traits[i].upper_bound.discrete, dec_var_traits[i].step_size.discrete };
				if (unif(gen) <= p_c) {
					((Gene_disc*)g)->val = ((Gene_disc*)ind1->genes[i])->val;
				}
				else {
					((Gene_disc*)g)->val = ((Gene_disc*)ind2->genes[i])->val;
				}
				break;
			case Gene_Traits::CONT:
				g = (Gene*) new Gene_cont{ dec_var_traits[i].lower_bound.continuous, dec_var_traits[i].upper_bound.continuous, dec_var_traits[i].step_size.continuous };
				if (unif(gen) <= p_c) {
					((Gene_cont*)g)->val = ((Gene_cont*)ind1->genes[i])->val;
				}
				else {
					((Gene_cont*)g)->val = ((Gene_cont*)ind2->genes[i])->val;
				}
				break;
			default:
				g = nullptr;
				break;
			}
			res->genes.push_back(g);
		}
		return res;
	}
	Individual* select_individual() {
		//Assumes that population is sorted
		size_t i = 0, sz = population.size();
		if (sz == 0) { return nullptr; }
		--sz;
		prob p = unif(gen);
		for (; i < sz && p <= ((Individual_ext*)population[i])->cumul_s; ++i);
		return population[i];
	}
	void process_population() {	//CALCULATES FIELDS OF INDIVIDUALS IN THE POPULATION
		std::vector<numeric_type> gene_avgs;
		size_t gsz = dec_var_values.size();
		size_t psz = population.size();
		gene_avgs.resize(gsz, 0.0);
		f_type fitness_max = -std::numeric_limits<f_type>::max();
		numeric_type diversity_max = -std::numeric_limits<numeric_type>::max();
		numeric_type x, y;
		prob cum, cur;

		for (size_t i = 0; i < psz; ++i) {
			for (size_t j = 0; j < gsz; ++j) {
				switch (dec_var_traits[j].type) {
				case Gene_Traits::DISC:
					x = static_cast<numeric_type>(((Gene_disc*)(population[i]->genes[j]))->val);
					break;
				case Gene_Traits::CONT:
					x = static_cast<numeric_type>(((Gene_cont*)(population[i]->genes[j]))->val);
					break;
				default:
					x = std::numeric_limits<numeric_type>::max();
					break;
				}
				dec_var_values[j] = x;
				gene_avgs[j] += (x - gene_avgs[j]) / (i + 1);
			}
			if (std::isnan(population[i]->fitness = static_cast<f_type>(expression.value()))) {
				population[i]->fitness = -std::numeric_limits<numeric_type>::max();
			}
			if (fitness_max < population[i]->fitness) {
				fitness_max = population[i]->fitness;
			}
			((Individual_ext*)population[i])->cons_viol = 0.0;
			((Individual_ext*)population[i])->num_of_cons_viol = 0;
			for (size_t j = 0, consz = constraint_expressions.size(); j < consz; ++j) {
				numeric_type cons_val = constraint_expressions[j].value();
				if (std::isnan(cons_val)) {
					cons_val = -std::numeric_limits<numeric_type>::max();
				}
				if (cons_val > 0.0) {	//CONSTRAINT VIOLATED!!!
					++((Individual_ext*)population[i])->num_of_cons_viol;
					((Individual_ext*)population[i])->cons_viol += cons_val;
				}
			}
		}
		cur_fit = fitness_max;
		for (size_t i = 0; i < psz; ++i) {
			((Individual_ext*)population[i])->diversity = 0.0;
			for (size_t j = 0; j < gsz; ++j) {
				switch (dec_var_traits[j].type) {
				case Gene_Traits::DISC:
					x = static_cast<numeric_type>(((Gene_disc*)(population[i]->genes[j]))->val - gene_avgs[j]);
					break;
				case Gene_Traits::CONT:
					x = static_cast<numeric_type>(((Gene_cont*)(population[i]->genes[j]))->val - gene_avgs[j]);
					break;
				default:
					x = std::numeric_limits<numeric_type>::max();
					break;
				}
				((Individual_ext*)population[i])->diversity += x*x;
			}
			if (diversity_max < ((Individual_ext*)population[i])->diversity) {
				diversity_max = ((Individual_ext*)population[i])->diversity;
			}
		}
		for (size_t i = 0; i < psz; ++i) {
			x = fitness_max - population[i]->fitness;
			y = diversity_max - ((Individual_ext*)population[i])->diversity;
			population[i]->criteria = x*x + y*y;
		}
		std::sort(population.begin(), population.end(), Comp_Ind_Ext_Ptr_Selection{});
		cum = 0.0;
		cur = p_s;
		--psz;
		for (size_t i = 0; i < psz; ++i) {
			cum += cur;
			((Individual_ext*)population[i])->cumul_s = cum;
			cur *= (1 - p_s);
		}
		((Individual_ext*)population[psz])->cumul_s = 1.0;
	}
	void set_individual(Individual* dest, Individual const * source) {
		for (size_t i = 0, gsz = dec_var_traits.size(); i < gsz; ++i) {
			switch (dec_var_traits[i].type) {
			case Gene_Traits::DISC:
				((Gene_disc*)dest->genes[i])->val = ((Gene_disc*)source->genes[i])->val;
				break;
			case Gene_Traits::CONT:
				((Gene_cont*)dest->genes[i])->val = ((Gene_cont*)source->genes[i])->val;
				break;
			default: break;
			}
		}
	}
public:
	Solution* solve() {
		std::unordered_set<size_t> elite_indexes;
		f_type fitness_max = -std::numeric_limits<f_type>::max();
		size_t fit_max_index = 0;
		Comp_Ind_Ext_Ptr_Elite_Selection elit_comp;
		size_t gsz = dec_var_values.size();
		size_t psz = population_coef * gsz;
		population.resize(psz);
		size_t npsz = psz - elit_num;
		std::vector<Individual*> newPopulation;
		newPopulation.resize(npsz);
		an->init();
		//INITIALIZE THE POPULATION RANDOMLY
		for (size_t i = 0; i < psz; ++i) {
			population[i] = get_individual(true);
		}
		for (size_t generation = 0; generation < generation_max; ++generation) {
			//PROCESS POPULATION TO FIND CUMULATIVE PROBABILITY OF SELECTION FOR EACH INDIVIDUAL
			process_population();

			//FIND 'elit_num' NUMBER OF ELITES WITH HIGHEST FITNESS VALUES
			//(fitness values are already calculated in 'process_population' function)
			binary_heap<Ind_ext_t, std::vector<Ind_ext_t>, Comp_Ind_Ext_t_Ptr_Elite_Selection> q;
			for (size_t i = 0; i < psz; ++i) {
				if (q.size() < elit_num) {
					q.push(Ind_ext_t{ population[i],i });
				}
				else if (elit_comp(population[q.top().second], population[i])) {
					q.pop(); q.push(Ind_ext_t{ population[i],i });
				}
			}
			std::vector<Ind_ext_t>& container = q.get_container();
			elite_indexes.clear();
			for (size_t i = 0; i < elit_num; ++i) {
				elite_indexes.insert(container[i].second);
			}

			if (cur_fit == prev_fit) {
				auto elite_not_found = elite_indexes.end();
				for (size_t i = psz / 2; i < psz; ++i) {
					if (elite_indexes.find(i) == elite_not_found) {
						delete population[i];
						population[i] = get_individual(true);
					}
				}
				//PROCESS POPULATION TO FIND CUMULATIVE PROBABILITY OF SELECTION FOR EACH INDIVIDUAL
				process_population();

				//FIND 'elit_num' NUMBER OF ELITES WITH HIGHEST FITNESS VALUES
				//(fitness values are already calculated in 'process_population' function)
				q.clear();
				for (size_t i = 0; i < psz; ++i) {
					if (q.size() < elit_num) {
						q.push(Ind_ext_t{ population[i],i });
					}
					else if (elit_comp(population[q.top().second], population[i])) {
						q.pop(); q.push(Ind_ext_t{ population[i],i });
					}
				}
				std::vector<Ind_ext_t>& container = q.get_container();
				elite_indexes.clear();
				for (size_t i = 0; i < elit_num; ++i) {
					elite_indexes.insert(container[i].second);
				}
			}
			prev_fit = cur_fit;
			/*
			//FIRST 'elit_num' INDIVIDUALS ARE MARKED AS ELITE
			elite_indexes.clear();
			for (size_t i = 0; i < elit_num; ++i) {
			elite_indexes.insert(i);
			}
			*/
			//CREATE NEW POPULATION TO REPLACE NON-ELITE INDIVIDUALS
			for (size_t i = 0; i < npsz; ++i) {
				Individual* ind1 = select_individual();
				Individual* ind2 = select_individual();
				Individual* ind;

				if (ind1 == ind2) {
					ind = get_individual();
					set_individual(ind, ind1);
				}
				else {
					ind = cross_over(ind1, ind2);
				}
				//MUTATING POSSIBLE ELITE !!!!!!!!!!!!!!!!!!!!!!!!!!
				mutate(ind);
				newPopulation[i] = ind;
			}
			//REPLACE NON-ELITES WITH NEW POPULATION
			auto elite_not_found = elite_indexes.end();
			for (size_t i = 0, index = 0; i < psz; ++i) {
				if (elite_indexes.find(i) == elite_not_found) {
					delete population[i];
					population[i] = newPopulation[index++];
				}
			}
			//ANNEALING COOL DOWN REDUCES MEAN STEP SIZE
			an->cool_down();
		}
		//FIND INDIVIDUAL WITH HIGHEST FITNESS
		auto elite_not_found = elite_indexes.end();
		for (size_t i = 0; i < psz; ++i) {
			if (elite_indexes.find(i) != elite_not_found) {
				if (fitness_max < population[i]->fitness) {
					fitness_max = population[i]->fitness;
					fit_max_index = i;
				}
			}
			else {
				for (size_t j = 0; j < gsz; ++j) {
					numeric_type x;
					switch (dec_var_traits[j].type) {
					case Gene_Traits::DISC:
						x = static_cast<numeric_type>(((Gene_disc*)(population[i]->genes[j]))->val);
						break;
					case Gene_Traits::CONT:
						x = static_cast<numeric_type>(((Gene_cont*)(population[i]->genes[j]))->val);
						break;
					default:
						x = std::numeric_limits<numeric_type>::max();
						break;
					}
					dec_var_values[j] = x;
				}
				if (std::isnan(population[i]->fitness = static_cast<f_type>(expression.value()))) {
					population[i]->fitness = -std::numeric_limits<numeric_type>::max();
				}
				if (fitness_max < population[i]->fitness) {
					fitness_max = population[i]->fitness;
					fit_max_index = i;
				}
			}
		}
		//CREATE NEW COPY OF THE FITTEST INDIVIDUAL AND MAKE NECESSARY CLEANUPS
		Solution* res = new Solution;
		Individual* source = population[fit_max_index];
		for (size_t i = 0; i < gsz; ++i) {
			switch (dec_var_traits[i].type) {
			case Gene_Traits::DISC:
				(*res)[dec_var_names[i]] = static_cast<numeric_type>(((Gene_disc*)source->genes[i])->val);
				break;
			case Gene_Traits::CONT:
				(*res)[dec_var_names[i]] = static_cast<numeric_type>(((Gene_cont*)source->genes[i])->val);
				break;
			default: break;
			}
		}
		for (size_t i = 0; i < psz; ++i) {
			delete population[i];
		}
		population.clear();
		return res;
	}

	//population coefficient has to be >= 1
	GA_model_solver(
		const std::string& objective_func, const std::unordered_map<std::string, Gene_Traits>& dec_vars,
		const std::unordered_map<std::string, numeric_type>& params, const std::vector<std::string>& constraints,
		size_t pop_coef = 10, double elit_ratio = 0.1, size_t generation = 1000, double cool_rate = 0.99, double conv_rate = 0.1,
		numeric_type cons_tol = numeric_type{ 0.001 }, prob mut = 1.0, prob cro = 0.5, prob select = 0.2,
		t_type temp_min = 0.0, t_type temp_max = 100000.0) {
		size_t vsz = dec_vars.size();
		if (vsz == 0) {
			throw "ERROR: dec_vars == 0";
		}
		prev_fit = -std::numeric_limits<numeric_type>::max();
		if (generation <= 10) {
			throw "ERROR: too few generations";
		}
		generation_max = generation;
		numeric_type cons_tolerance = std::abs(cons_tol);
		if (!check_prob(p_m) || !check_prob(cro) || !check_prob(p_s)) {
			throw "ERROR: invalid probability";
		}
		if (vsz == 1) {
			p_m = 1.0;
		}
		else {
			p_m = mut;
		}
		p_c = cro;
		p_s = select;
		if (pop_coef == 0) {
			throw "ERROR: invalid pop_coef";
		}
		population_coef = pop_coef;
		if (elit_ratio <= 0 || elit_ratio >= 1.0) {
			throw "ERROR: invalid elit_n";
		}
		elit_num = static_cast<size_t>(std::round(elit_ratio * pop_coef * vsz));
		if (elit_num == pop_coef * vsz) {
			--elit_num;
		}
		else if (elit_num == 0) {
			++elit_num;
		}
		an = nullptr;
		an = (Annealing*) new Anneal_impl{ temp_min, temp_max, cool_rate, conv_rate };
		symbol_table.add_constants();
		for (auto itr = params.begin(), end = params.end(); itr != end; ++itr) {
			if (symbol_table.add_constant(itr->first, itr->second) == false) {
				throw "ERROR: invalid constant";
			}
		}
		size_t i = 0;
		dec_var_values.resize(vsz);
		for (auto itr = dec_vars.begin(), end = dec_vars.end(); itr != end; ++itr) {
			if (symbol_table.add_variable(itr->first, dec_var_values[i++]) == false) {
				throw "ERROR: invalid variable";
			}
			dec_var_names.push_back(itr->first);
			dec_var_traits.push_back(itr->second);
		}
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
	~GA_model_solver() {
		if (an) { delete an; }
		for (auto itr = population.begin(), end = population.end(); itr != end; ++itr) {
			delete *itr;
		}
	}

	friend std::ostream& operator<<(std::ostream& out, const Solution* soln);
};

std::ostream& operator<<(std::ostream& out, const GA_model_solver<>::Solution* soln) {
	if (soln) {
		for (auto itr = soln->begin(), end = soln->end(); itr != end; ++itr) {
			out << itr->first << ": " << itr->second << " ";
		}
	}
	return out;
}