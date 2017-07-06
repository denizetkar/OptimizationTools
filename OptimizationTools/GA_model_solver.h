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
static std::uniform_int_distribution<> dis{ 0, 1 };

class GA_model_solver : public GA_solver {

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
public:
	using disc = long long;
	using cont = double;
	struct Gene_disc : Gene {
		disc val, lower, upper, step_size;
		Gene_disc(disc l, disc u, disc s_s) : lower{ l }, upper{ u }, step_size{ s_s } {
			if (l >= u) {
				throw "ERROR: l >= u";
			}
			if (s_s <= 0) {
				throw "ERROR: s_s <= 0";
			}
		}
		~Gene_disc() {}
		void init_rand() {
			val = lower + (disc)(unif(gen)*(upper - lower));
		}
		void mutate(double annealing) {
			disc delta = (disc)std::llround(2.0 * step_size * unif(gen) * annealing);
			if (dis(gen) == 0) {
				if (delta >= upper - val) {
					delta = upper - val;
				}
				val += delta;
			}
			else {
				if (delta >= val - lower) {
					delta = val - lower;
				}
				val -= delta;
			}
		}
	};
	struct Gene_cont : Gene {
		cont val, lower, upper, step_size;
		Gene_cont(cont l, cont u, cont s_s) : lower{ l }, upper{ u }, step_size{ s_s } {
			if (l >= u) {
				throw "ERROR: l >= u";
			}
			if (s_s <= 0.0) {
				throw "ERROR: s_s <= 0";
			}
		}
		~Gene_cont() {}
		void init_rand() {
			val = lower + unif(gen)*(upper - lower);
		}
		void mutate(double annealing) {
			cont delta = (cont)(2.0 * step_size * unif(gen) * annealing);
			if (dis(gen) == 0) {
				if (delta >= upper - val) {
					delta = upper - val;
				}
				val += delta;
			}
			else {
				if (delta >= val - lower) {
					delta = val - lower;
				}
				val -= delta;
			}
		}
	};
	struct Gene_Traits : Gene_Traits_Base {
		enum { DISC, CONT } type;
		union Lower {
			disc discrete;
			cont continuous;
		} lower_bound;
		union Upper {
			disc discrete;
			cont continuous;
		} upper_bound;
		union Stpsz {
			disc discrete;
			cont continuous;
		} step_size;
	};
private:
	using Ind_ext_t = std::pair<Individual_ext*, size_t>;
	struct Comp_Ind_Ext_Ptr_Fitness {
		bool operator()(const Ind_ext_t& o1, const Ind_ext_t& o2) {
			return o1.first->fitness > o2.first->fitness;
		}
	};
	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;
	typedef exprtk::parser_error::type error_t;

	symbol_table_t symbol_table;
	expression_t expression;
	parser_t parser;
	std::vector<double> dec_var_values;
	std::vector<std::string> dec_var_names;
	std::vector<Gene_Traits> dec_var_traits;
	std::unordered_set<size_t> elite_indexes;	//MAKE SURE IT IS USED!!!!!
	prob p_s;
	size_t population_coef, elit_num, generation_max;

	bool check_prob(prob p) {
		return (p >= 0.0 && p <= 1.0);
	}
	Individual* get_individual(bool isRandom = false) {
		Individual* res = new Individual;
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
		double x;
		for (size_t i = 0, gsz = dec_var_traits.size(); i < gsz; ++i) {
			switch (dec_var_traits[i].type) {
			case Gene_Traits::DISC:
				x = (double)((Gene_disc*)(ind->genes[i]))->val;
				break;
			case Gene_Traits::CONT:
				x = (double)((Gene_cont*)(ind->genes[i]))->val;
				break;
			default:
				x = std::numeric_limits<double>::max();
				break;
			}
			dec_var_values[i] = x;
		}
		return (f_type)expression.value();
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
		Individual* res = new Individual;
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
		for (; i < sz && p <= population[i]->cumul_s; ++i);
		return population[i]->ind;
	}
	void process_population() {
		std::vector<double> gene_avgs;
		size_t gsz = dec_var_values.size();
		size_t psz = population.size();
		gene_avgs.resize(gsz, 0.0);
		f_type fitness_max = -std::numeric_limits<f_type>::max();
		double diversity_max = -std::numeric_limits<double>::max();
		double x, y;
		prob cum, cur;

		for (size_t i = 0; i < psz; ++i) {
			for (size_t j = 0; j < gsz; ++j) {
				switch (dec_var_traits[j].type) {
				case Gene_Traits::DISC:
					x = (double)((Gene_disc*)(population[i]->ind->genes[j]))->val;
					break;
				case Gene_Traits::CONT:
					x = (double)((Gene_cont*)(population[i]->ind->genes[j]))->val;
					break;
				default:
					x = std::numeric_limits<double>::max();
					break;
				}
				dec_var_values[j] = x;
				gene_avgs[j] += (x - gene_avgs[j]) / (i + 1);
			}
			if (fitness_max < (population[i]->fitness = (f_type)expression.value())) {
				fitness_max = population[i]->fitness;
			}
		}
		for (size_t i = 0; i < psz; ++i) {
			population[i]->diversity = 0.0;
			for (size_t j = 0; j < gsz; ++j) {
				switch (dec_var_traits[j].type) {
				case Gene_Traits::DISC:
					x = (double)((Gene_disc*)(population[i]->ind->genes[j]))->val - gene_avgs[j];
					break;
				case Gene_Traits::CONT:
					x = (double)((Gene_cont*)(population[i]->ind->genes[j]))->val - gene_avgs[j];
					break;
				default:
					x = std::numeric_limits<double>::max();
					break;
				}
				population[i]->diversity += x*x;
			}
			if (diversity_max < population[i]->diversity) {
				diversity_max = population[i]->diversity;
			}
		}
		for (size_t i = 0; i < psz; ++i) {
			x = fitness_max - population[i]->fitness;
			y = diversity_max - population[i]->diversity;
			population[i]->criteria = x*x + y*y;
		}
		std::sort(population.begin(), population.end(), Individual_ext::Comp_Ptr_Smaller{});
		cum = 0.0;
		cur = p_s;
		--psz;
		for (size_t i = 0; i < psz; ++i) {
			cum += cur;
			population[i]->cumul_s = cum;
			cur *= (1 - p_s);
		}
		population[psz]->cumul_s = 1.0;
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
	Individual* solve() {
		f_type fitness_max = -std::numeric_limits<f_type>::max();
		size_t fit_max_index = 0;
		size_t gsz = dec_var_values.size();
		size_t psz = population_coef * gsz;
		population.resize(psz);
		size_t npsz = psz - elit_num;
		std::vector<Individual_ext*> newPopulation;
		newPopulation.resize(npsz);
		an->init();
		//INITIALIZE THE POPULATION RANDOMLY
		for (size_t i = 0; i < psz; ++i) {
			Individual_ext* ind_ext = new Individual_ext;
			ind_ext->ind = get_individual(true);
			population[i] = ind_ext;
		}
		for (size_t generation = 0; generation < generation_max; ++generation) {
			//PROCESS POPULATION TO FIND CUMULATIVE PROBABILITY OF SELECTION FOR EACH INDIVIDUAL
			process_population();
			//FIND 'elit_num' NUMBER OF ELITES WITH HIGHEST FITNESS VALUES
			//(fitness values are already calculated in 'process_population' function)
			binary_heap<Ind_ext_t, std::vector<Ind_ext_t>, Comp_Ind_Ext_Ptr_Fitness> q;
			for (size_t i = 0; i < psz; ++i) {
				if (q.size() < elit_num) {
					q.push(Ind_ext_t{ population[i],i });
				}
				else if (population[q.top().second]->fitness < population[i]->fitness) {
					q.pop(); q.push(Ind_ext_t{ population[i],i });
				}
			}
			std::vector<Ind_ext_t>& container = q.get_container();
			elite_indexes.clear();
			for (size_t i = 0; i < elit_num; ++i) {
				elite_indexes.insert(container[i].second);
			}
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
				newPopulation[i] = new Individual_ext{ ind };
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
					double x;
					switch (dec_var_traits[j].type) {
					case Gene_Traits::DISC:
						x = (double)((Gene_disc*)(population[i]->ind->genes[j]))->val;
						break;
					case Gene_Traits::CONT:
						x = (double)((Gene_cont*)(population[i]->ind->genes[j]))->val;
						break;
					default:
						x = std::numeric_limits<double>::max();
						break;
					}
					dec_var_values[j] = x;
				}
				if (fitness_max < (population[i]->fitness = (f_type)expression.value())) {
					fitness_max = population[i]->fitness;
					fit_max_index = i;
				}
			}
		}
		//CREATE NEW COPY OF THE FITTEST INDIVIDUAL AND MAKE NECESSARY CLEANUPS
		Individual* res = get_individual();
		set_individual(res, population[fit_max_index]->ind);
		for (size_t i = 0; i < psz; ++i) {
			delete population[i];
		}
		elite_indexes.clear();
		population.clear();
		return res;
	}
	
	//population coefficient has to be >= 1
	GA_model_solver(
		const std::string& objective_func, const std::unordered_map<std::string, Gene_Traits>& dec_vars, 
		const std::unordered_map<std::string, cont>& params, size_t pop_coef, double elit_ratio, size_t generation = 1000, 
		double cool_rate = 0.99, double conv_rate = 0.1, prob mut = 0.1, prob cro = 0.5, prob select = 0.2,
		t_type temp_min = 0.0, t_type temp_max = 100000.0) {
		size_t vsz = dec_vars.size();
		if (vsz == 0) {
			throw "ERROR: dec_vars == 0";
		}
		if (generation <= 10) {
			throw "ERROR: too few generations";
		}
		generation_max = generation;
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
		elit_num = (size_t)std::llround(elit_ratio * pop_coef * vsz);
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

	void print_individual(const Individual * ind, std::ostream& out = std::cout) const {
		for (size_t i = 0, gsz = dec_var_traits.size(); i < gsz; ++i) {
			switch (dec_var_traits[i].type) {
			case Gene_Traits::DISC:
				out << dec_var_names[i] << ": " << ((Gene_disc*)(ind->genes[i]))->val << " ";
				break;
			case Gene_Traits::CONT:
				out << dec_var_names[i] << ": " << ((Gene_cont*)(ind->genes[i]))->val << " ";
				break;
			default: break;
			}
		}
	}
};