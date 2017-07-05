#include<vector>
#include<limits>
#include<utility>

class GA_solver {
public:
	using t_type = double;
	using prob = double;
protected:
	using f_type = double;
	class Annealing {
	protected:
		t_type temperature;
		t_type temp_min, temp_max;
		double cool_rate;
	public:
		Annealing(t_type min, t_type max, double c_rate) : temperature{ max }, temp_min{ min }, temp_max{ max }, cool_rate{ c_rate } {
			if (min >= max) {
				throw "ERROR: min >= max";
			}
			if (c_rate < 0.0 || c_rate > 1.0) {
				throw "ERROR: c_rate is invalid";
			}
		}
		virtual ~Annealing() {}
		virtual void init() = 0;
		virtual double get_annealing() = 0;
		virtual void cool_down() = 0;
	};
public:
	struct Gene {
		virtual ~Gene() {}
		virtual void init_rand() = 0;
		virtual void mutate(double) = 0;
	};
	struct Gene_Traits_Base {};
	struct Individual {
		std::vector<Gene*> genes;
		virtual ~Individual() {
			for (size_t i = 0, sz = genes.size(); i < sz; ++i) {
				if (genes[i]) {
					delete genes[i];
				}
			}
		}
	};
	virtual Individual* solve() = 0;
protected:
	struct Individual_ext {
		Individual* ind;
		prob cumul_s;
		double criteria;
		double diversity;
		f_type fitness;
		Individual_ext(Individual* _ind = nullptr) : ind{ _ind } {}
		virtual ~Individual_ext() {
			if (ind) {
				delete ind;
			}
		}
		bool operator<(const Individual_ext& rhs) const {
			return criteria < rhs.criteria;
		}
		bool operator>(const Individual_ext& rhs) const {
			return criteria > rhs.criteria;
		}
		struct Comp_Ptr_Smaller {
			bool operator()(const Individual_ext * i1, const Individual_ext * i2) {
				return *i1 < *i2;
			}
		};
		struct Comp_Ptr_Bigger {
			bool operator()(const Individual_ext * i1, const Individual_ext * i2) {
				return *i1 > *i2;
			}
		};
	};

	virtual f_type fitness(Individual const *) = 0;
	virtual void mutate(Individual*) = 0;
	virtual Individual* cross_over(Individual*, Individual*) = 0;
	virtual Individual* select_individual() = 0;

	prob p_m, p_c;
	Annealing* an; // for dynamic step size
	std::vector<Individual_ext*> population;
};