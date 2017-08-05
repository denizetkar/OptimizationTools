#include<vector>
#include<limits>
#include<utility>
#include<unordered_map>
#include<string>

template <typename numeric_type>
class GA_solver {
public:
	using t_type = double;
	using prob = double;
	using Solution = std::unordered_map<std::string, numeric_type>;
protected:
	using f_type = numeric_type;
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
	struct Gene {
		virtual ~Gene() {}
		virtual void init_rand() = 0;
		virtual void mutate(double) = 0;
	};
	struct Gene_Traits_Base {};
	struct Individual {
		std::vector<Gene*> genes;
		f_type fitness;
		numeric_type criteria;
		virtual bool operator<(const Individual& rhs) const {
			return criteria < rhs.criteria;
		}
		virtual bool operator>(const Individual& rhs) const {
			return criteria > rhs.criteria;
		}
		virtual ~Individual() {
			for (size_t i = 0, sz = genes.size(); i < sz; ++i) {
				if (genes[i]) {
					delete genes[i];
				}
			}
		}
	};
public:
	virtual Solution* solve() = 0;
protected:
	virtual f_type fitness(Individual const *) = 0;
	virtual void mutate(Individual*) = 0;
	virtual Individual* cross_over(Individual*, Individual*) = 0;
	virtual Individual* select_individual() = 0;

	prob p_m, p_c;
	Annealing* an; // for dynamic step size
	std::vector<Individual*> population;
};
