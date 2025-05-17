#include "ogdf_instance.hpp"
#include "ogdf_solver.hpp"
#include "ogdf_util.hpp"

uint64_t hash_clauses(const std::vector<std::vector<int>>& clauses) {
	uint64_t hash = FNV1a_64_SEED;
	OGDF_ASSERT(clauses.size() > 0);
	FNV1a_64_update(hash, clauses.size());
	for (auto& clause : clauses) {
		OGDF_ASSERT(clause.size() > 0);
		FNV1a_64_update(hash, clause.size());
		for (auto v : clause) {
			FNV1a_64_update(hash, v);
		}
	}
	return hash;
}

void dump_sat(const std::string& file, const std::vector<std::vector<int>>& clauses) {
	std::ofstream f(file);
	for (auto& clause : clauses) {
		bool first = true;
		for (auto x : clause) {
			if (first) {
				first = false;
			} else {
				f << " ";
			}
			f << x;
		}
		f << "\n";
	}
}

bool is_same_sat(const std::string& file, const std::vector<std::vector<int>>& clauses) {
	std::ifstream f(file);
	int a;
	for (auto& clause : clauses) {
		for (auto x : clause) {
			if (!(f >> a) || a != x) {
				return false;
			}
		}
	}
	if (f >> a) {
		return false;
	} else {
		return true;
	}
}

std::string get_filename(uint64_t hash, const std::string& ext, const std::string& dir) {
	std::stringstream fnstr;
	fnstr << dir;
	fnstr << std::setfill('0') << std::setw(sizeof(uint64_t) * 2) << std::hex << hash;
	fnstr << ext;
	return fnstr.str();
}

bool try_load_solution(Instance& I, std::vector<std::vector<int>>& hclauses, std::string& filename) {
	std::sort(hclauses.begin(), hclauses.end());
	uint64_t hash = hash_clauses(hclauses);
	filename = get_filename(hash, ".sol");
	std::string filename_sat = get_filename(hash, ".sat");

	int before = I.DS.size();
	if (std::filesystem::exists(filename)) {
		log << "Found cached solution " << filename << std::endl;
		bool can_load = true;
		if (!is_same_sat(filename_sat, hclauses)) {
			log << "Hash collision! Solution " << filename << " was generated for different input!"
				<< std::endl;
			can_load = false;
			for (int i = 0; i < 100; ++i) {
				filename = get_filename(hash, ".sol.col" + std::to_string(i));
				filename_sat = get_filename(hash, ".sat.col" + std::to_string(i));
				if (!std::filesystem::exists(filename)) {
					break;
				} else if (is_same_sat(filename_sat, hclauses)) {
					log << "Will load solution from " << filename
						<< " that was generated for the same input." << std::endl;
					can_load = true;
					break;
				}
			}
		}
		if (can_load) {
			auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
			std::ifstream f(filename);
			int id;
			while (f >> id) {
				I.DS.insert(id);
				l << " " << id;
			}
			l << "\n";
			if (I.DS.size() > before) {
				log << "Updated DS (cached): " << before << "+" << (I.DS.size() - before) << "="
					<< I.DS.size() << std::endl;
				return true;
			} else {
				log << "Cached file seems empty, discarding!" << std::endl;
			}
		}
	}

	log << "Will cache solution in " << filename << std::endl;
	std::filesystem::create_directory("cache");
	dump_sat(filename_sat, hclauses);
	return false;
}
