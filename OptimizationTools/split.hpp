#pragma once
#include<string>
#include<cctype>

struct __isspace {
	bool operator()(int ch) const {
		return (std::isspace(ch) != 0);
	}
};

namespace std {
	template<class Container, typename Delimiters = ::__isspace>
	void split(const std::string& str, Container& splits, Delimiters delim = Delimiters{}, const bool trim = true) {
		static_assert(std::is_same<typename Container::value_type, std::string>::value, "Container type is not compatible!");
		using iterator = std::string::const_iterator;

		iterator itr = str.begin(), end = str.end();
		std::string temp_split;
		if (trim) {
			while (true) {
				while (true) {
					if (itr == end) return;
					if (!delim(*itr)) break;
					++itr;
				}
				while (true) {
					temp_split += *itr;
					if (++itr == end) {
						splits.push_back(temp_split);
						return;
					}
					if (delim(*itr)) break;
				}
				splits.push_back(temp_split);
				temp_split.clear();
			}
		}
		else {
			while (true) {
				while (true) {
					if (itr == end) {
						splits.push_back(temp_split);
						return;
					}
					if (delim(*itr)) break;
					else
						temp_split += *itr;
					++itr;
				}
				splits.push_back(temp_split);
				++itr;
				temp_split.clear();
			}
		}
	}
}
