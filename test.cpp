// Author: Bernhard Liebl, 2020
// Released under a MIT license.

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "aligner.h"

#include <iostream>
#include <unordered_map>

class BinarySimilarity {
private:
	const std::string _a;
	const std::string _b;
	const int _match;
	const int _mismatch;

public:
	BinarySimilarity(
		const std::string &a,
		const std::string &b,
		int match = 1,
		int mismatch = -1) :

		_a(a), _b(b), _match(match), _mismatch(mismatch) {

	}

	int operator()(int i, int j) const {
		return _a[i] == _b[j] ? _match : _mismatch;
	}
};

class DNAFull {
private:
	// see http://rosalind.info/glossary/dnafull/

	const char _alphabet[15][2] = {
		"A", "T", "G", "C", "S", "W", "R", "Y",
		"K", "M", "B", "V", "H", "D", "N"};

	std::unordered_map<char, int> _alphabet_lut;

	const int _scores[15][15] = {
		{5,-4,-4,-4,-4,1,1,-4,-4,1,-4,-1,-1,-1,-2},
		{-4,5,-4,-4,-4,1,-4,1,1,-4,-1,-4,-1,-1,-2},
		{-4,-4,5,-4,1,-4,1,-4,1,-4,-1,-1,-4,-1,-2},
		{-4,-4,-4,5,1,-4,-4,1,-4,1,-1,-1,-1,-4,-2},
		{-4,-4,1,1,-1,-4,-2,-2,-2,-2,-1,-1,-3,-3,-1},
		{1,1,-4,-4,-4,-1,-2,-2,-2,-2,-3,-3,-1,-1,-1},
		{1,-4,1,-4,-2,-2,-1,-4,-2,-2,-3,-1,-3,-1,-1},
		{-4,1,-4,1,-2,-2,-4,-1,-2,-2,-1,-3,-1,-3,-1},
		{-4,1,1,-4,-2,-2,-2,-2,-1,-4,-1,-3,-3,-1,-1},
		{1,-4,-4,1,-2,-2,-2,-2,-4,-1,-3,-1,-1,-3,-1},
		{-4,-1,-1,-1,-1,-3,-3,-1,-1,-3,-1,-2,-2,-2,-1},
		{-1,-4,-1,-1,-1,-3,-1,-3,-3,-1,-2,-1,-2,-2,-1},
		{-1,-1,-4,-1,-3,-1,-3,-1,-3,-1,-2,-2,-1,-2,-1,},
		{-1,-1,-1,-4,-3,-1,-1,-3,-1,-3,-2,-2,-2,-1,-1},
		{-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}};

	const std::string _a;
	const std::string _b;

public:
	DNAFull(
		const std::string &a,
		const std::string &b) : _a(a), _b(b) {

		for (int i = 0; i < 15; i++) {
			_alphabet_lut[_alphabet[i][0]] = i;
		}
	}

	int operator()(int i, int j) const {
		return _scores[_alphabet_lut.at(_a[i])][_alphabet_lut.at(_b[j])];
	}
};

class AffineGap {
private:
	const int _opening;
	const int _extension;

public:
	AffineGap(int opening, int extension) :
		_opening(opening), _extension(extension) {

	}

	int operator()(int n) const {
		if (n > 0) {
			return _opening + (n - 1) * _extension;
		} else {
			return 0;
		}
	}
};

std::string trim(const std::string &s) {
	int i = 0, j = s.length();

	while (i < s.length() && isspace(s[i])) {
		i++;
	}
	while (j > i && isspace(s[j - 1])) {
		j--;
	}

	return s.substr(i, j - i);
}

std::string normalize(const std::string &s) {
	std::istringstream stream(s);
	std::string line;
	std::ostringstream r;

	int i = 0;
	while (std::getline(stream, line)) {
		std::string s = trim(line);
		if (s.length() > 0) {
			if (i > 0) {
				r << "\n";
			}
			r << s;
			i++;
		}
	}

	return r.str();
}

TEST_CASE( "Needleman-Wunsch", "[nw]" ) {

	SECTION( "AAGDAXSFXAF, GDSXFF" ) {

		Aligner aligner(20, 20);

		const std::string s("AAGDAXSFXAF");
		const std::string t("GDSXFF");

		aligner.needleman_wunsch(
			BinarySimilarity(s, t),
			2,
			s.length(),
			t.length());

		REQUIRE( normalize(aligner.pretty_printed(s, t)) == normalize(R"ALIGNMENT(
			AAGDAXSFXAF
			||  | |||
			--GD--S-XFF
			)ALIGNMENT") );

		REQUIRE( aligner.score() == -6 );

	}
}

TEST_CASE( "Smith-Waterman", "[sw]" ) {

	SECTION( "AAGDAXSFXAF, GDSXFF" ) {

		Aligner aligner(20, 20);

		const std::string s("AAGDAXSFXAF");
		const std::string t("GDSXFF");

		aligner.smith_waterman(
			BinarySimilarity(s, t, 2, -2),
			1,
			s.length(),
			t.length());

		REQUIRE( normalize(aligner.pretty_printed(s, t)) == normalize(R"ALIGNMENT(
			AAGDAXSFXAF-
			||  | | | 
			--GD--S-X-FF
			)ALIGNMENT") );

		REQUIRE( aligner.score() == 6);

	}

	SECTION( "en.wikipedia.org, DNAfull, linear gap 1" ) {

		Aligner aligner(20, 20);

		const std::string s("TACGGGCCCGCTAC");
		const std::string t("TAGCCCTATCGGTCA");

		aligner.smith_waterman(
			DNAFull(s, t),
			1,
			s.length(),
			t.length());

			REQUIRE( normalize(aligner.pretty_printed(s, t)) == normalize(R"ALIGNMENT(
				TACGGGCCCGCTAAC-----
				||   | || |||||     
				TA---G-CC-CTATCGGTCA
				)ALIGNMENT") );

    }

	SECTION( "en.wikipedia.org, DNAfull, affine gap 5, 1" ) {

		Aligner aligner(20, 20);

		const std::string s("TACGGGCCCGCTAC");
		const std::string t("TAGCCCTATCGGTCA");

		aligner.waterman_smith_beyer(
			DNAFull(s, t),
			AffineGap(5, 1),
			s.length(),
			t.length());

		REQUIRE( normalize(aligner.pretty_printed(s, t)) == normalize(R"ALIGNMENT(
			TACGGGCCCGCTA-------C
			||   |||  |||        
			TA---GCC--CTATCGGTCA-	
			)ALIGNMENT") );

	}

}
