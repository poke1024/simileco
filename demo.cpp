// Author: Bernhard Liebl, 2020
// Released under a MIT license.

#include <iostream>

#include "aligner.h"

int main() {

    // the example from the README.

    Aligner aligner(20, 20);

    const std::string s("CHOCOLATEISTHEANSWER");
    const std::string t("LATETHAW");

    const auto similarity = [&s, &t] (int i, int j) {
        return s[i] == t[j] ? 1 : -1;
    };

    const auto gap_cost = [] (int n) {
        return std::pow(1.25, n);
    };

    aligner.waterman_smith_beyer(
        similarity,
        gap_cost,
        s.length(),
        t.length());

    std::cout << aligner.pretty_printed(s, t);

    return 0;
}