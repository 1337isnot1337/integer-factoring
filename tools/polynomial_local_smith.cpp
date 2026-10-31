// Local Smith orders over F2[[t]], using exact truncated polynomial elimination.
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using boost::multiprecision::cpp_int;

cpp_int precision_mask(unsigned precision) {
    cpp_int result = 1;
    result <<= precision;
    --result;
    return result;
}

cpp_int multiply(cpp_int a, cpp_int b, const cpp_int &mask) {
    cpp_int result = 0;
    while (b != 0) {
        unsigned bit = boost::multiprecision::lsb(b);
        result ^= a << bit;
        b &= b-1;
    }
    return result & mask;
}

cpp_int unit_inverse(const cpp_int &unit, unsigned precision) {
    if ((unit & 1) == 0) throw std::runtime_error("nonunit pivot");
    cpp_int inverse = 1;
    for (unsigned known = 1; known < precision; known *= 2) {
        unsigned next = std::min(2*known, precision);
        cpp_int mask = precision_mask(next);
        // In characteristic two, Newton's inverse update is u * inverse^2.
        inverse = multiply(unit & mask, multiply(inverse, inverse, mask), mask);
    }
    return inverse;
}

int main() {
    unsigned rows, columns, precision;
    while (std::cin >> rows >> columns >> precision) {
        if (rows < columns || columns == 0 || precision == 0 || precision > 8192)
            throw std::runtime_error("invalid dimensions or precision");
        const cpp_int mask = precision_mask(precision);
        std::vector<std::vector<cpp_int>> matrix(rows, std::vector<cpp_int>(columns));
        for (auto &row : matrix) for (auto &entry : row) {
            std::string token;
            if (!(std::cin >> token) || token.empty()) throw std::runtime_error("invalid polynomial input");
            entry = 0;
            for (char digit : token) {
                if (digit < '0' || digit > '9') throw std::runtime_error("invalid polynomial input");
                entry *= 10;
                entry += digit-'0';
            }
            entry &= mask;
        }
        std::vector<unsigned> orders;
        for (unsigned k = 0; k < columns; ++k) {
            unsigned minimum = precision, pivot_row = k, pivot_column = k;
            for (unsigned i = k; i < rows; ++i) for (unsigned j = k; j < columns; ++j) {
                if (matrix[i][j] == 0) continue;
                unsigned valuation = boost::multiprecision::lsb(matrix[i][j]);
                if (valuation < minimum) {
                    minimum = valuation;
                    pivot_row = i;
                    pivot_column = j;
                }
            }
            if (minimum == precision) break;
            std::swap(matrix[k], matrix[pivot_row]);
            for (auto &row : matrix) std::swap(row[k], row[pivot_column]);
            orders.push_back(minimum);
            cpp_int inverse = unit_inverse(matrix[k][k] >> minimum, precision-minimum);
            cpp_int quotient_mask = precision_mask(precision-minimum);
            for (unsigned i = k+1; i < rows; ++i) {
                if (matrix[i][k] == 0) continue;
                if (boost::multiprecision::lsb(matrix[i][k]) < minimum)
                    throw std::runtime_error("pivot did not have minimum valuation");
                cpp_int quotient = multiply(matrix[i][k] >> minimum, inverse, quotient_mask);
                if (multiply(quotient, matrix[k][k], mask) != matrix[i][k])
                    throw std::runtime_error("pivot division failed");
                for (unsigned j = k+1; j < columns; ++j)
                    matrix[i][j] ^= multiply(quotient, matrix[k][j], mask);
                matrix[i][k] = 0;
            }
        }
        std::cout << "{\"precision\":" << precision << ",\"complete\":"
                  << (orders.size() == columns ? "true" : "false") << ",\"orders\":[";
        for (unsigned i = 0; i < orders.size(); ++i) std::cout << (i ? "," : "") << orders[i];
        std::cout << "]}\n";
    }
}
