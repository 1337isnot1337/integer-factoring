// Exact finite F4 coefficient enumeration at inherited integer physical heights.
// The included polynomial engine has a fixed precision window. An unresolved
// promised rank is recorded, never counted as a successful profile.
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include "native_weighted_smith_f4.cpp"

namespace {
constexpr std::uint32_t total_codes = 1u << 22;
struct Record {
    std::uint64_t count = 0;
    std::uint32_t representative = 0;
    std::uint64_t lo = 0, hi = 0;
};

std::uint32_t number(const char* value) {
    std::string text(value);
    std::size_t consumed = 0;
    auto result = std::stoull(text, &consumed);
    if (consumed != text.size() || result > total_codes)
        throw std::runtime_error("Invalid coefficient enumeration argument");
    return std::uint32_t(result);
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::runtime_error("Use START STOP ODD_MULTIPLIER OFFSET and the PENCIL1 input on stdin");
        auto start = number(argv[1]), stop = number(argv[2]);
        auto multiplier = number(argv[3]), offset = number(argv[4]);
        if (start > stop || stop > total_codes || !(multiplier & 1)
            || multiplier >= total_codes || offset >= total_codes)
            throw std::runtime_error("Invalid coefficient permutation or interval");
        std::string magic;
        std::cin >> magic;
        if (magic != "PENCIL1") throw std::runtime_error("Unknown pencil input format");
        std::array<std::uint64_t, 1250> base{}, values{};
        std::array<unsigned, 625> first{}, second{};
        for (int i = 0; i < 625; ++i) {
            std::cin >> base[2*i] >> base[2*i+1] >> first[i] >> second[i];
            if (first[i] > 1 || second[i] > 1) throw std::runtime_error("Invalid pencil direction bit");
        }
        std::array<double, 33> parameters{};
        for (auto& value : parameters) {
            std::cin >> value;
            if (!std::isfinite(value) || value != std::round(value))
                throw std::runtime_error("The enumeration requires integer profile weights");
        }
        if (!std::cin || parameters[32] != 1.) throw std::runtime_error("Truncated pencil input or nonunit scale");
        std::map<std::array<int, 50>, Record> classes;
        std::uint64_t unresolved = 0, endpoint_failures = 0;
        for (std::uint32_t position = start; position < stop; ++position) {
            auto code = std::uint32_t((std::uint64_t(multiplier)*position + offset) & (total_codes-1));
            std::uint64_t lo = 0, hi = 0;
            for (int i = 0; i < 11; ++i) {
                lo |= std::uint64_t((code >> (2*i)) & 1) << i;
                hi |= std::uint64_t((code >> (2*i+1)) & 1) << i;
            }
            lo ^= lo << 1;
            hi ^= hi << 1;
            for (int i = 0; i < 625; ++i) {
                values[2*i] = base[2*i] ^ (first[i] ? lo : 0) ^ (second[i] ? lo << 4 : 0);
                values[2*i+1] = base[2*i+1] ^ (first[i] ? hi : 0) ^ (second[i] ? hi << 4 : 0);
            }
            int status = gf4_set_seed(5, values.data(), 1, 1);
            if (status) { ++endpoint_failures; continue; }
            std::array<int, 50> orders{};
            bool resolved = true;
            for (int side = 0; side < 2; ++side) {
                auto pivots = f4::profile(parameters.data() + 16*side,
                                          side == 0 ? 1 : 0, side == 0, 1.);
                if (pivots.size() != 25) { resolved = false; break; }
                for (int i = 0; i < 25; ++i) {
                    double value = pivots[i].value;
                    if (std::abs(value - std::round(value)) > 1e-9)
                        throw std::runtime_error("A supposedly integer profile acquired a fractional order");
                    orders[25*side+i] = int(std::llround(value));
                    if (i && orders[25*side+i] < orders[25*side+i-1])
                        throw std::runtime_error("A profile was not sorted");
                }
            }
            if (!resolved) { ++unresolved; continue; }
            int zero_sum = 0;
            for (int i = 25; i < 50; ++i) zero_sum += orders[i];
            if (zero_sum != 444) throw std::runtime_error("The conjugation determinant invariant failed");
            auto& record = classes[orders];
            if (!record.count) { record.representative = code; record.lo = lo; record.hi = hi; }
            ++record.count;
        }
        std::cout << "{\"start\":" << start << ",\"stop\":" << stop
                  << ",\"multiplier\":" << multiplier << ",\"offset\":" << offset
                  << ",\"visited\":" << stop-start << ",\"unresolved\":" << unresolved
                  << ",\"endpoint_failures\":" << endpoint_failures << ",\"classes\":[";
        bool first_record = true;
        for (const auto& [orders, record] : classes) {
            if (!first_record) std::cout << ',';
            first_record = false;
            std::cout << "{\"count\":" << record.count << ",\"representative\":" << record.representative
                      << ",\"q_planes\":[" << record.lo << ',' << record.hi << "],\"orders\":[";
            for (int i = 0; i < 50; ++i) { if (i) std::cout << ','; std::cout << orders[i]; }
            std::cout << "]}";
        }
        std::cout << "]}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
