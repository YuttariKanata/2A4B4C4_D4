#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <chrono>

using json = nlohmann::json;
using u16  = std::uint16_t;
using u64  = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 DMAX = 1'000'000ULL;

// constexpr std::array<u64, 36> MODS = {
//     28561, 83521, 12167, 24389, 29791, 50653,
//     68921,  5041,  6241,  9409, 10201, 10609,
//     11881, 16129, 18769, 22201, 22801, 24649,
//     27889, 29929, 32761, 36481, 37249, 38809,
//     39601, 49729, 52441, 57121, 58081, 69169,
//     72361, 73441, 76729, 85849, 96721, 97969
// };

constexpr std::array<u64, 41> MODS = {
  15625, 16807, 28561, 83521, 12167,
  24389, 29791, 50653, 68921,103823,  2809,
   3721,  5041,  6241,  9409, 10201, 10609,
  11881, 16129, 18769, 22201, 22801, 24649,
  27889, 29929, 32761, 36481, 37249, 38809,
  39601, 49729, 52441, 57121, 58081, 69169,
  72361, 73441, 76729, 85849, 96721, 97969
};

constexpr std::array<u64, 4> RHO = {
    1, 182, 443, 624
};


// ============================================================
// x^4 mod m
// ============================================================

constexpr u64 pow4_mod(u64 x, u64 m) noexcept {
    const u64 x2 = ((x * x) % m);
    return ((x2 * x2) % m);
}


// ============================================================
// 32a^4 + b^4 の取り得る剰余テーブル
// ============================================================

void save_left_residue_tables(const std::vector<std::vector<bool>>& tables, const std::string& filename) {
    if (tables.size() != MODS.size()) {
        throw std::runtime_error("table count does not match MODS");
    }

    json root;
    root["version"] = 1;
    root["type"] = "left_residue_tables";

    root["moduli"] = json::array();

    for (u64 mod : MODS) {
        root["moduli"].push_back(mod);
    }

    root["tables"] = json::array();

    for (std::size_t i = 0; i < MODS.size(); ++i) {
        const u64 mod = MODS[i];
        const auto& table = tables[i];

        if (table.size() != mod) {
            throw std::runtime_error("table size does not match modulus");
        }

        json residues = json::array();

        for (u64 r = 0; r < mod; ++r) {
            if (table[r]) residues.push_back(r);
        }

        root["tables"].push_back({
            {"modulus", mod},
            {"residues", std::move(residues)}
        });
    }

    std::ofstream out(filename);

    if (!out) {
        throw std::runtime_error("failed to open " + filename);
    }

    out << root.dump(2) << '\n';
}

std::vector<std::vector<bool>> load_left_residue_tables(const std::string& filename) {
    std::ifstream in(filename);

    if (!in) {
        throw std::runtime_error("failed to open " + filename);
    }

    const json root = json::parse(in);

    if (!root.contains("version") || root["version"] != 1) {
        throw std::runtime_error("invalid residue table version");
    }

    if (!root.contains("type") || root["type"] != "left_residue_tables") {
        throw std::runtime_error("invalid residue table type");
    }

    if (!root.contains("moduli") || !root["moduli"].is_array()) {
        throw std::runtime_error("missing moduli");
    }

    if (root["moduli"].size() != MODS.size()) {
        throw std::runtime_error("MODS count mismatch");
    }

    for (std::size_t i = 0; i < MODS.size(); ++i) {
        const u64 file_mod = root["moduli"][i].get<u64>();

        if (file_mod != MODS[i]) {
            throw std::runtime_error("MODS mismatch at index " + std::to_string(i));
        }
    }

    if (!root.contains("tables") || !root["tables"].is_array()) {
        throw std::runtime_error("missing tables");
    }

    if (root["tables"].size() != MODS.size()) {
        throw std::runtime_error("table count mismatch");
    }

    std::vector<std::vector<bool>> tables;
    tables.reserve(MODS.size());

    for (std::size_t i = 0; i < MODS.size(); ++i) {
        const u64 mod = MODS[i];
        const auto& entry = root["tables"][i];

        if (!entry.contains("modulus") || entry["modulus"].get<u64>() != mod) {
            throw std::runtime_error("table modulus mismatch at index " + std::to_string(i));
        }

        if (!entry.contains("residues") || !entry["residues"].is_array()) {
            throw std::runtime_error("missing residues at index " + std::to_string(i));
        }

        std::vector<bool> table(mod, false);

        for (const auto& value : entry["residues"]) {
            const u64 r = value.get<u64>();

            if (r >= mod) {
                throw std::runtime_error("residue out of range at index " + std::to_string(i));
            }
            table[r] = true;
        }
        tables.push_back(std::move(table));
    }
    return tables;
}

std::vector<std::vector<bool>> make_left_residue_tables() {

    std::vector<std::vector<bool>> tables;
    tables.reserve(MODS.size());

    for (std::size_t i = 0; i < MODS.size(); ++i) {
        const u64 m = MODS[i];

        std::vector<bool> fourth(m, false);

        for (u64 x = 0; x < m; ++x) {
            fourth[pow4_mod(x, m)] = true;
        }

        std::vector<bool> left(m, false);

        for (u64 a4 = 0; a4 < m; ++a4) {
            if (!fourth[a4]) continue;

            for (u64 b4 = 0; b4 < m; ++b4) {
                if (!fourth[b4]) continue;

                const u64 r = (static_cast<u128>(32) * a4 + b4) % m;

                left[r] = true;
            }
        }

        tables.push_back(std::move(left));

        std::cout << m << " OK" << std::endl;
    }

    return tables;
}

std::vector<std::vector<bool>> get_left_residue_tables() {
    constexpr const char* filename = "left_residue_tables.json";

    {
        std::ifstream in(filename);

        if (in) {
            std::cout << "loading residue tables..." << std::endl;

            try {
                return load_left_residue_tables(filename);
            } catch (const std::exception& e) {
                std::cout << "invalid residue table: " << e.what() << std::endl;
                std::cout << "rebuilding residue tables..." << std::endl;
            }
        }
    }

    auto tables = make_left_residue_tables();

    std::cout << "saving residue tables..." << std::endl;

    save_left_residue_tables(tables, filename);

    return tables;
}

std::vector<bool> make_left_residue_table_65536() {

    const u64 m = 65536;

    std::vector<bool> fourth(m, false);

    for (u64 x = 0; x < m; ++x) {
        fourth[pow4_mod(x, m)] = true;
    }

    std::vector<bool> left(m, false);

    for (u64 a4 = 0; a4 < m; ++a4) {
        if (!fourth[a4]) continue;

        for (u64 b4 = 0; b4 < m; ++b4) {
            if (!fourth[b4]) continue;

            const u64 r = (static_cast<u128>(32) * a4 + b4) % m;

            left[r] = true;
        }
    }


    return left;
}

// ============================================================
// q mod 65536
//
// 65536 = 2^16 なので剰余は下位16bitだけ見ればよい。
// ============================================================

[[nodiscard]] inline constexpr u16 q_mod_65536(u128 q) noexcept {
    return static_cast<u16>(q);
}


// ============================================================
// q が mod m において左辺として可能でないか
// ============================================================

[[nodiscard]] inline bool is_not_q_residue_mod(u128 q, u64 mod, const std::vector<bool>& table) noexcept {
    return table[static_cast<u64>(q % mod)] == false;
}


// ============================================================
// q が全 modulus で左辺として可能か
// ============================================================

[[nodiscard]] inline bool is_q_residue(u128 q, const std::vector<std::vector<bool>>& tables, const std::vector<bool>& table_65536) noexcept {
    // mod 65536
    if (!table_65536[q_mod_65536(q)]) return false;

    // その他の素数
    for (std::size_t i = 0; i < MODS.size(); ++i) {
        if (is_not_q_residue_mod(q, MODS[i], tables[i])) return false;
    }

    return true;
}

[[nodiscard]] inline bool is_q_residue_65536(u128 q,const std::vector<bool>& table_65536) noexcept {
    // mod 65536
    if (!table_65536[q_mod_65536(q)]) return false;

    // その他の素数
    // for (std::size_t i = 0; i < MODS.size(); ++i) {
    //     if (is_not_q_residue_mod(q, MODS[i], tables[i])) return false;
    // }

    return true;
}


// ============================================================
// main
// ============================================================

int main()
{
    auto start = std::chrono::steady_clock::now();

    std::cout << "building residue tables..." << std::endl;

    const auto left_table_65536 = make_left_residue_table_65536();
    std::cout << "65536 OK" << std::endl;
    const auto left_tables = get_left_residue_tables();

    std::cout << "residue tables ready." << std::endl;

    std::ofstream out("q_table.bin", std::ios::binary);

    if (!out) {
        std::cerr << "failed to open q_table.bin" << std::endl;
        return 1;
    }

    u64 candidate_count = 0;
    u64 survived_count  = 0;

    for (u64 d = 1; d <= DMAX; d += 2) {
        if (d % 5 == 0) continue;

        const u128 d2 = static_cast<u128>(d) * d;
        const u128 d4 = d2 * d2;

        for (u64 rho : RHO) {
            const u64 c0 = (static_cast<u128>(d) * rho) % 625;

            for (u64 c = c0; c < d; c += 625) {
                if (c == 0) continue;

                ++candidate_count;

                const u128 c2 = static_cast<u128>(c) * c;
                const u128 c4 = c2 * c2;

                const u128 q = (d4 - c4) / 625;

                if (!is_q_residue(q, left_tables, left_table_65536)) continue;
                //if (!is_q_residue_65536(q, left_table_65536)) continue;

                ++survived_count;

                out.write(
                    reinterpret_cast<const char*>(&q),
                    sizeof(q)
                );
            }
        }

        if ((d+1) % 10000 == 0) {
            std::cout << "d = " << d
                      << " per = " << (static_cast<double>(survived_count) / candidate_count)
                    //   << "  candidates = " << candidate_count
                    //   << "  survived = " << survived_count
                      << std::endl;
        }
    }

    std::cout << "\n========== result ==========" << std::endl;
    std::cout << "candidate count: " << candidate_count << std::endl;
    std::cout << "survived count:  " << survived_count << std::endl;

    if (candidate_count != 0) {
        std::cout << "survival rate:   "
                  << static_cast<double>(survived_count) / static_cast<double>(candidate_count)
                  << std::endl;
    }

    std::cout << "file size:       "
              << static_cast<double>(survived_count * sizeof(u128))/(1024 * 1024) << "MB"
              << std::endl;

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "time: " << ms / 1000.0 << "seconds" << std::endl;

    return 0;
}