#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>

using u16  = std::uint16_t;
using u64  = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 DMAX = 400'000ULL;

// constexpr std::array<u64, 36> MODS = {
//     28561, 83521, 12167, 24389, 29791, 50653,
//     68921,  5041,  6241,  9409, 10201, 10609,
//     11881, 16129, 18769, 22201, 22801, 24649,
//     27889, 29929, 32761, 36481, 37249, 38809,
//     39601, 49729, 52441, 57121, 58081, 69169,
//     72361, 73441, 76729, 85849, 96721, 97969
// };

constexpr std::array<u64, 42> MODS = {
    625,  2197,  2401,  2197,  4913,   529,
    841,   961,  1369,  1681,  2209,  2809,
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

void save_left_residue_tables(
    const std::vector<std::vector<bool>>& tables,
    const std::string& filename
) {
    std::ofstream out(filename, std::ios::binary);

    if (!out) {
        throw std::runtime_error("failed to open " + filename);
    }

    const std::uint64_t count = tables.size();
    out.write(
        reinterpret_cast<const char*>(&count),
        sizeof(count)
    );

    for (const auto& table : tables) {
        const std::uint64_t size = table.size();

        out.write(
            reinterpret_cast<const char*>(&size),
            sizeof(size)
        );

        const std::size_t byte_count = (size + 7) / 8;

        std::vector<std::uint8_t> bytes(byte_count, 0);

        for (std::size_t i = 0; i < size; ++i) {
            if (table[i]) {
                bytes[i >> 3] |= static_cast<std::uint8_t>(1u << (i & 7));
            }
        }

        out.write(
            reinterpret_cast<const char*>(bytes.data()),
            bytes.size()
        );
    }
}

std::vector<std::vector<bool>> load_left_residue_tables(
    const std::string& filename
) {
    std::ifstream in(filename, std::ios::binary);

    if (!in) {
        throw std::runtime_error("failed to open " + filename);
    }

    std::uint64_t count;

    in.read(
        reinterpret_cast<char*>(&count),
        sizeof(count)
    );

    std::vector<std::vector<bool>> tables;
    tables.reserve(count);

    for (std::uint64_t j = 0; j < count; ++j) {
        std::uint64_t size;

        in.read(
            reinterpret_cast<char*>(&size),
            sizeof(size)
        );

        std::vector<bool> table(size, false);

        const std::size_t byte_count = (size + 7) / 8;
        std::vector<std::uint8_t> bytes(byte_count);

        in.read(
            reinterpret_cast<char*>(bytes.data()),
            bytes.size()
        );

        for (std::size_t i = 0; i < size; ++i) {
            table[i] = (bytes[i >> 3] >> (i & 7)) & 1;
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

std::vector<bool> make_left_residue_table_65536() {

    const u64 m = 32768;

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

std::vector<std::vector<bool>> get_left_residue_tables() {
    constexpr const char* filename = "left_residue_tables.bin";

    {
        std::ifstream in(filename, std::ios::binary);

        if (in) {
            std::cout << "loading residue tables...\n";
            return load_left_residue_tables(filename);
        }
    }

    std::cout << "building residue tables...\n";

    auto tables = make_left_residue_tables();

    std::cout << "saving residue tables...\n";

    save_left_residue_tables(tables, filename);

    return tables;
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
              << survived_count * sizeof(u128)
              << std::endl;

    return 0;
}