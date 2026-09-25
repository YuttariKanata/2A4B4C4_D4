#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using u64  = std::uint64_t;
using u128 = __uint128_t;

struct Solution {
    u64 A;
    u64 B;
    u64 C;
    u64 D;
};

// rho^4 % 625 = 1
constexpr u64 RHO[] = {
    1, 182, 443, 624
};

// ------------------------------------------------------------
// integer fourth root
// returns floor(x^(1/4))
// ------------------------------------------------------------

u64 fourth_root(u128 x)
{
    // 簡易版。
    // Phase 1.0 では correctness を優先する。
    u64 lo = 0;
    u64 hi = 1ULL << 32;

    while (lo + 1 < hi) {
        u64 mid = lo + (hi - lo) / 2;

        u128 m = u128(mid);
        u128 m4 = m * m * m * m;

        if (m4 <= x)
            lo = mid;
        else
            hi = mid;
    }

    return lo;
}

bool is_fourth_power(u128 x, u64& root)
{
    root = fourth_root(x);

    u128 r = u128(root);
    return r * r * r * r == x;
}

// ------------------------------------------------------------
// Solve
//
// 32 a^4 + b^4 = q
//
// ------------------------------------------------------------

bool solve_q(u128 q, u64& a_sol, u64& b_sol)
{
    // 32 a^4 <= q - 1
    //
    // ここも最初は単純な二分探索。
    u64 lo = 0;
    u64 hi = 100'000'000ULL + 1;

    while (lo + 1 < hi) {
        u64 mid = lo + (hi - lo) / 2;

        u128 a = u128(mid);
        u128 a4 = a * a * a * a;

        if (u128(32) * a4 <= q - 1)
            lo = mid;
        else
            hi = mid;
    }

    const u64 amax = lo;

    for (u64 a = 1; a <= amax; ++a) {
        u128 aa = u128(a);
        u128 a4 = aa * aa * aa * aa;

        u128 r = q - u128(32) * a4;

        u64 b;
        if (is_fourth_power(r, b)) {
            a_sol = a;
            b_sol = b;
            return true;
        }
    }

    return false;
}

// ------------------------------------------------------------
// Search
// ------------------------------------------------------------

void search(u64 Dmax)
{

    for (u64 d = 1; d <= Dmax; ++d) {

        // D is odd and 5 ∤ D.
        if ((d & 1) == 0) continue;

        if (d % 5 == 0) continue;

        for (u64 rho : RHO) {

            // c ≡ d*rho (mod 625)
            u64 residue = (u64)((u128(d) * rho) % 625);

            // c = residue + 625*k
            //
            // residue = 0 はこの場合発生しないが、
            // 一般形として処理しておく。
            u64 c = residue;

            if (c == 0) c = 625;

            for (; c < d; c += 625) {

                // q = (d^4-c^4)/625

                u128 dd = u128(d);
                u128 cc = u128(c);

                u128 d4 = dd * dd * dd * dd;
                u128 c4 = cc * cc * cc * cc;

                u128 diff = d4 - c4;
                u128 q = diff / 625;

                // 32a^4+b^4=q
                u64 a;
                u64 b;

                if (!solve_q(q, a, b)) continue;

                const u64 A = 10 * a;
                const u64 B = 5 * b;
                const u64 C = c;

                // 念のため元の式を直接検証
                u128 AA = u128(A);
                u128 BB = u128(B);
                u128 CC = u128(C);
                u128 DD = u128(d);

                u128 lhs = u128(2) * AA * AA * AA * AA + BB * BB * BB * BB + CC * CC * CC * CC;

                u128 rhs = DD * DD * DD * DD;

                if (lhs != rhs) std::terminate();

                const u64 g1 = std::gcd(A, B);
                const u64 g2 = std::gcd(C, d);
                const u64 g  = std::gcd(g1, g2);

                if (g != 1) continue;

                std::cout
                    << "solution: "
                    << "A=" << A
                    << ", B=" << B
                    << ", C=" << C
                    << ", D=" << d
                    << '\n';
            }
        }
    }
}

int main() {
    std::cout << "go!" << std::endl;
    search(10000ULL);
    std::cout << "\nend!" << std::endl;

}