#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

// エラトステネスの篩でN以下の素数を列挙する関数
std::vector<int> get_primes(int n) {
    std::vector<bool> is_prime(n + 1, true);
    std::vector<int> primes;
    if (n >= 2) is_prime[0] = is_prime[1] = false;
    for (int p = 2; p * p <= n; p++) {
        if (is_prime[p]) {
            for (int i = p * p; i <= n; i += p)
                is_prime[i] = false;
        }
    }
    for (int p = 2; p <= n; p++) {
        if (is_prime[p]) primes.push_back(p);
    }
    return primes;
}

// 繰り返し二乗法によるベキ剰余計算 (x^4 % p)
long long powermod(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return res;
}

int main() {
    const int N = 65536; // 上限値

    // 各素数 p に対し、p^k <= N である限り p^k を生成
    std::vector<int> prime_powers;
    std::vector<int> primes = get_primes(N);

    for (int p : primes) {
        long long val = p;
        while (val <= N) {
            prime_powers.push_back(static_cast<int>(val));
            // オーバーフロー防止対策を含めた安全な乗算
            if (val > N / p) {
                break;
            }
            val *= p;
        }
    }

    // 昇順に並び替え
    std::sort(prime_powers.begin(), prime_powers.end());

    int cnt = 0;

    for (int p : prime_powers) {
        // Juliaの falses(p) に相当 (0 から p-1 までの要素)
        std::vector<bool> r4(p, false);
        for (int x = 0; x <= p / 2; ++x) {
            long long rem = powermod(x, 4, p);
            r4[rem] = true; // C++は0始まりなのでそのままインデックスに使える
        }

        // findall(r4) の代わりとして、trueのインデックス（0始まり）を抽出
        std::vector<int> r4_indices;
        for (int i = 0; i < p; ++i) {
            if (r4[i]) {
                r4_indices.push_back(i);
            }
        }

        std::vector<bool> rq(p, false);
        for (int a : r4_indices) {
            for (int b : r4_indices) {
                // Julia: ((32*a + b) % p) + 1
                // C++は0始まりなので、余りをそのままインデックスとして使用
                long long idx = (32LL * a + b) % p;
                rq[idx] = true;
            }
        }

        // count(rq) の計算
        int rq_count = 0;
        for (bool val : rq) {
            if (val) rq_count++;
        }

        if (rq_count != p) {
            // 出力フォーマットの調整 (lpadの再現)
            double ratio = std::round((static_cast<double>(rq_count) / p) * 100000.0) / 100000.0;
            
            std::cout << std::setw(6) << p 
                      << std::setw(10) << std::fixed << std::setprecision(5) << ratio 
                      << " | ";

            cnt++;

            if (cnt % 11 == 0) {
                std::cout << "\n";
            }
        }
    }
    
    // 最後に改行が残っていれば出力
    if (cnt % 11 != 0) {
        std::cout << "\n";
    }

    return 0;
}
