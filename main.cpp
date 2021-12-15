#include <iostream>
#include <cmath>
#include <vector>


struct tridiagonal_matrix {
private:
    std::vector<double> top, mid, bot;
public:
    tridiagonal_matrix() {}

    void resize(size_t n) {
        mid.resize(n, 0);
        bot.resize(n-1, 0);
        top.resize(n-1, 0);
    }

    size_t size() {
        return mid.size();
    }

    double &operator()(size_t i, size_t j) {
        if (i == j) {
            return mid[i];
        }
        else if (j == i-1) {
            return bot[j];
        }
        else if (j == i+1) {
            return top[i];
        }
        return *(double *)nullptr;
    }

    double get(size_t i, size_t j) {
        if (i - j + 1 > 2) {  // unsigned
            return 0;
        }
        return (*this)(i, j);
    }
};


namespace var {
    double u1 = 0, u2 = 1;
    double x0 = 0.525;
    double eps = 0.01;

    double k(double x, bool first = true) {
        if (x < x0 || (x == x0 && first)) {
            return std::exp(-x*x);
        }
        return x;
    }

    double q(double x) {
        return x*x;
    }

    double f(double x) {
        return std::sin(x);
    }
}


// using namespace var implicitly -- not pure
void create_mat_vec(size_t n, tridiagonal_matrix &mat, std::vector<double> &vec) {
    // mat @ solution = vec
}


void solve(tridiagonal_matrix &mat, std::vector<double> &vec,
        std::vector<double> &solution) {
    size_t n = vec.size();
    solution.resize(n, 0);

    std::vector<double> alpha(n, 0.0/0.0), beta(n, 0.0/0.0);

    alpha[1] = -mat(0, 1) / mat(0, 0);
    beta[1] = vec[0] / mat(0, 0);
    for (size_t i = 1; i < n-1; ++i) {
        alpha[i+1] = -mat(i, i+1) / (mat(i, i-1) * alpha[i] + mat(i, i));
        beta[i+1] = (vec[i] - mat(i, i-1) * beta[i]) / (mat(i, i-1) * alpha[i] + mat(i, i));
    }

    solution[n-1] = (vec[n-1] - mat(n-1, n-2) * beta[n-1])
        / (mat(n-1, n-1) + mat(n-1, n-2) * alpha[n-1]);

    for (size_t i = n-2; i < n-1; --i) {
        solution[i] = alpha[i+1] * solution[i+1] + beta[i+1];
    }
}


void test_solve() {
    tridiagonal_matrix mat;
    std::vector<double> vec, solution;

    size_t n = 4;
    std::cin >> n;
    create_mat_vec(n, mat, vec);
    solve(mat, vec, solution);

    for (size_t i = 0; i < n; ++i) {
        std::cout << "[";
        for (size_t j = 0; j < n; ++j) {
            std::cout << mat.get(i, j) << ", ";
        }
        std::cout << "],\n";
    }

    for (size_t i = 0; i < n; ++i) {
        std::cout << vec[i] << ", ";
    }
    std::cout << "\n";

    for (size_t i = 0; i < n; ++i) {
        std::cout << solution[i] << ", ";
    }
    std::cout << std::endl;
}

int main() {
    test_solve();
}
