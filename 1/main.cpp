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
#ifdef TEST
    bool model = true;
#else
    bool model = false;
#endif

    double k(double x) {
        if (model) {
            return 0.5*0.5 + 1;
        }
        return x*x + 1;
    }

    double q(double x) {
        if (model) {
            return 0.5;
        }
        return x;
    }

    double f(double x) {
        if (model) {
            return std::exp(-0.5);
        }
        return std::exp(-x);
    }

    double model_u(double x) {
        double s10 = std::sqrt(10);
        return ((8 - 2*s10) * std::exp(s10 * x/5)
                + (2*s10 + 8) * std::exp(s10*(x + 2)/5)
                - 8 * std::exp(s10 * (2*x + 1)/5) - 8 * std::exp(s10/5))
            * std::exp(-s10*x/5 - 0.5)
            / (-s10 + 4 + (s10 + 4) * std::exp(2*s10/5));
    }

    double beta_1 = 0, mu_1 = 0;
    double beta_2 = 1, mu_2 = 0;
    double eps = 0.01;
}


void create_mat_vec(size_t n, tridiagonal_matrix &mat, std::vector<double> &vec) {
    mat.resize(n);
    vec.resize(n);
    double h = 1.0 / (n-1);

    mat(0, 0) = -var::k(h/2) / h - (var::beta_1 + h * var::q(0) / 2);
    mat(0, 1) = var::k(h/2) / h;
    vec[0] = -(var::mu_1 + h * var::f(0) / 2);

    mat(n-1, n-1) = -var::k(1-h/2) / h - (var::beta_2 + h * var::q(1) / 2);
    mat(n-1, n-2) = var::k(1-h/2) / h;
    vec[n-1] = -(var::mu_2 + h * var::f(1) / 2);

    for (size_t i = 1; i < n-1; ++i) {
        mat(i, i-1) = var::k(i*h - h/2) / (h*h);
        mat(i, i) = -(var::k(i*h+h/2) + var::k(i*h-h/2)) / (h*h) - var::q(i*h);
        mat(i, i+1) = var::k(i*h + h/2) / (h*h);
        vec[i] = -var::f(i*h);
    }
}


// works only if vec.size() > 2
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


// vn.size() * 2 should be equal to v2n.size()
double calc_error(std::vector<double> &vn, std::vector<double> &v2n) {
    double ans = 0;
    for (size_t i = 0; i < vn.size(); ++i) {
        double err = std::abs(vn[i] - v2n[2*i]);
        ans = std::max(err, ans);
        //ans += err*err;
    }
    return ans;
}


void runge_rule(double eps, std::vector<double> &sol2) {
    size_t n = 3;
    double err = -1;
    tridiagonal_matrix mat;
    std::vector<double> vec, solution;

    create_mat_vec(n, mat, vec);
    solve(mat, vec, sol2);

    do {
        if (err > 0) {
            std::cerr << "[Runge error " << err << " " << sol2.size() << "]" << std::endl;
        }
        n *= 2;
        solution.swap(sol2);
        create_mat_vec(n, mat, vec);
        solve(mat, vec, sol2);
    } while ((err = calc_error(solution, sol2)) > eps);
    std::cerr << "[Runge error " << err << " " << sol2.size() << "]" << std::endl;
}


int test_solve(size_t n, bool print) {
    tridiagonal_matrix mat;
    std::vector<double> vec, solution;

    if (n < 2) {
        std::cerr << "n<2 not supported" << std::endl;
        return 1;
    }
    create_mat_vec(n, mat, vec);
    solve(mat, vec, solution);

    if (print) {
        // print matrix
        for (size_t i = 1; i < n; ++i) {
            std::cout << mat.get(i-1, i) << ", ";
        }
        std::cout << "\n";

        for (size_t i = 0; i < n; ++i) {
            std::cout << mat.get(i, i) << ", ";
        }
        std::cout << "\n";

        for (size_t i = 1; i < n; ++i) {
            std::cout << mat.get(i, i-1) << ", ";
        }
        std::cout << "\n";

        // print right side
        for (size_t i = 0; i < n; ++i) {
            std::cout << vec[i] << ", ";
        }
        std::cout << "\n";

        // print solution
        for (size_t i = 0; i < n; ++i) {
            std::cout << solution[i] << ", ";
        }
        std::cout << std::endl;
    }
    else {
        double me = 0;
        for (size_t i = 0; i < n; ++i) {
            double err = std::abs(var::model_u(i * 1.0 / (n-1)) - solution[i]);
            if (me < err) {
                std::cerr << "{{" << i << " " << var::model_u(i / (n-1)) << " " << solution[i] << "}}\n";
                me = err;
            }
        }
        std::cerr << "[error = " << me << "]" << std::endl;
    }
    return 0;
}

int main(int argc, char **args) {
    std::cout.precision(16);
#ifdef TEST
    size_t n;
    std::cin >> n;
    if (args == nullptr) {
        return 1;
    }
    return test_solve(n, argc <= 1);
#endif

    std::vector<double> ans;
    runge_rule(var::eps, ans);

    for (auto i : ans) {
        std::cout << i << ", ";
    }
    std::cout << std::endl;
}
