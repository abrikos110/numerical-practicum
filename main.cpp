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
    double k(double x) {
        return 0.5*0.5 + 1;
        return x*x + 1;
    }

    double q(double x) {
        return 0.5;
        return x;
    }

    double f(double x) {
        return std::exp(-0.5);
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

    // k(0) u'(0) = beta_1 u(0) - mu_1
    // -k(1) u'(1) = beta_2 u(1) - mu_2
    // u'(1) = -0.5 u(1)

    // model : k = k(0.5), q = q(0.5), f = f(0.5)
}


void create_mat_vec(size_t n, tridiagonal_matrix &mat, std::vector<double> &vec) {
    // mat @ solution = vec
    mat.resize(n);
    vec.resize(n);
    double h = 1.0 / (n-1);

    for (size_t i = 1; i < n-1; ++i) {
        double kk = var::k(i*h);
        double kd = (var::k(i*h + h) - var::k(i*h - h)) / (2*h);
        mat(i, i-1) = -kd * h / 2 + kk;
        mat(i, i) = -2 * kk - h*h * var::q(i*h);
        mat(i, i+1) = kd * h / 2 + kk;
        vec[i] = -var::f(i*h) * h*h;
        /*
        mat(i, i-1) = kk;
        mat(i, i) = -kk - var::k(i*h + h) - var::q(i*h) * h*h;
        mat(i, i+1) = var::k(i*h + h);
        vec[i] = -var::f(i*h)*h*h;*/
    }

    // first order
    /*mat(0, 0) = -var::k(0) / h - var::beta_1;
    mat(0, 1) = var::k(0) / h;
    vec[0] = -var::mu_1;

    mat(n-1, n-1) = -var::k(1) / h - var::beta_2;
    mat(n-1, n-2) = var::k(1) / h;
    vec[n-1] = -var::mu_2;*/


    // second order
    mat(0, 0) = -var::k(0) - 2 * (var::k(h)-var::k(0)) - 2 * h*h * var::q(0) - h * var::beta_1;
    mat(0, 1) = 2 * var::k(h) - var::k(0);
    vec[0] = -h * var::mu_1 - 4 * h*h * var::f(0);

    mat(n-1, n-1) = -var::k(1) - 2 * (var::k(1-h)-var::k(1)) - 2 * h*h * var::q(1) - h * var::beta_2;
    mat(n-1, n-2) = 2 * var::k(1-h) - var::k(1);
    vec[n-1] = -h * var::mu_2 - 4 * h*h * var::f(1);
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
        /*for (size_t i = 0; i < n; ++i) {
            std::cout << "[";
            for (size_t j = 0; j < n; ++j) {
                std::cout << mat.get(i, j) << ", ";
            }
            std::cout << "],\n";
        }*/
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
    size_t n;
    std::cin >> n;
    if (args == nullptr) {
        return 1;
    }
    return test_solve(n, argc <= 1);
}
