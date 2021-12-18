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
        double k = std::sqrt(0.4),
            C = -2 * std::exp(-0.5) / ((1+2*k) * std::exp(k) + (1-2*k) * std::exp(-k));
        return C * (std::exp(k * x) + std::exp(-k * x)) + 2 * std::exp(-0.5);
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
        mat(i, i-1) = var::k(i*h);
        mat(i, i) = -var::k(i * h) - var::k(i*h + h) - var::q(i*h) * h*h;
        mat(i, i+1) = var::k(i*h + h);
        vec[i] = -var::f(i*h) * h*h;
    }

    mat(0, 0) = -(var::k(0) + var::k(h)) / 2 - h * (var::beta_1 + h * var::q(0) / 2);
    mat(0, 1) = (var::k(0) + var::k(h)) / 2;
    mat(n-1, n-2) = (var::k(1) + var::k(1-h)) / 2;
    mat(n-1, n-1) = -(var::k(1) + var::k(1-h)) / 2 - h * (var::beta_2 + h * var::q(1) / 2);
    vec[0] = -h * (var::mu_1 - h * var::f(0) / 2);
    vec[n-1] = -h * (var::mu_2 - h * var::f(1) / 2);
    /*
    s = !echo 99 | ./main
    s = s[:-1]
    m = numpy.array(eval('['+''.join(s[:-2])+']'))
    v = numpy.array(eval('[' + s[-1] + ']'))
    a = numpy.array(eval('[' + s[-2] + ']'))
    xx = numpy.linspace(0, 1, len(v))
    plt.plot(xx, v-v[0]+g([0])[0])
    plt.plot(xx, g(xx))
    plt.show()
    */
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

    double me = 0;
    for (size_t i = 0; i < n; ++i) {
        double err = std::abs(var::model_u(i / (n-1)) - solution[i]);
        if (me < err) {
            me = err;
        }
    }
    std::cerr << "[error = " << me << "]" << std::endl;
    return 0;
}

int main() {
    std::cout.precision(16);
    size_t n;
    std::cin >> n;
    return test_solve(n, 1);
}
