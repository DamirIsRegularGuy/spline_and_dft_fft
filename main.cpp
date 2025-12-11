#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <chrono>

//Функция для ДПФ и БПФ
double func_dft(double x){
    return x * (1-x) * std::log1p(-x);
}

//Реализация прямого и обратного ДПФ в одной функции, прямое и обратное преобразование определяются через булевую переменную dir
std::vector<std::complex<double>> dft(std::vector<std::complex<double>>& f, bool dir = true){
    std::vector<std::complex<double>> A_m(f.size());
    std::complex<double> i(0,1);
    int sign = dir? -1: 1;
    for(int m=0; m<f.size(); m++){
        std::complex<double> sum(0,0);
        for(int j=0; j<f.size(); j++){
            sum += f[j] * std::exp(sign * 2.0 * M_PI * m * i * double(j)/double(f.size()));
        }
        if(dir) sum /= double(f.size());
        A_m[m]=sum;
    }

    return A_m;
}

//Функция по вычислению погрешности через норму супремом
double sup_norm(std::vector<std::complex<double>>& f, std::vector<std::complex<double>>& f_j){
    double err=0.0;
    for(int i=0; i<f.size(); i++){
        if(std::abs(f[i].real() - f_j[i].real())>err) err=std::abs(f[i].real() - f_j[i].real());
    }
    return err;
}

//Функция по сохранению ДПФ и БПФ для дальнейшей визуализации
void save_to_file(std::string file, std::vector<std::complex<double>>& analitic, std::vector<std::complex<double>>& f_j, int N){
    std::ofstream filename(file);
    filename << "x analitic dft/fft\n";
    for(int i=0; i<N; i++){
        filename << std::scientific << std::setprecision(15) << double(i)/double(N) << " " << analitic[i].real() << " " << f_j[i].real() << "\n";
    }
    filename.close();
}

//Метод прогонки для вычисления вторых производных для общей формулы кубического сплайна
//так как сетка равномерная, то нет необходимости в трех векторах, а просто в постоянных a, b, c
std::vector<double> progon_method(double h, const std::vector<double>& matrix_val) {
    int n = matrix_val.size() + 1;
    double a = h/6;
    double b = 2*h/3;
    double c = h/6;

    std::vector<double> cp(n-1);
    std::vector<double> dp(n-1);

    cp[0] = c / b;
    dp[0] = matrix_val[0] / b;

    for(int i = 1; i < n-1; i++){
        double denom = b - a * cp[i-1];
        cp[i] = c / denom;
        dp[i] = (matrix_val[i] - a * dp[i-1]) / denom;
    }

    std::vector<double> M(n+1, 0.0);
    M[n-1] = dp[n-2];

    for(int i = n-3; i >= 0; i--){
        M[i+1] = dp[i] - cp[i] * M[i+2];
    }

    return M;
}

//Функция по сохранению данных для визуализации сплайна
void save_spline_to_file(const std::vector<double>& X,
                         const std::vector<double>& S,
                         const std::vector<double>& Sin)
{
    std::ofstream file("spline_graphics.dat");
    file << "x spline sin\n";

    for(size_t i = 0; i < X.size(); i++){
        file << std::scientific << std::setprecision(15) << X[i] << " " << S[i] << " " << Sin[i] << "\n";
    }

    file.close();
}

//Функия БПФ так же в одной функции реализовано прямое и обратное преобразование через invert
void fft(std::vector<std::complex<double>>& a, bool invert = true) {
    int n = a.size();
    for (int i=1, j=0; i<n; i++){
        int bit = n >> 1;
        for (; j>=bit; bit>>=1) j -= bit;
        j += bit;
        if(i<j) swap(a[i], a[j]);
    }

    for (int len=2; len<=n; len<<=1){
        double ang = 2*M_PI/len * (invert ? -1 : 1);
        std::complex<double> wlen(cos(ang), sin(ang));
        for (int i=0; i<n; i+=len){
            std::complex<double> w(1,0);
            for (int j=0; j<len/2; j++){
                std::complex<double> u = a[i+j];
                std::complex<double> v = a[i+j+len/2]*w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert)
        for (auto &x : a) x /= n;
}


int main() {
    // ДПФ
    int N = 10000;
    //Все вектора сделаны комплексными во избежание ошибок и общей логики реализции функций
    std::vector<std::complex<double>> func_f(N);
    std::vector<std::complex<double>> A_m(N);
    std::vector<std::complex<double>> f_j(N);

    //Определение значений функции на равномерной сетке
    for(int i=0; i<N; i++){
        func_f[i]= std::complex<double>(func_dft(double(i)/double(N)), 0.0);
    }

    auto start_dft = std::chrono::high_resolution_clock::now();
    //прямое ДПФ
    A_m = dft(func_f, true);
    //Обратное ДПФ
    f_j = dft(A_m, false);
    auto end_dft = std::chrono::high_resolution_clock::now();
    auto elapsed_dft = std::chrono::duration<double, std::milli>(end_dft - start_dft);
    std::cout << "Time for dft " << elapsed_dft.count() << " milliseconds\n";


    //Вывод в консоль значений для самой функции на узлах, прямого и обратного ДПФ
//    std::cout << std::fixed << std::setprecision(6);
//    std::cout << std::setw(13) << " f(x)"
//              << std::setw(13) << " A_m"
//              << std::setw(13) << " f_j"
//              << "\n";
//    for(int i = 0; i < N; i++){
//        std::cout << std::setw(13) << " " << func_f[i].real()
//                  << std::setw(13) << " " << std::complex<double>(A_m[i])
//                  << std::setw(13) << " " << f_j[i].real()
//                  << "\n";
//    }

    //вычисление ошибки и сохранение значений в файл
    std::cout << "\n";
    std::cout << "dft error: ";
    std::cout << std::scientific << std::setprecision(15) << sup_norm(func_f, f_j) << "\n";
    save_to_file("dft_graphics.dat",func_f,f_j,N);


    //Кубический сплайн
    int n=10000;
    double x0=0, xn=1;
    double h=(xn-x0)/n;
    std::vector<double> arg_x(n);
    std::vector<double> func_val(n);
    std::vector<double> matrix_val(n);

    std::vector<double> a(n);
    std::vector<double> b(n);
    std::vector<double> c(n);
    std::vector<double> d(n);

    //Определение узлов и значений функции в этих узнах
    for(int i=0; i<n; i++){
        arg_x[i]=x0 + i*h;
        func_val[i] = std::sin(arg_x[i]);
    }

    //Вычисляем правую часть от равенства для 3-хдиагональной матрицы
    for(int i=1; i<n; i++){
        matrix_val[i-1] = (std::sin(arg_x[i+1]) - std::sin(arg_x[i]))/h - (std::sin(arg_x[i]-std::sin(arg_x[i-1])))/h;
    }

    //Решение 3-хдиагональной матрицы
    std::vector<double> M = progon_method(h, matrix_val);

//    for(int i = 0; i < n; i++){
//        std::cout << "x=" << arg_x[i] << " M=" << M[i] << "\n";
//    }

    //Вычисление коэффицентов для общей формулы кубического сплайна
    for(int i=0; i<n; i++){
        a[i] = func_val[i];
        b[i] = (func_val[i+1] - func_val[i])/h - (2*M[i]+M[i+1])*h/6;
        c[i] = M[i]/2;
        d[i] = (M[i+1]-M[i])/(6*h);
    }

    std::vector<double> plot_x;
    std::vector<double> plot_spline;
    std::vector<double> plot_sin;

    //Деление каждого отрезка для создания более гладкого сплайна, а не только в заданных узлах
    //Выыод значений и запись их в файл для дальнейшей визуализации
    int pointsPerInterval = 10;
//    std::cout << std::setw(13) << "x"
//              << std::setw(13) << "Sx"
//              << std::setw(13) << "sin(x)"
//              << "\n";
    for(int i = 0; i < N-1; i++){
        for(int j = 0; j <= pointsPerInterval; j++){
            double t = j * h / pointsPerInterval;
            double Sx = a[i] + b[i]*t + c[i]*t*t + d[i]*t*t*t;
            double xx = arg_x[i] + t;

            plot_x.push_back(xx);
            plot_spline.push_back(Sx);
            plot_sin.push_back(std::sin(xx));

//            std::cout << std::setw(12) << " " << xx
//                      << std::setw(12) << " " << Sx
//                      << std::setw(12) << " " << std::sin(xx)
//                      << "\n";
        }
    }

    save_spline_to_file(plot_x, plot_spline, plot_sin);



    //БПФ
    int n_fft = 10000;
    int n2_fft = 1;

    //Изменение размера числа узлов, так как БПФ реализовано на побайтовом смещении, соответсвенно число узлов, должно быть какой-нибудь степенью двойки
    if(!((n_fft>0)&&((n_fft&(n_fft-1))==0))){
        while(n2_fft<n_fft){
            n2_fft <<= 1;
        }
    } else{
        n2_fft=n_fft;
    }

    //Заполнение дополнительных узлов 0
    std::vector<std::complex<double>> func_f_fft(n2_fft, 0.0);

    //Вычисление значений функции в узлах
    for(int i = 0; i < n_fft; i++){
        func_f_fft[i] = std::complex<double>(func_dft(double(i)/double(n_fft)), 0.0);
    }

    auto start_fft = std::chrono::high_resolution_clock::now();
    //Прямое БПФ
    std::vector<std::complex<double>> A_m_fft = func_f_fft;
    fft(A_m_fft, true);

    //Обратное БПФ
    std::vector<std::complex<double>> f_j_fft = A_m_fft;
    fft(f_j_fft, false);

    auto end_fft = std::chrono::high_resolution_clock::now();
    auto elapsed_fft = std::chrono::duration<double, std::milli>(end_fft - start_fft);
    std::cout << "Time for fft " << elapsed_fft.count() << " milliseconds\n";



    //Выыод прямого, обратного и аналатических значений БПФ
//    std::cout << std::fixed << std::setprecision(6);
//    std::cout << std::setw(13) << " f(x)"
//              << std::setw(13) << " A_m_fft"
//              << std::setw(13) << " f_j_fft"
//              << "\n";
//    for(int i = 0; i < n_fft; i++){
//        std::cout << std::setw(13) << " " << func_f_fft[i].real()
//                  << std::setw(13) << " " << std::complex<double>(A_m_fft[i])
//                  << std::setw(13) << " " << f_j_fft[i].real()
//                  << "\n";
//    }

    //Вычисление ошибки и сохранение в файл значений для дальнейшей визуализации
    std::cout << "\n";
    std::cout << "fft error: ";
    std::cout << std::scientific << std::setprecision(15) << sup_norm(func_f_fft, f_j_fft) << "\n";
    save_to_file("fft_graphics.dat",func_f_fft,f_j_fft,n_fft);





}

/*
//Алгоритм реализации кубического сплайна
//Создать равномерную сетку
h = (xN - x0)/N
for i-N:
    x[i] = x0 + i*h
    f[i] = f(x[i])

// Вычисляем вектор RHS для трёхдиагональной системы
for i-N-1:
    RHS[i] = (f[i+1]-f[i])/h - (f[i]-f[i-1])/h

// Создаём 3-хдиагональную матрицу
for i-N-1:
    A[i][i-1] = h/6
    A[i][i]   = 2*h/3
    A[i][i+1] = h/6

// Решаем 3-хдиагональную систему A * M = RHS для M[1..N-1]

// Вычисляем коэффициенты S_i
for i-N-1:
    a[i] = f[i]
    b[i] = (f[i+1]-f[i])/h - (2*M[i]+M[i+1])*h/6
    c[i] = M[i]/2
    d[i] = (M[i+1]-M[i])/(6*h)

// Теперь можно вычислять S(x) для любого x на собственном отрезке:
// S(x) = a[i] + b[i]*(x-x_i) + c[i]*(x-x_i)^2 + d[i]*(x-x_i)^3
*/