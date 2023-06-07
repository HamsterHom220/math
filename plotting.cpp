#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

#include <fstream>
#include <cstdio>
#include <string>

#ifndef MATRICES_IMPORTED
#include "matrices.cpp"
#define MATRICES_IMPORTED 1
#endif

using namespace std;

void polynomial_regression_plot(int pts_num, int deg, ColumnVector& x, ColumnVector& y){
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
    VandermondeMatrix F(x, deg);
    ColumnVector* alpha = F.fit(y); // polynom coefs

    ofstream commands_file("commands.txt");
    commands_file << "set term wxt persist" << endl;
    commands_file << "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5" << endl;
    commands_file << "set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7 ps 1.5" << endl;

    // Construct plot command for polynomial and data points
    string full_cmd = "plot ";
    for (int i = deg; i >= 0; i--)
    {
        full_cmd += to_string(alpha->arr[i][0]) + "*x**" + to_string(i);
        if (i != 0)
            full_cmd += " + ";
    }
    full_cmd += " with linespoints ls 1 title 'Polynomial', '-' with points ls 2 title 'Data Points'\n";
    for (size_t i = 0; i < pts_num; i++)
        full_cmd += to_string(x.arr[i][0]) + " " + to_string(y.arr[i][0]) + "\n";
    full_cmd += "e\n";

    // Plot polynomial and data points
    commands_file << full_cmd << endl;

    commands_file.close();
    system("gnuplot commands.txt");
}

// builds a graph by connecting the points from a given set
void plot(bool parametric, ColumnVector* f1, ColumnVector* f2, ColumnVector* t=nullptr){
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
    ofstream commands_file("commands.txt");
    commands_file << "set term wxt persist" << endl;
    string cmd, data="";

    if (!parametric) {
        // 2 graphs in tOy plane: y = f1(t), y = f2(t)
        commands_file << "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5" << endl;
        commands_file << "set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 7 ps 1.5" << endl;
        cmd = "plot '-' using 1:2 with linespoints ls 1 title 'f2(t)', '-' using 1:2 with linespoints ls 2 title 'f1(t)'\n";

        // Join data for f2(t) and f1(t) into a single data block
        for (int i = 0; i < f1->get_rows(); i++) {
            data += to_string(t->arr[i][0]) + " " + to_string(f2->arr[i][0]) + "\n";
        }
        data += "e\n";
        for (int i = 0; i < f1->get_rows(); i++) {
            data += to_string(t->arr[i][0]) + " " + to_string(f1->arr[i][0]) + "\n";
        }
        data += "e\n";
    }
    else {
        // 1 graph in xOy plane: x = f1(t), y = f2(t)
        commands_file << "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5" << endl;
        cmd = "plot '-' using 1:2 with linespoints ls 1 title 'x = f1(t), y = f2(t)'\n";
        for (int i = 0; i < f1->get_rows(); i++) {
            data += to_string(f1->arr[i][0]) + " " + to_string(f2->arr[i][0]) + "\n";
        }
        data += "e\n";
    }
    // Send plot commands and data block to gnuplot
    commands_file << cmd << endl;
    commands_file << data << endl;
    commands_file.close();
    system("gnuplot commands.txt");
}
