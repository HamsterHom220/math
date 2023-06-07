#ifndef MATRICES_IMPORTED
#include "matrices.cpp"
#define MATRICES_IMPORTED 1
#endif

using namespace std;

class PreyPredatorModel{
private:
    int v0, k0, n; // initial number of preys, initial number of predators, number of time moments
    double t_max, a1, b1, a2, b2; // max time moment, positive constants:
    //                                              prey birthrate,
    //                                              prey deathrate with respect to the number of predators,
    //                                              predator deathrate,
    //                                              predator birthrate with respect to the number of preys
public:
    PreyPredatorModel(int v_init,int k_init,int moments_num,double max_moment,double v_birth,double v_death,double k_death,double k_birth): 
    v0(v_init), k0(k_init), n(moments_num), t_max(max_moment), a1(v_birth), b1(v_death), a2(k_death), b2(k_birth),
    v(n+1),k(n+1),t(n+1){};
    
    ColumnVector v;
    ColumnVector k;
    ColumnVector t;
    
    void build() {
        vector<double> v_temp(n+1);
        vector<double> k_temp(n+1);
        double incr = t_max / n;
        v0 -= a2 / b2;
        k0 -= a1 / b1;
        for (int i = 0; i <= n; i++) {
            t_max = incr * i;
            t.arr[i][0] = t_max;
            v_temp[i] = a2 / b2 + v0 * cos(t_max * sqrt(a1 * a2)) - k0 * (sqrt(a2) * b1 / b2 / sqrt(a1)) * sin(t_max * sqrt(a1 * a2));
            k_temp[i] = a1 / b1 + v0 * (sqrt(a1) * b2 / b1 / sqrt(a2)) * sin(t_max * sqrt(a1 * a2)) + k0 * cos(t_max * sqrt(a1 * a2));
        }

        for (int i = 0; i <= n; i++) {
            v.arr[i][0] = v_temp[i];
            k.arr[i][0] = k_temp[i];
        }
    }
};


