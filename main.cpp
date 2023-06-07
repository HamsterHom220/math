#include "plotting.cpp"
#include "prey_predator.cpp"
using namespace std;

// prey-predator plotting test
int main() {
    int v0, k0, n; // initial number of preys, initial number of predators, number of time moments
    double t_max, a1, b1, a2, b2; // max time moment, positive constants:
    //                                              prey birthrate,
    //                                              prey deathrate with respect to the number of predators,
    //                                              predator deathrate,
    //                                              predator birthrate with respect to the number of preys

    v0 = 300;
    k0 = 55;
    n = 500;
    t_max = 850;
    a1 = 0.5;
    b1 = 0.35;
    a2 = 0.48;
    b2 = 0.3;

    PreyPredatorModel model(v0,k0,n,t_max,a1,b1,a2,b2);
    model.build();
    plot(false,&(model.v),&(model.k),&(model.t));
    return 0;
}
