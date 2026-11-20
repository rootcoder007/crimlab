#include <cmath>
#include <cstdio>

#include <crimlab/crimlab.h>

int main() {
    double a = 0, b = 0;
    int rc = clab_orrisk(2.0, 0.1, 0.3, &a, &b);
    std::printf("clab_orrisk: rc %d, exposed %.7f, unexposed %.8f\n", rc, a, b);
    return (rc == 0 && std::fabs(a - 0.1475205) < 1e-6) ? 0 : 1;
}
