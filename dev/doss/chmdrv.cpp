// Driver for dev/doss/toolmaRk.R: reads two profiles (one value a line)
// and prints clab_chmb's n_same n_diff U p. Build from the package root:
//   g++ -std=c++17 -O2 -Isrc src/chm.cpp src/lws.cpp dev/doss/chmdrv.cpp -o chmdrv
#include <cstdio>
#include <cstdlib>

#include "crimlab.h"

static double a[200000], b[200000], w[4000000];

static int rd(const char *f, double *v) {
  FILE *fp = std::fopen(f, "r");
  int n = 0;
  while (fp != nullptr && n < 200000 && std::fscanf(fp, "%lf", &v[n]) == 1) ++n;
  if (fp != nullptr) std::fclose(fp);
  return n;
}

int main(int argc, char **argv) {
  if (argc < 6) return 2;
  const int na = rd(argv[1], a), nb = rd(argv[2], b);
  int ns, nd;
  double u, p;
  const int rc = clab_chmb(a, na, b, nb, std::atoi(argv[3]), std::atoi(argv[4]), std::atof(argv[5]),
                           &ns, &nd, &u, &p, w);
  std::printf("%d %d %.17g %.17g rc=%d\n", ns, nd, u, p, rc);
  return 0;
}
