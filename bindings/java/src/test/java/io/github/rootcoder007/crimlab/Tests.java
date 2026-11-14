// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
package io.github.rootcoder007.crimlab;

/** The tests, with no test framework: java ... Tests; exit status 1 on a failure. */
public final class Tests {
  private Tests() {}

  private static int bad = 0;

  private static void ck(String what, boolean ok) {
    if (!ok) {
      bad++;
      System.out.println("FAIL " + what);
    }
  }

  private static String err(String f, String a) {
    try {
      Crimlab.jrun(f, a);
      return "";
    } catch (IllegalArgumentException e) {
      return e.getMessage();
    }
  }

  public static void main(String[] args) {
    String r = Crimlab.jrun("losbnd", "{\"done\": [30, 45, 60, 90, 120], \"open\": [100, 150, 200]}");
    ck("losbnd " + r, r.equals("{\"mean\":69,\"age\":150,\"lower\":99.375,\"bias\":30.375}"));
    r = Crimlab.jrun("sntbnd", "{\"y\": [1, 0, 1, 1, 0, 0], \"z\": [\"a\", \"a\", \"b\", \"b\", \"b\", \"a\"], \"treat\": \"b\"}");
    ck("sntbnd " + r, r.startsWith("{\"term\":[\"a\",\"b\",\"effect\"],\"share\":[0.5,0.5,null]"));
    ck("selftest", Crimlab.jrun("selftest", "{}").contains("\"nfail\":0,\"report\":\"\""));
    ck("jfns", Crimlab.jfns().contains("\"kinidx\":[\"g1\",\"g2\",\"freq\",\"rel?\"]"));
    ck("refused", err("orrisk", "{\"or\": 2, \"rate\": 1, \"share\": 0.3}")
        .equals("orrisk: `rate` must be one number strictly between 0 and 1."));
    ck("no function", err("nope", "{}").contains("no function nope"));
    r = Crimlab.jrun("ntdeg", "{\"from\": [\"a\", \"a\"], \"to\": [\"b\", \"c\"]}");
    ck("ntdeg " + r, r.startsWith("{\"node\":[\"a\",\"b\",\"c\"],\"deg\":[2,1,1]"));
    for (String f : new String[] {"cgstat", "bdig", "txcal", "crsev", "cconc", "clrilr", "lgtreg", "jpipe", "evlx"}) {
      ck(f, err(f, "{}").endsWith("is missing."));
    }
    StringBuilder b = new StringBuilder("{\"same\": [");
    for (int i = 1; i < 20000; i++) {
      b.append(i == 1 ? "" : ",").append(i / 7.0);
    }
    String v = b.substring(b.indexOf("[") + 1);
    r = Crimlab.jrun("lrcal", b + "], \"diff\": [" + v + "]}");
    ck("large result", r.length() > 1 << 16);
    System.out.println("Tests: " + (bad == 0 ? "all passed" : bad + " failed"));
    if (bad != 0) {
      System.exit(1);
    }
  }
}
