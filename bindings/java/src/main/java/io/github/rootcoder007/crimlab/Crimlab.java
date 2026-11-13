// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
package io.github.rootcoder007.crimlab;

import static java.lang.foreign.ValueLayout.ADDRESS;
import static java.lang.foreign.ValueLayout.JAVA_INT;

import java.io.IOException;
import java.io.InputStream;
import java.lang.foreign.Arena;
import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.Linker;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.SymbolLookup;
import java.lang.invoke.MethodHandle;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Criminology and sociolegal statistics: every function of the R and Python packages crimlab, by
 * its R name, through the Foreign Function and Memory API (JDK 22 or later). A function takes its
 * named arguments as a JSON object and returns its result as JSON: a data frame is an object of
 * columns, a named vector an object of numbers, a missing value {@code null}, and a number JSON
 * cannot hold {@code "NaN"}, {@code "Inf"} or {@code "-Inf"}.
 *
 * <pre>{@code
 * String r = Crimlab.jrun("orrisk", "{\"or\": 2, \"rate\": 0.1, \"share\": 0.3}");
 * }</pre>
 *
 * <p>The shared library is the file named by the system property {@code crimlab.library} or the
 * environment variable {@code CRIMLAB_LIBRARY}; else the copy the jar carries for this operating
 * system and CPU ({@code native/<os>-<arch>/}), unpacked to a private temporary folder; else
 * {@code libcrimlab} on the system's library path.
 */
@SuppressWarnings("restricted") // the FFM calls this class exists to make
public final class Crimlab {
  private Crimlab() {}

  private static final Linker L = Linker.nativeLinker();
  private static final SymbolLookup LIB = lookup();
  private static final MethodHandle JRUN =
      fn("clab_jrun", FunctionDescriptor.of(JAVA_INT, ADDRESS, ADDRESS, ADDRESS, JAVA_INT, ADDRESS));
  private static final MethodHandle JFNS = fn("clab_jfns", FunctionDescriptor.of(ADDRESS));

  /**
   * Runs a function.
   *
   * @param name the R name of the function, such as {@code "orrisk"}
   * @param args its named arguments, a JSON object
   * @return the result, JSON
   * @throws IllegalArgumentException when the function refuses its arguments; the message names
   *     the function and what is wrong
   */
  public static String jrun(String name, String args) {
    try (Arena a = Arena.ofConfined()) {
      MemorySegment f = a.allocateFrom(name), x = a.allocateFrom(args);
      MemorySegment need = a.allocate(JAVA_INT);
      int cap = 1 << 16;
      for (; ; ) {
        MemorySegment out = a.allocate(cap);
        int rc = (int) JRUN.invokeExact(f, x, out, cap, need);
        int n = need.get(JAVA_INT, 0);
        if (rc == 3 && n > cap) {
          cap = n;
          continue;
        }
        String s = out.getString(0);
        if (rc != 0) {
          throw new IllegalArgumentException(s);
        }
        return s;
      }
    } catch (RuntimeException | Error e) {
      throw e;
    } catch (Throwable e) {
      throw new IllegalStateException(e);
    }
  }

  /**
   * The functions and their arguments.
   *
   * @return a JSON object from function name to argument names ({@code ?} marks an optional one)
   */
  public static String jfns() {
    try {
      MemorySegment s = (MemorySegment) JFNS.invokeExact();
      return s.reinterpret(Long.MAX_VALUE).getString(0);
    } catch (Throwable e) {
      throw new IllegalStateException(e);
    }
  }

  private static MethodHandle fn(String name, FunctionDescriptor d) {
    MemorySegment s =
        LIB.find(name).orElseThrow(() -> new UnsatisfiedLinkError("crimlab: no symbol " + name));
    return L.downcallHandle(s, d);
  }

  private static SymbolLookup lookup() {
    List<String> tried = new ArrayList<>();
    String p = System.getProperty("crimlab.library");
    if (p == null || p.isEmpty()) {
      p = System.getenv("CRIMLAB_LIBRARY");
    }
    if (p != null && !p.isEmpty()) {
      return SymbolLookup.libraryLookup(Path.of(p).toAbsolutePath(), Arena.global());
    }
    String file = file();
    String pl = platform();
    if (pl != null) {
      String res = "/native/" + pl + "/" + file;
      try (InputStream in = Crimlab.class.getResourceAsStream(res)) {
        if (in != null) {
          // a fresh folder only this user can open; removed at exit
          Path d = Files.createTempDirectory("crimlab");
          Path f = d.resolve(file);
          Files.copy(in, f, StandardCopyOption.REPLACE_EXISTING);
          d.toFile().deleteOnExit();
          f.toFile().deleteOnExit();
          return SymbolLookup.libraryLookup(f, Arena.global());
        }
      } catch (IOException | IllegalArgumentException e) {
        tried.add("jar:" + res);
      }
    }
    for (String n : new String[] {file, "libcrimlab.so.0", "libcrimlab.0.dylib"}) {
      try {
        return SymbolLookup.libraryLookup(n, Arena.global());
      } catch (IllegalArgumentException e) {
        tried.add(n);
      }
    }
    throw new UnsatisfiedLinkError(
        "crimlab: cannot load the library (tried "
            + String.join(", ", tried)
            + "); set -Dcrimlab.library=<path> or CRIMLAB_LIBRARY");
  }

  /** The library's file name in the jar. */
  static String file() {
    return System.mapLibraryName("crimlab");
  }

  /** The platform's name in the jar: linux, macos or windows, then x86_64 or arm64. */
  static String platform() {
    String os = System.getProperty("os.name", "").toLowerCase(Locale.ROOT);
    String arch = System.getProperty("os.arch", "").toLowerCase(Locale.ROOT);
    String o =
        os.startsWith("linux")
            ? "linux"
            : os.startsWith("mac") ? "macos" : os.startsWith("windows") ? "windows" : null;
    String a =
        arch.equals("amd64") || arch.equals("x86_64")
            ? "x86_64"
            : arch.equals("aarch64") || arch.equals("arm64") ? "arm64" : null;
    return o == null || a == null ? null : o + "-" + a;
  }
}
