# Power readers. The figures over snapshots are checked by hand in
# src/t_pwr.cpp; here the readers on fake sysfs and proc trees in temporary
# folders, the parsers on recorded output, and the R entry. No test reads a
# real device except the IOReport one, on macOS only.

off <- c("nvml", "smi", "ioreport", "ioreg", "emi")
none <- c(off, "rapl", "gpu", "battery", "proc")
snap <- function(root, skip = off, pid = 0L, nvml = "", smi = "") {
  .Call(C_clab_r_pwrs, root, nvml, smi, "", skip, as.integer(pid))
}
rows <- function(s, rated = c(NA_real_, NA_real_)) {
  .Call(C_clab_r_pwrd, s, rated)
}
put <- function(root, path, x) {
  f <- file.path(root, path)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(x), f)
}
tree <- function() {
  d <- tempfile("pw")
  dir.create(d)
  d
}
row <- function(r, name, scope = "device") {
  r[r$name == name & r$scope == scope, , drop = FALSE]
}
# a snapshot as text: t, then one counter device per element of ctr
stxt <- function(t, ctr, aon = 1L, kind = "cpu") {
  paste0(
    "hwpwr1\t",
    t,
    "\t-1\t-1\t-1\n",
    paste0(
      "ior:",
      kind,
      "\t",
      kind,
      "\tapple\t",
      kind,
      "\t\t1\t",
      ctr,
      "\t0\t-1\t-1\t-1\t-1\t0\t",
      aon,
      "\t0\n",
      collapse = ""
    )
  )
}
rapl <- function(root, z, name, uj, max = 262143328850, mmio = FALSE) {
  p <- file.path("sys/class/powercap", z)
  put(root, file.path(p, "name"), name)
  if (!is.na(uj)) {
    put(root, file.path(p, "energy_uj"), uj)
  }
  put(root, file.path(p, "max_energy_range_uj"), max)
}

test_that("RAPL: packages only, DRAM as memory, wrap-around", {
  skip_on_os("windows") # a colon is not allowed in a Windows file name
  d <- tree()
  put(d, "proc/cpuinfo", "vendor_id\t: GenuineIntel")
  rapl(d, "intel-rapl:0", "package-0", 262143000000, max = 262143328850)
  rapl(d, "intel-rapl:0:0", "core", 1e6)
  rapl(d, "intel-rapl:0:1", "uncore", 1e6)
  rapl(d, "intel-rapl:0:2", "dram", 5e6)
  rapl(d, "intel-rapl:1", "psys", 9e6)
  put(d, "sys/class/powercap/intel-rapl:0/constraint_0_max_power_uw", 15e6)
  s1 <- snap(d)
  # the package wraps: 328850 uJ to the end of the range, then 671150
  rapl(d, "intel-rapl:0", "package-0", 671150)
  rapl(d, "intel-rapl:0:0", "core", 9e8)
  rapl(d, "intel-rapl:0:2", "dram", 7e6)
  rapl(d, "intel-rapl:1", "psys", 9e9)
  r <- rows(c(s1, snap(d)))
  expect_identical(sum(r$kind == "cpu"), 1L)
  expect_equal(row(r, "package-0")$joules, 1, tolerance = 1e-12)
  expect_equal(row(r, "dram of package-0")$joules, 2, tolerance = 1e-12)
  expect_identical(row(r, "package-0")$vendor, "intel")
  expect_equal(r$joules[r$kind == "total"], 3, tolerance = 1e-12)
  expect_false(any(grepl("core|psys", r$name)))
})

test_that("RAPL: an MMIO copy and per-die mirrors are counted once", {
  skip_on_os("windows")
  d <- tree()
  rapl(d, "intel-rapl:0", "package-0", 1e9)
  rapl(d, "intel-rapl-mmio:0", "package-0", 1e9)
  # two dies exposing the same counter
  rapl(d, "intel-rapl:1", "package-1", 5e9)
  rapl(d, "intel-rapl:2", "package-1-die-1", 5e9)
  s1 <- snap(d)
  rapl(d, "intel-rapl:0", "package-0", 1.004e9)
  rapl(d, "intel-rapl-mmio:0", "package-0", 1.004e9)
  rapl(d, "intel-rapl:1", "package-1", 5.010e9)
  rapl(d, "intel-rapl:2", "package-1-die-1", 5.010e9)
  r <- rows(c(s1, snap(d)))
  expect_identical(sum(r$name == "package-0"), 1L)
  expect_true(is.na(row(r, "package-1-die-1")$joules))
  expect_match(row(r, "package-1-die-1")$why, "per-die copy")
  expect_equal(r$joules[r$kind == "total"], 14, tolerance = 1e-12)
})

test_that("RAPL: unreadable, absent, malformed, and an estimate", {
  skip_on_os("windows")
  d <- tree()
  rapl(d, "intel-rapl:0", "package-0", 1e6)
  f <- file.path(d, "sys/class/powercap/intel-rapl:0/energy_uj")
  Sys.chmod(f, "000")
  on.exit(Sys.chmod(f, "644"))
  put(d, "sys/class/powercap/intel-rapl:0/constraint_0_max_power_uw", 15e6)
  put(d, "proc/stat", "cpu  100 0 100 800 0 0 0 0 0 0")
  s1 <- snap(d)
  put(d, "proc/stat", "cpu  200 0 200 1400 0 0 0 0 0 0")
  s2 <- snap(d)
  if (file.access(f, 4L) != 0L) {
    r <- rows(c(s1, s2))
    p <- row(r, "package-0")
    expect_match(p$why, "not readable: needs read access to the powercap")
    # rated 15 W x the machine's busy share 200 / 800, as an estimate
    expect_identical(p$how, "rated x utilisation")
    expect_equal(p$util, 0.25)
    expect_equal(p$watts, 3.75)
    expect_identical(c(p$flo, p$fhi), c(0.4, 2.9))
  }
  Sys.chmod(f, "644")
  writeLines("12x", f)
  expect_match(row(rows(snap(d)), "package-0")$why, "malformed")
  e <- tree()
  r <- rows(snap(e))
  expect_match(r$why[r$kind == "cpu"], "no RAPL powercap domains")
  expect_match(r$why[r$kind == "total"], "not counted: CPU")
})

test_that("GPUs: several cards of mixed vendors, with and without power", {
  skip_on_os("windows")
  d <- tree()
  card <- function(n, ven, dev, slot) {
    p <- file.path("sys/class/drm", n, "device")
    put(d, file.path(p, "vendor"), ven)
    put(d, file.path(p, "device"), dev)
    put(d, file.path(p, "uevent"), paste0("PCI_SLOT_NAME=", slot))
    p
  }
  a <- card("card0", "0x1002", "0x73bf", "0000:03:00.0")
  put(d, file.path(a, "gpu_busy_percent"), 40)
  put(d, file.path(a, "hwmon/hwmon3/name"), "amdgpu")
  put(d, file.path(a, "hwmon/hwmon3/power1_average"), 100e6)
  b <- card("card1", "0x1002", "0x6900", "0000:04:00.0")
  put(d, file.path(b, "gpu_busy_percent"), 10)
  i <- card("card2", "0x8086", "0x9b41", "0000:00:02.0")
  x <- card("card3", "0x8086", "0x56a0", "0000:05:00.0")
  put(d, file.path(x, "hwmon/hwmon5/name"), "i915")
  put(d, file.path(x, "hwmon/hwmon5/energy1_input"), 1e6)
  card("card4", "0x10de", "0x2230", "0000:06:00.0")
  put(d, "sys/class/drm/card2-eDP-1/status", "connected")
  lib <- file.path(d, "no-such-lib")
  s1 <- snap(d, nvml = lib, skip = setdiff(off, "nvml"))
  put(d, file.path(a, "hwmon/hwmon3/power1_average"), 120e6)
  put(d, file.path(x, "hwmon/hwmon5/energy1_input"), 31e6)
  Sys.sleep(0.05)
  s2 <- snap(d, nvml = lib, skip = setdiff(off, "nvml"))
  r <- rows(c(s1, s2), c(NA, 200))
  expect_identical(sum(r$kind == "gpu" & r$scope == "device"), 5L)
  g <- row(r, "card0 (1002:73bf)")
  expect_identical(g$how, "measured power integrated")
  expect_equal(g$watts, 110, tolerance = 1e-9)
  expect_equal(g$util, 0.4)
  # present without power, estimated from the rated power given
  h <- row(r, "card1 (1002:6900)")
  expect_identical(h$how, "rated x utilisation")
  expect_equal(h$watts, 20, tolerance = 1e-12)
  expect_equal(row(r, "card3 (8086:56a0)")$joules, 30, tolerance = 1e-12)
  ig <- row(r, "card2 (8086:9b41)")
  expect_true(is.na(ig$joules))
  expect_match(ig$why, "integrated")
  expect_match(row(r, "card4 (10de:2230)")$why, "NVIDIA driver library")
  tot <- r[r$kind == "total", ]
  expect_match(tot$why, "card4")
  expect_false(grepl("card2", tot$why, fixed = TRUE))
})

test_that("nvidia-smi: output, a failing program, an absent one", {
  skip_on_os("windows")
  d <- tree()
  prog <- function(body) {
    f <- tempfile("smi", tmpdir = d)
    writeLines(c("#!/bin/sh", body), f)
    Sys.chmod(f, "755")
    f
  }
  ok <- prog("printf '41.20, 7, NVIDIA A40\\n[N/A], 3, Tesla K80\\n'")
  sk <- setdiff(off, c("nvml", "smi"))
  lib <- file.path(d, "no-such-lib")
  r <- rows(snap(d, skip = sk, nvml = lib, smi = ok))
  expect_identical(r$name[r$kind == "gpu"], c("NVIDIA A40", "Tesla K80"))
  expect_equal(r$watts[r$name == "NVIDIA A40"], 41.2)
  # "gpu" leaves out NVIDIA's readers too, not only the kernel's cards
  r <- rows(snap(d, skip = c(sk, "gpu"), nvml = lib, smi = ok))
  expect_identical(sum(r$kind == "gpu"), 0L)
  r <- rows(snap(d, skip = sk, nvml = lib, smi = prog("exit 3")))
  expect_match(r$why[r$kind == "gpu"], "exit status 3")
  r <- rows(snap(d, skip = sk, nvml = lib, smi = file.path(d, "none")))
  expect_identical(sum(r$kind == "gpu"), 0L)
})

test_that("the battery counts only while it discharges", {
  skip_on_os("windows")
  d <- tree()
  b <- "sys/class/power_supply/BAT0"
  put(d, file.path(b, "type"), "Battery")
  put(d, file.path(b, "status"), "Not charging")
  put(d, file.path(b, "power_now"), 0)
  put(d, "sys/class/power_supply/AC/type", "Mains")
  r <- rows(snap(d))
  expect_match(row(r, "battery BAT0", "machine")$why, "while it discharges")
  put(d, file.path(b, "status"), "Discharging")
  put(d, file.path(b, "power_now"), 9.5e6)
  expect_equal(row(rows(snap(d)), "battery BAT0", "machine")$watts, 9.5)
  unlink(file.path(d, b, "power_now"))
  put(d, file.path(b, "current_now"), 1e6)
  put(d, file.path(b, "voltage_now"), 12e6)
  r <- rows(snap(d))
  expect_equal(row(r, "battery BAT0", "machine")$watts, 12, tolerance = 1e-12)
  # the machine is beside the devices, not in their total
  expect_true(is.na(r$watts[r$kind == "total"]))
})

test_that("the process tree's share of the package", {
  skip_on_os("windows")
  d <- tree()
  st <- function(pid, comm, ppid, ut, sy, cu = 0, cs = 0) {
    put(
      d,
      file.path("proc", pid, "stat"),
      sprintf(
        "%d (%s) S %d 1 1 0 -1 0 0 0 0 0 %d %d %d %d 20 0 1 0",
        pid,
        comm,
        ppid,
        ut,
        sy,
        cu,
        cs
      )
    )
  }
  rapl(d, "intel-rapl:0", "package-0", 1e6)
  put(d, "proc/stat", "cpu  100 0 100 800 0 0 0 0 0 0")
  st(100, "R", 1, 10, 0)
  st(101, "a (b) c", 100, 5, 0)
  st(200, "other", 1, 50, 0)
  s1 <- snap(d, pid = 100L)
  rapl(d, "intel-rapl:0", "package-0", 11e6)
  put(d, "proc/stat", "cpu  150 0 150 900 0 0 0 0 0 0")
  st(100, "R", 1, 20, 5, 5, 0)
  st(101, "a (b) c", 100, 15, 0)
  st(200, "other", 1, 90, 0)
  r <- rows(c(s1, snap(d, pid = 100L)))
  # the tree used 30 of the 100 busy ticks
  p <- row(r, "package-0", "process")
  expect_equal(p$joules, 3, tolerance = 1e-12)
  expect_equal(p$util, 0.3)
  expect_equal(row(r, "package-0")$joules, 10, tolerance = 1e-12)
})

test_that("the NVIDIA library: absent is no reading, not an error", {
  d <- tree()
  r <- rows(snap(d, skip = setdiff(off, "nvml"), nvml = file.path(d, "x.so")))
  expect_identical(sum(r$kind == "gpu"), 0L)
})

test_that("Apple counters: steady, paused, frozen, caught up", {
  sq <- function(ctr) {
    rows(vapply(seq_along(ctr), function(k) stxt(0.25 * k, ctr[k]), ""))
  }
  expect_equal(sq(0:4)$joules[1], 4, tolerance = 1e-12)
  r <- sq(c(0, 1, 1, 3, 4))
  expect_equal(r$joules[1], 4, tolerance = 1e-12)
  expect_match(r$why[1], "spread over the pause")
  r <- sq(c(0, 1, 2, 3, 3))
  expect_true(is.na(r$joules[1]))
  expect_match(r$why[1], "stopped changing")
  expect_true(is.na(sq(c(0, 30, 31, 32, 33))$joules[1]))
  expect_true(is.na(sq(c(5, 5))$joules[1]))
  # a counter that may rest at zero (the neural engine) is not a fault
  z <- vapply(1:3, stxt, "", ctr = 0, aon = 0L, kind = "ane")
  expect_equal(rows(z)$joules[1], 0, tolerance = 1e-12)
})

test_that("snapshot text that is not a reading is refused", {
  errs(
    "hwpwr",
    list(
      list("`since` must be an earlier result", since = data.frame()),
      list("`dt` must be one number from 0", dt = -1),
      list("`step` must be one number above 0", step = 0),
      list("`skip` must name sources", skip = "x"),
      list("`paths` must be named strings", paths = "a"),
      list("`rated` must be watts", rated = c(cpu = -1)),
      list("`pid` must be a process id", pid = 1.5)
    )
  )
  bad <- data.frame()
  attr(bad, "snap") <- "hwpwr1\t1\n"
  expect_error(hwpwr(since = bad, skip = none), "not a reading of hwpwr")
  late <- data.frame()
  attr(late, "snap") <- stxt(1e15, 1)
  expect_error(hwpwr(since = late, skip = none), "clock went back")
})

test_that("parsers: nvidia-smi, ioreg and the EMI metadata", {
  txt <- "41.20, 7, NVIDIA A40\n[N/A], [N/A], K80\n"
  p <- .Call(C_clab_r_pwrp, "smi", txt, NULL)
  expect_equal(p$watts, c(41.2, NA))
  expect_equal(p$util, c(0.07, NA))
  expect_match(p$why[2], "no power")
  expect_identical(nrow(.Call(C_clab_r_pwrp, "smi", "", NULL)), 0L)
  expect_match(.Call(C_clab_r_pwrp, "smi", "garbage\n", NULL)$why, "not power")
  acc <- paste0(
    "+-o AGXAcceleratorG14G  <class AGXAcceleratorG14G>\n",
    "    \"PerformanceStatistics\" = {\"Tiler Utilization %\"=4,",
    "\"Device Utilization %\"=14,\"Renderer Utilization %\"=4}\n"
  )
  expect_equal(.Call(C_clab_r_pwrp, "ioacc", acc, NULL), 0.14)
  bt <- function(ext, amp) {
    paste0(
      "  \"ExternalConnected\" = ",
      ext,
      "\n",
      "  \"BatteryData\" = {\"CellVoltage\"=(4100,4100),\"Voltage\"=1}\n",
      "  \"InstantAmperage\" = ",
      amp,
      "\n",
      "  \"Voltage\" = 12302\n"
    )
  }
  b <- .Call(C_clab_r_pwrp, "iobat", bt("No", "18446744073709550616"), NULL)
  expect_equal(b$watts, 12.302)
  b <- .Call(C_clab_r_pwrp, "iobat", bt("Yes", "0"), NULL)
  expect_true(is.na(b$watts))
  expect_match(b$why, "mains")
  expect_match(.Call(C_clab_r_pwrp, "iobat", "", NULL)$why, "no battery")
  expect_match(
    .Call(C_clab_r_pwrp, "iobat", bt("No", "-12x"), NULL)$why,
    "malformed"
  )
  w16 <- function(s) {
    c(iconv(s, to = "UTF-16LE", toRaw = TRUE)[[1]], as.raw(c(0, 0)))
  }
  ch <- function(s, unit = 0L) {
    n <- w16(s)
    c(as.raw(c(unit, 0, 0, 0, length(n), 0)), n)
  }
  v2 <- c(raw(66), as.raw(c(2, 0)), ch("RAPL_Package0_PKG"), ch("Other", 1L))
  m <- .Call(C_clab_r_pwrp, "emi", v2, 2L)
  expect_identical(m$name, c("RAPL_Package0_PKG", "Other"))
  expect_identical(m$pwh, c(TRUE, FALSE))
  expect_length(.Call(C_clab_r_pwrp, "emi", v2[1:80], 2L)$name, 0L)
  v1 <- c(raw(70), as.raw(c(length(w16("CPU")), 0)), w16("CPU"))
  expect_identical(.Call(C_clab_r_pwrp, "emi", v1, 1L)$name, "CPU")
})

test_that("hwpwr(): the shape, and an interval since a reading", {
  a <- hwpwr(0, skip = none)
  expect_named(
    a,
    c(
      "kind",
      "vendor",
      "name",
      "scope",
      "watts",
      "joules",
      "util",
      "how",
      "flo",
      "fhi",
      "why"
    )
  )
  expect_identical(a$kind, "total")
  b <- hwpwr(since = a, skip = none)
  expect_gte(attr(b, "dt"), 0)
  expect_type(attr(b, "snap"), "character")
})

test_that("IOReport on this Mac: CPU and GPU energy, or the reason", {
  skip_on_cran()
  skip_if_not(identical(Sys.info()[["sysname"]], "Darwin"), "macOS only")
  r <- hwpwr(0.5, skip = c("battery", "ioreg"))
  cpu <- r[r$kind == "cpu" & r$scope == "device", ]
  expect_identical(nrow(cpu), 1L)
  expect_true(cpu$joules > 0 || nzchar(cpu$why))
})

test_that("rated power: the long-term limit, the one given, through hwpwr()", {
  skip_on_os("windows")
  d <- tree()
  rapl(d, "intel-rapl:0", "package-0", NA)
  z <- "sys/class/powercap/intel-rapl:0"
  put(d, file.path(z, "constraint_0_name"), "long_term")
  put(d, file.path(z, "constraint_0_power_limit_uw"), 10e6)
  put(d, "proc/stat", "cpu  100 0 100 800 0 0 0 0 0 0")
  s1 <- snap(d)
  Sys.sleep(0.02)
  put(d, "proc/stat", "cpu  200 0 200 1400 0 0 0 0 0 0")
  s2 <- snap(d)
  # the zone has no energy counter: the limit 10 W x the busy share 0.25
  p <- row(rows(c(s1, s2)), "package-0")
  expect_identical(p$how, "rated x utilisation")
  expect_equal(p$watts, 2.5, tolerance = 1e-12)
  expect_match(p$why, "the rated power the machine reports \\(10 W\\)")
  g <- row(rows(c(s1, s2), c(20, NA)), "package-0")
  expect_equal(g$watts, 5, tolerance = 1e-12)
  expect_match(g$why, "the rated power given \\(20 W\\)")
  # the same through hwpwr(), the tree given by paths: busy 200 of 1000
  a <- hwpwr(0, skip = off, paths = c(root = d))
  Sys.sleep(0.02)
  put(d, "proc/stat", "cpu  300 0 300 2200 0 0 0 0 0 0")
  b <- hwpwr(since = a, skip = off, paths = c(root = d), rated = c(cpu = 20))
  h <- row(b, "package-0")
  expect_equal(h$util, 0.2, tolerance = 1e-12)
  expect_equal(h$watts, 4, tolerance = 1e-12)
  expect_match(h$why, "the rated power given \\(20 W\\)")
})

test_that("hwpwr() reads every step over dt", {
  r <- hwpwr(0.05, step = 0.01, skip = none)
  expect_gt(attr(r, "dt"), 0)
  expect_identical(r$kind, "total")
})

test_that("GPUs and a battery that give no number", {
  skip_on_os("windows")
  d <- tree()
  card <- function(n, ven, dev, slot) {
    p <- file.path("sys/class/drm", n, "device")
    put(d, file.path(p, "vendor"), ven)
    put(d, file.path(p, "device"), dev)
    put(d, file.path(p, "uevent"), paste0("PCI_SLOT_NAME=", slot))
  }
  card("card0", "0x8086", "0x56a0", "0000:05:00.0")
  card("card1", "0x1a03", "0x2000", "0000:07:00.0")
  b <- "sys/class/power_supply/BAT0"
  put(d, file.path(b, "type"), "Battery")
  put(d, file.path(b, "status"), "Discharging")
  r <- rows(snap(d))
  # a discrete Intel card without hwmon energy, a vendor with no reader
  expect_match(row(r, "card0 (8086:56a0)")$why,
               "the kernel exposes no energy file for it")
  expect_match(row(r, "card1 (1a03:2000)")$why, "no reader for this vendor")
  expect_match(row(r, "battery BAT0", "machine")$why,
               "the battery reports no discharge power")
})

test_that("the NVIDIA library: one without the NVML functions", {
  # the package's own library loads but exports no nvml symbol
  lib <- normalizePath(getLoadedDLLs()[["crimlab"]][["path"]],
                       winslash = "\\", mustWork = FALSE)
  sk <- setdiff(off, "nvml")
  r <- rows(snap(tree(), skip = sk, nvml = lib))
  expect_match(r$why[r$kind == "gpu"],
               "the NVIDIA driver library lacks the functions to enumerate")
  # the default names: no library here, or the GPUs of a real one
  r0 <- rows(snap(tree(), skip = sk, nvml = ""))
  expect_true(all(r0$vendor[r0$kind == "gpu"] == "nvidia"))
})

test_that("parsers: a charging battery and an unknown parser", {
  bt <- paste0(
    "  \"ExternalConnected\" = No\n",
    "  \"InstantAmperage\" = 1500\n",
    "  \"Voltage\" = 12302\n"
  )
  b <- .Call(C_clab_r_pwrp, "iobat", bt, NULL)
  expect_true(is.na(b$watts))
  expect_identical(b$why, "the battery is charging")
  expect_error(.Call(C_clab_r_pwrp, "nope", "", NULL), "unknown parser")
})

test_that("snapshot rows: other kinds, process GPU shares, gaps", {
  # one device line of the snapshot text: id, kind, vendor, name, why,
  # src, ctr, wrap, w, util, rated, pown, inpkg, aon, bad
  dv <- function(id, kind, src, ctr = 0, w = -1, pown = -1, bad = 0) {
    paste(id, kind, "v", id, "", src, ctr, 0, w, -1, -1, pown, 0, 0, bad,
          sep = "\t")
  }
  sn <- function(t, ...) {
    paste0("hwpwr1\t", t, "\t-1\t-1\t-1\n", paste0(c(...), "\n",
                                                   collapse = ""))
  }
  two <- function(t) {
    sn(t, dv("zz", "zz", 2, w = 5), dv("g", "gpu", 2, w = 10, pown = 0.5),
       dv("c", "cpu", 2, w = 1))
  }
  r <- rows(c(two(1), two(3)))
  # known kinds first, in their order; an unknown kind after them
  expect_identical(r$kind, c("cpu", "gpu", "zz", "gpu", "total"))
  p <- row(r, "g", "process")
  expect_equal(p$joules, 10, tolerance = 1e-12)
  expect_equal(p$util, 0.5)
  expect_identical(
    p$why,
    "the process's share of the GPU's compute samples, 50.0%"
  )
  expect_equal(r$joules[r$kind == "total"], 22, tolerance = 1e-12)
  one <- rows(sn(1, dv("c", "cpu", 1, ctr = 5)))
  expect_match(one$why[1], "a counter needs an interval: dt > 0 or since")
  gap <- rows(c(sn(1, dv("c", "cpu", 1, ctr = 1)),
                sn(2, dv("c", "cpu", 1, ctr = 2, bad = 1)),
                sn(3, dv("c", "cpu", 1, ctr = 3))))
  expect_true(is.na(gap$joules[1]))
  expect_match(gap$why[1], "not readable at every snapshot")
  pw <- rows(c(sn(1, dv("g", "gpu", 2, w = 5)), sn(2, dv("g", "gpu", 2)),
               sn(3, dv("g", "gpu", 2, w = 5))))
  expect_true(is.na(pw$joules[1]))
  expect_match(pw$why[1], "power was not read at every snapshot")
})
