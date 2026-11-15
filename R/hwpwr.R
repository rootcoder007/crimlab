#' Power and energy of the machine's devices
#'
#' Reads what the CPU packages, the memory and the GPUs of this machine draw,
#' from the operating system's own energy counters and power readings, over
#' an interval: `dt` seconds from now, or since an earlier reading
#' (`since`). Every figure says how it was obtained, and a device that
#' cannot be measured is listed with the reason, never counted as zero. It
#' is the measured layer under an energy or carbon footprint of a
#' computation.
#'
#' Each row is one device (`scope = "device"`), the share of a device that
#' belongs to the process `pid` (`scope = "process"`), the whole machine
#' (`scope = "machine"`, the battery) or the total over the devices
#' (`kind = "total"`). `how` is one of
#' * `"measured counter"`: the difference of an energy counter;
#' * `"measured power integrated"`: power read at every sample, integrated
#'   by the trapezoid rule;
#' * `"rated x utilisation"`: an estimate, a rated power times the measured
#'   utilisation, with the factors `flo` and `fhi` between which the true
#'   figure may lie (0.4 and 2.9: the power-from-load model of CodeCarbon
#'   was off by -60% to +190% against RAPL in its own comparison);
#' * `"unknown"`, with the reason in `why`.
#'
#' The sources, by system:
#' * Linux, CPU and memory: the RAPL counters of the powercap tree
#'   (`/sys/class/powercap`, `energy_uj` in microjoules). Only the top-level
#'   `package` domains are counted for the CPU and the `dram` domain for
#'   memory; `core`, `uncore` and `psys` are never added (they are parts of
#'   a package or more than the CPU). A counter that passes its range
#'   (`max_energy_range_uj`) is unwrapped, a package seen through both the
#'   MSR and the MMIO interface is counted once, and per-die copies of one
#'   counter (multi-die AMD parts) are counted once. Most distributions let
#'   only the administrator read `energy_uj`; the row then says "not
#'   readable: needs read access to the powercap files", and the estimate,
#'   when one is made, uses the package's rated power from its
#'   `constraint_0_max_power_uw`.
#' * Linux, GPUs: AMD (`amdgpu`) from the hwmon `power1_average` or
#'   `power1_input` (microwatts) and `gpu_busy_percent`; Intel (`i915`,
#'   `xe`) from the hwmon `energy1_input` (microjoules) where the kernel
#'   has it. Intel's integrated graphics are part of the CPU package
#'   figure and are not added again.
#' * NVIDIA GPUs (Linux, Windows): the driver's management library
#'   (`libnvidia-ml.so.1`, `nvml.dll`), loaded at run time: the total energy
#'   counter (millijoules), else the power (milliwatts) integrated; the GPU's
#'   utilisation and the process's share of its compute samples. Without the
#'   library, on Linux, `nvidia-smi` (power, utilisation, name).
#' * macOS on Apple silicon: the energy counters of the IOReport "Energy
#'   Model" group (CPU, GPU, neural engine and DRAM), read without
#'   administrator rights. IOReport is not a documented public interface of
#'   Apple's, so it may change with any release; `skip = "ioreport"` turns
#'   it off. A counter that stops changing before the end of the interval,
#'   or whose first step carries energy from before it (a fault reported on
#'   newer macOS), is not reported; a pause caught up within the interval
#'   is spread over it. The GPU's utilisation comes from `ioreg`.
#' * Windows, CPU and memory: the Energy Meter Interface, the channels
#'   named PKG (CPU packages) and DRAM (memory); other channels are listed
#'   and not added. Whether a machine has such a meter, and whether an
#'   ordinary user may open it, varies; the row says which.
#' * The whole machine: the battery's discharge (Linux
#'   `/sys/class/power_supply`, macOS `ioreg`, Windows the battery's
#'   discharge rate), only while it runs on battery. It stands beside the
#'   devices and is not added to them.
#' * The process's share of a CPU package: the CPU time of `pid` and its
#'   descendants (Linux; on macOS and Windows the calling process with its
#'   finished children) over the busy CPU time of the machine, times the
#'   package's energy.
#'
#' Not measurable anywhere here: CPU energy on Linux without read access to
#' the powercap files; most virtual machines and cloud instances (no
#' counters are passed through); memory on CPUs without a DRAM domain;
#' integrated GPUs apart from the CPU package; Intel GPUs on Windows and
#' macOS, and GPUs of other vendors on Windows; CPUs of Intel Macs.
#'
#' Two readings: `a <- hwpwr(0)` before a computation and `hwpwr(since = a)`
#' after it give the energy of the computation. Counters are exact over any
#' interval (but for a counter that could wrap more than once between the
#' readings, which the row notes); a device that gives only power is
#' integrated from the two readings alone. Over `dt` the devices are read
#' every `step` seconds. Nothing is written, nothing runs in the background
#' and no network is used.
#'
#' @param dt Seconds to measure over: one number from 0. With 0, one reading
#'   (power where a device gives it, for a later `since`). Not used with
#'   `since`.
#' @param since `NULL`, or an earlier result of `hwpwr()` in this session:
#'   the interval runs from that reading to now.
#' @param step Seconds between readings within `dt`: one number above 0.
#' @param rated `NULL`, or rated powers in watts by kind, `c(cpu = , gpu =
#'   )`: each CPU package's and each GPU's, used for an estimate where a
#'   device gives a utilisation but no power.
#' @param pid The process whose share is reported: a process id.
#' @param skip Sources not to use: any of `"rapl"`, `"gpu"`, `"battery"`,
#'   `"proc"`, `"nvml"`, `"smi"`, `"ioreport"`, `"ioreg"`, `"emi"`.
#' @param paths Locations to read instead of the system's, by name: `root`
#'   (a folder holding a copy of `sys` and `proc`), `nvml` (the NVIDIA
#'   library), `smi` (the `nvidia-smi` program), `ioreg` (the `ioreg`
#'   program); `""` turns a program off.
#' @return A data frame: `kind` (`cpu`, `memory`, `gpu`, `accelerator`,
#'   `meter`, `machine` or `total`), `vendor`, `name`, `scope`, `watts`
#'   (the mean over the interval), `joules`, `util` (busy share), `how`,
#'   `flo`, `fhi` and `why`. Attributes: `dt`, the interval in seconds, and
#'   `snap`, the last reading, for `since`.
#' @references
#'   Linux kernel documentation: Power Capping Framework; sysfs ABI of the
#'   power supply class; AMDgpu thermal and power interfaces; sysfs ABI of
#'   the i915 and xe hwmon interfaces.
#'
#'   NVIDIA. NVML API Reference Guide.
#'
#'   Microsoft. Energy Meter Interface; IOCTL_BATTERY_QUERY_STATUS.
#'
#'   CodeCarbon documentation: Accuracy.
#' @examples
#' r <- hwpwr(0.5)
#' r[, c("kind", "name", "scope", "watts", "joules", "how")]
#' a <- hwpwr(0)
#' x <- sum(sqrt(seq_len(1e6)))
#' b <- hwpwr(since = a)
#' b[b$kind == "total", c("joules", "why")]
#' @export
hwpwr <- function(
  dt = 1,
  since = NULL,
  step = 0.25,
  rated = NULL,
  pid = Sys.getpid(),
  skip = character(),
  paths = character()
) {
  bad <- function(...) stop(paste0("hwpwr: ", ...), call. = FALSE)
  one <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!one(dt) || dt < 0) {
    bad("`dt` must be one number from 0.")
  }
  if (!one(step) || step <= 0) {
    bad("`step` must be one number above 0.")
  }
  if (!one(pid) || pid != round(pid) || pid < 0) {
    bad("`pid` must be a process id.")
  }
  src <- c(
    "rapl",
    "gpu",
    "battery",
    "proc",
    "nvml",
    "smi",
    "ioreport",
    "ioreg",
    "emi"
  )
  if (!is.character(skip) || anyNA(skip) || !all(skip %in% src)) {
    bad("`skip` must name sources among ", paste(src, collapse = ", "), ".")
  }
  lnx <- identical(Sys.info()[["sysname"]], "Linux")
  pth <- c(
    root = "",
    nvml = "",
    smi = if (lnx) "nvidia-smi" else "",
    ioreg = if (file.exists("/usr/sbin/ioreg")) "/usr/sbin/ioreg" else ""
  )
  if (length(paths)) {
    nm <- names(paths)
    if (
      !is.character(paths) ||
        anyNA(paths) ||
        is.null(nm) ||
        !all(nm %in% names(pth)) ||
        anyDuplicated(nm)
    ) {
      bad("`paths` must be named strings among root, nvml, smi and ioreg.")
    }
    pth[nm] <- paths
  }
  rt <- c(cpu = NA_real_, gpu = NA_real_)
  if (!is.null(rated)) {
    nm <- names(rated)
    if (
      !is.numeric(rated) ||
        is.null(nm) ||
        !all(nm %in% names(rt)) ||
        anyDuplicated(nm) ||
        any(!is.finite(rated)) ||
        any(rated < 0)
    ) {
      bad("`rated` must be watts from 0 named cpu or gpu.")
    }
    rt[nm] <- rated
  }
  if (!is.null(since)) {
    sn <- attr(since, "snap")
    if (
      !is.data.frame(since) ||
        !is.character(sn) ||
        length(sn) != 1L ||
        is.na(sn)
    ) {
      bad("`since` must be an earlier result of hwpwr().")
    }
  }
  rd <- function() {
    .Call(
      C_clab_r_pwrs,
      pth[["root"]],
      pth[["nvml"]],
      pth[["smi"]],
      pth[["ioreg"]],
      skip,
      as.integer(pid)
    )
  }
  if (is.null(since)) {
    s <- rd()
    if (dt > 0) {
      end <- proc.time()[[3L]] + dt
      repeat {
        left <- end - proc.time()[[3L]]
        if (left <= 0) {
          break
        }
        Sys.sleep(min(step, left))
        s <- c(s, rd())
      }
    }
  } else {
    s <- c(sn, rd())
  }
  r <- .Call(C_clab_r_pwrd, s, unname(rt))
  attr(r, "snap") <- s[length(s)]
  r
}
