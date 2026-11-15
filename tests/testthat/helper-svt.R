# Helpers of the scaling tests (test-svt.R, test-svb.R).

# The columns of x with the signs of the columns of ref: the sign of a
# scaling dimension is arbitrary.
colsgn <- function(x, ref) {
  x <- as.matrix(x)
  ref <- matrix(ref, nrow(x))
  for (k in seq_len(ncol(x))) {
    i <- which.max(abs(ref[, k]))
    if (sign(x[i, k]) != sign(ref[i, k])) x[, k] <- -x[, k]
  }
  x
}

# A simulated roll call of 30 legislators and 20 votes (1 yea, 0 nay, n
# absent), from a two-dimensional probit model with 5% absences; the
# reference values of the EM and Clinton-Jackman-Rivers tests were
# computed on it.
em_votes <- function() {
  rows <- c(
    "101110111101n1011110",
    "11111001110111111010",
    "1010n00n11111101n001",
    "00001100011n00100110",
    "n0100000000101011010",
    "11000010011100110011",
    "11101101101110010011",
    "111100111111110n1100",
    "00001001n00000010000",
    "111110011111111111n0",
    "000010n01000001n0000",
    "001000n0011001001010",
    "n0000100010010100n01",
    "00001000011010nn0001",
    "0000000n000010000001",
    "10010000101101n10011",
    "00000100100011100100",
    "11111111100111111110",
    "11110101111011111110",
    "1001101110011000n010",
    "00001011000101000110",
    "00001100000110100000",
    "00000100001000000010",
    "11110101100111111111",
    "11111100111111011001",
    "11n10011111111011111",
    "01001000001100100000",
    "0n100101n11111110010",
    "1010011n111010001110",
    "n0010100001n01000100"
  )
  v <- do.call(rbind, strsplit(rows, ""))
  matrix(suppressWarnings(as.numeric(v)), 30, 20)
}
