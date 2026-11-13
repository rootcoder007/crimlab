# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#
# Installs the crimlab command from a GitHub release of
# rootcoder007/crimlab on Windows, into $HOME\.local (no administrator).
#
#   irm https://github.com/rootcoder007/crimlab/releases/latest/download/install.ps1 | iex
#   powershell -ExecutionPolicy Bypass -File install.ps1
#
# In order, stopping at the first failure: downloads the release's
# SHA256SUMS (with CRIMLAB_SIGNER set, checks its SSH signature
# SHA256SUMS.sig), downloads crimlab-<version>-windows-x86_64.zip and
# checks its SHA-256 against SHA256SUMS, checks every path in it stays
# inside its top folder, and copies bin\ into the prefix. Windows on ARM
# runs the x86-64 build under emulation.
#
# Environment: CRIMLAB_PREFIX (default $HOME\.local), CRIMLAB_VERSION
# (default: the latest release), CRIMLAB_SIGNER (the release key as an
# allowed-signers line; the release signing hook, as in install.sh),
# CRIMLAB_BASE (default https://github.com/rootcoder007/crimlab/releases;
# https, or http to localhost for tests).
#
# Errors are thrown, never `exit`, so `irm ... | iex` does not close the
# window it runs in.

function Install-Crimlab {
  $ErrorActionPreference = 'Stop'
  Set-StrictMode -Version 3
  $repo = 'rootcoder007/crimlab'
  $defsigner = ''

  function Say([string]$m) { [Console]::Error.WriteLine("crimlab install: $m") }
  function Fail([string]$m) { throw "crimlab install: error: $m; nothing was installed" }
  function Get-File([string]$url, [string]$dest) {
    $old = $ProgressPreference
    $ProgressPreference = 'SilentlyContinue'
    try { Invoke-WebRequest -UseBasicParsing -Uri $url -OutFile $dest } catch { return $false }
    finally { $ProgressPreference = $old }
    return $true
  }
  try {
    [Net.ServicePointManager]::SecurityProtocol =
      [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
  } catch { }

  $prefix = $env:CRIMLAB_PREFIX
  if (-not $prefix) { $prefix = Join-Path $HOME '.local' }
  $ver = $env:CRIMLAB_VERSION
  if ($ver) { $ver = $ver.TrimStart('v') }
  if ($ver -and $ver -notmatch '^[0-9]+\.[0-9]+\.[0-9]+$') { Fail "not a release version: $ver" }
  $signer = $env:CRIMLAB_SIGNER
  if (-not $signer) { $signer = $defsigner }
  $base = $env:CRIMLAB_BASE
  if (-not $base) { $base = "https://github.com/$repo/releases" }
  $base = $base.TrimEnd('/')
  if ($base -notmatch '^https://' -and
      $base -notmatch '^http://(127\.0\.0\.1|localhost)(:[0-9]{1,5})?(/[^@]*)?$') {
    Fail "CRIMLAB_BASE must be https (or http to localhost): $base"
  }
  $target = 'windows-x86_64'
  if ($ver) { $url = "$base/download/v$ver" } else { $url = "$base/latest/download" }

  $tmp = Join-Path ([IO.Path]::GetTempPath()) ('crimlab-' + [guid]::NewGuid().ToString('N'))
  New-Item -ItemType Directory -Path $tmp | Out-Null
  try {
    $sums = Join-Path $tmp 'SHA256SUMS'
    Say "fetching $url/SHA256SUMS"
    if (-not (Get-File "$url/SHA256SUMS" $sums)) { Fail 'cannot download SHA256SUMS' }
    $table = @{}
    foreach ($line in Get-Content -LiteralPath $sums) {
      if ($line -match '^([0-9A-Fa-f]{64}) [ *]?(\S+)$') {
        if ($table.ContainsKey($Matches[2])) { Fail "SHA256SUMS lists $($Matches[2]) twice" }
        $table[$Matches[2]] = $Matches[1].ToLowerInvariant()
      }
    }
    if ($ver) {
      $name = "crimlab-$ver-$target.zip"
      if (-not $table.ContainsKey($name)) { Fail "SHA256SUMS of v$ver lists no $name" }
    } else {
      $names = @($table.Keys | Where-Object { $_ -match ('^crimlab-[0-9]+\.[0-9]+\.[0-9]+-' + $target + '\.zip$') })
      if ($names.Count -ne 1) { Fail "the latest release lists $($names.Count) archives for $target, not one" }
      $name = $names[0]
      $ver = $name.Substring(8, $name.Length - 8 - $target.Length - 5)
    }
    $want = $table[$name]
    if ($signer) {
      $kg = Get-Command ssh-keygen -ErrorAction SilentlyContinue
      if (-not $kg) { Fail 'CRIMLAB_SIGNER is set but ssh-keygen (OpenSSH) is not installed' }
      $sig = Join-Path $tmp 'SHA256SUMS.sig'
      if (-not (Get-File "$url/SHA256SUMS.sig" $sig)) { Fail 'cannot download SHA256SUMS.sig' }
      $allowed = Join-Path $tmp 'signers'
      [IO.File]::WriteAllText($allowed, "$signer`n")
      $pr = Start-Process -FilePath $kg.Source -NoNewWindow -Wait -PassThru `
        -ArgumentList @('-Y', 'verify', '-f', "`"$allowed`"", '-I', 'crimlab', '-n', 'file', '-s', "`"$sig`"") `
        -RedirectStandardInput $sums -RedirectStandardOutput (Join-Path $tmp 'ssh.out') `
        -RedirectStandardError (Join-Path $tmp 'ssh.err')
      if ($pr.ExitCode -ne 0) { Fail 'the SSH signature of SHA256SUMS does not verify' }
      Say 'SSH signature of SHA256SUMS: good'
    } else {
      Say 'no release key set (CRIMLAB_SIGNER): who made SHA256SUMS is not checked'
    }

    $zip = Join-Path $tmp $name
    Say "fetching $url/$name"
    if (-not (Get-File "$url/$name" $zip)) { Fail "cannot download $name" }
    $got = (Get-FileHash -Algorithm SHA256 -LiteralPath $zip).Hash.ToLowerInvariant()
    if ($got -ne $want) { Fail "SHA-256 mismatch for ${name}: SHA256SUMS says $want, the download is $got" }
    Say "SHA-256 of $name matches SHA256SUMS"

    $top = "crimlab-$ver-$target"
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $z = [IO.Compression.ZipFile]::OpenRead($zip)
    try {
      foreach ($e in $z.Entries) {
        $p = $e.FullName.Replace('\', '/')
        if (-not ($p -eq "$top/" -or $p.StartsWith("$top/"))) { Fail "$name has entries outside $top/" }
        if ($p -match '(^|/)\.\.(/|$)' -or $p -match '^[A-Za-z]:') { Fail "$name has an unsafe path: $p" }
      }
    } finally { $z.Dispose() }
    $x = Join-Path $tmp 'x'
    Expand-Archive -LiteralPath $zip -DestinationPath $x
    $exe = Join-Path (Join-Path $x $top) 'bin\crimlab.exe'
    if (-not (Test-Path -LiteralPath $exe)) { Fail "$name has no bin\crimlab.exe" }

    $bin = Join-Path $prefix 'bin'
    try { New-Item -ItemType Directory -Force -Path $bin | Out-Null } catch { Fail "cannot write to $prefix; choose another with CRIMLAB_PREFIX" }
    Copy-Item -LiteralPath $exe -Destination (Join-Path $bin 'crimlab.exe') -Force
    Say "installed crimlab $ver ($target) into $prefix"
    & (Join-Path $bin 'crimlab.exe') -V
    if ($LASTEXITCODE -ne 0) { Fail 'the installed crimlab does not run' }
    $path = @($env:PATH -split ';' | ForEach-Object { $_.TrimEnd('\') })
    if ($path -notcontains $bin.TrimEnd('\')) {
      Say "add $bin to PATH (for your user, new windows):"
      Say "  [Environment]::SetEnvironmentVariable('PATH', `"$bin;`" + [Environment]::GetEnvironmentVariable('PATH', 'User'), 'User')"
    }
  } finally {
    Remove-Item -LiteralPath $tmp -Recurse -Force -ErrorAction SilentlyContinue
  }
}

Install-Crimlab
