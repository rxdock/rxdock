Target platforms
================

Primary and secondary target platforms
--------------------------------------

Primary target platforms are:

* Linux on AMD64 (x86-64) with the GCC compiler
* FreeBSD on AMD64 with the Clang compiler

These platforms are used for |Dock| development on an ongoing basis.
Specifically, the target operating system and compiler versions are the latest
long-term supported releases of the main Linux distributions (e.g. Red Hat
Enterprise Linux and CentOS, SUSE Linux Enterprise and openSUSE Leap, Ubuntu)
and the latest point release of FreeBSD-STABLE.

Secondary target platforms are:

* macOS on x86-64 with Xcode (Clang compiler)
* Windows on x86-64 with MinGW (GCC compiler)

The feature set on the secondary target platforms will match the primary
platforms and bugs that affect these platforms will be treated as release
blockers.

Tertiary and quaternary target platforms
----------------------------------------

Tertiary target platforms are:

* Windows on x86-64 with MSVC compiler
* Linux on POWER8 and POWER9 (both big- and little-endian: power64, ppc64,
  power64le, ppc64le) with GCC and Clang compilers
* other modern Unix-like operating systems (e.g. Illumos) with the GCC and
  Clang compilers
* other compilers on Linux, FreeBSD, macOS, and Windows (e.g. PGI)

Bugs that affect the tertiary target platforms will be fixed if feasible. The
care will be taken to make sure that the code compiles and that the basic
functionality works. Realistically, most of the features available on the
primary and the secondary target platforms will be provided on these platforms
as well, but no promises will be made in advance.

Quaternary target platforms are:

* Linux, FreeBSD, macOS, and Windows on 32-bit x86
* Linux and FreeBSD on ARMv7 (armhf, armhfp) and ARMv8 (arm64, aarch64)
* Linux and FreeBSD on POWER7 and older processors

Bugs that affect the quaternary target platforms will be considered for fixing
if an interested user contributes a patch.
