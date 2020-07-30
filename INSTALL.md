# Installing RxDock

RxDock uses [Meson](https://mesonbuild.com/) for building.

Generic instructions: you'll need Meson and Ninja. Other dependencies
([Eigen](http://eigen.tuxfamily.org/), [PCG](https://www.pcg-random.org/),
[fmt](https://fmt.dev/), [Loguru](https://emilk.github.io/loguru/),
[nlohmann_json](https://nlohmann.github.io/json/),
[cxxopts](https://github.com/jarro2783/cxxopts), and optional
[Google Test](https://github.com/google/googletest)) will be downloaded from
[Wrap DB](https://wrapdb.mesonbuild.com/) if not available on the system.

If your platform and compiler are listed mentioned below as supported, but
you're experiencing issues,
[please let us know](https://gitlab.com/rxdock/rxdock/issues).

## Linux, FreeBSD, macOS, illumos, and Windows with MinGW/Cygwin

C++11 support in the compiler is mandatory, so
[GCC 4.8+](https://www.gnu.org/software/gcc/projects/cxx-status.html) and
[Clang 3.3+](https://clang.llvm.org/cxx_status.html)
([Xcode 5+](https://developer.apple.com/library/archive/releasenotes/DeveloperTools/RN-Xcode/Chapters/Introduction.html#//apple_ref/doc/uid/TP40001051-CH1-SW1)
on macOS) are supported.

Download a release, a snapshot, or clone the repository and `cd` into the
directory. Then

```
$ cd /wherever/your/rxdock/directory/is
$ meson -Dbuildtype=release -Dtests=true builddir
$ ninja -C builddir
$ ninja -C builddir test
```

At the moment installation is only supported on Unix-like operating systems.
Either use `sudo` to do `ninja install` such as

```
$ sudo ninja -C builddir install
```

or `su` to root and then run `ninja install` like

```
# ninja -C builddir install
```

## Windows with Visual Studio

C++11 support in the compiler is mandatory, so
[MSVC 2015+](https://docs.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance?view=vs-2015)
is supported.

Meson has a Visual Studio backend for generating Visual Studio projects that can
be built with `msbuild`. The commands are

```
> meson builddir --backend=vs
> cd builddir
> msbuild RxDock.sln /m
> msbuild RUN_TESTS.vcxproj
```
