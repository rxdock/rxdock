.. _build-system:

Build system
============

|Dock| uses `the Meson build system <https://mesonbuild.com/>`__.

Rationale
---------

`Meson <https://mesonbuild.com/>`__ is a fast and user-friendly build system,
and therefore building and testing of |Dock| is now straightforward. But it was
not always so simple.

rDock built using tmake
^^^^^^^^^^^^^^^^^^^^^^^

When we stumbled upon rDock back in April 2019, the first thing we noticed was
its build process that required ``tmake`` command. As you might know,
`tmake <http://tmake.sourceforge.net/>`__ was developed by Trolltech (now
`Qt Company <https://www.qt.io/company>`__) during the 1990s as a build system
for cross-platform applications, specifically
`Qt Toolkit <https://www.qt.io/product>`__ and applications using it. It was
still somewhat popular in the early 2000s, but afterward, it got superseded by
`qmake <https://en.wikipedia.org/wiki/Qmake>`__ and later
`Qbs <https://en.wikipedia.org/wiki/Qbs_(build_tool)>`__. Well, as you might
observe, rDock never migrated away from tmake. Luckily, tmake could still run
on modern Linux distributions so building rDock was possible.

Once we got rDock to build, albeit, with thousands of lines of warnings (more
on that in the future posts), we decided that the first step is replacing tmake.
The replacement should be some actively maintained build system, but we also
wanted it to be cross-platform just like tmake as we intended to expand the
number of :ref:`supported platforms <target-platforms>`. rDock as we found it
ran only on Linux, but after forking as RxDock we ported it to FreeBSD, Solaris,
and Windows over a few months.

Aside from portability, we also wanted our build system to be easy to use; one
wants to spend time developing software, not writing build scripts. In the
process of choosing the replacement for tmake, we evaluated
`CMake <https://cmake.org/>`__, `Waf <https://waf.io/>`__, and
`Meson <https://mesonbuild.com/>`__.

Choosing the build system
^^^^^^^^^^^^^^^^^^^^^^^^^

There is the old saying that `Linux is user friendly, but picky who its friends
are <https://devrant.com/rants/818086/linux-is-user-friendly-its-just-very-picky-about-who-its-friends-are>`__.
CMake is similar, and unless you can sing precisely quoted passages from
`CGold <https://cgold.readthedocs.io/>`__ when woken up in the middle of the
night, you are not in that circle of friends (sorry to disappoint). To be fair,
the syntax of the `CMakeLists.txt` build configuration files has gotten simpler
over time, but CMake is still a lot more painful to use compared to Waf and
Meson. Furthermore, guides like
`An Introduction to Modern CMake <https://cliutils.gitlab.io/modern-cmake/>`__
are a great resource if, for whatever reason, CMake is the way you want to go.

`Waf <https://waf.io/>`__ is a mature build system, but it's nowhere near as big
and popular as CMake. Vedran Miletić worked with it quite a bit while developing
`ns-3 <https://www.nsnam.org/>`__ during his Ph.D. days so it was an attractive
option for |Dock| as well. Waf's build scripts are Python scripts so it's very
flexible in terms of what it can do while being easy to use. In summary, we
would likely have picked Waf over CMake if Meson didn't exist. Despite being
somewhat niche, Waf is still used by important projects such as
`mpv <https://mpv.io/>`__ and `Samba <https://www.samba.org/>`__.

Finally, `Meson <https://mesonbuild.com/>`__, a build system originally
developed by `Jussi Pakkanen <https://twitter.com/jpakkane>`__
`as a hobby since 2012 <https://youtu.be/gHdTzdXkhRY>`__. Jussi carefully
evaluated the existing solutions and concluded that building software can be
done better (`xkcd about standards <https://xkcd.com/927/>`__ comes to our mind
and `critics of Meson agree in that regard <https://www.rojtberg.net/1481/do-not-use-meson/>`__).
A good overview of Meson is
`the presentation at Linux.conf.au in 2015 <https://youtu.be/KPi0AuVpxLI>`__,
where the author credits `SCons <https://scons.org/>`__,
`GYP <https://gyp.gsrc.io/>`__, and qmake/Qbs for inspiration. The usage of Meson
started to grow steadily after that presentation:
`it started building GStreamer in 2016 <https://gstconf.ubicast.tv/videos/done-in-60-seconds-a-new-build-system-for-gstreamer/>`__
and `it was adopted by GNOME in 2017 <https://blogs.gnome.org/mclasen/2017/04/20/meson-considerations/>`__.
It's no longer a niche tool as
`it's used by 100+ other projects today <https://mesonbuild.com/Users.html>`__,
including some big names such as `Mesa <https://www.mesa3d.org/>`__,
`Linux Vendor Firmware Service <https://fwupd.org/>`__, `systemd <https://systemd.io/>`__,
and `DXVK <https://github.com/doitsujin/dxvk>`__. DOSBox is
`switching to Meson at the time of writing <https://github.com/dosbox-staging/dosbox-staging/issues/854>`__
in the staging branch.

The joy of working with Meson
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reason why we picked Meson among the three build systems was the ease of
use. It took just a few hours and 120 lines of the ``meson.build`` build
configuration file (more than a third of that being the names of the source
files) to get |Dock| building and running, and tests came soon afterward.
`This tweet from Vedran Miletić <https://twitter.com/VedranMiletic/status/1121310731184611328>`__
captures some of the initial excitement coming from replacing tmake and custom
``Makefile``\ s with Meson.

Unlike Waf, the `syntax isn't Python <https://mesonbuild.com/Syntax.html>`__ (there's
`a good reason for that <https://mesonbuild.com/FAQ.html#why-is-meson-not-just-a-python-module-so-i-could-code-my-build-setup-in-python>`__),
but it's very much Python-like and therefore easy to get started. At
`CppCon 2018 <https://cppcon.org/cppcon-2018-program/>`__, Jussi told the
audience to just
`look at the build configuration file for 10 seconds <https://youtu.be/SCZLnopmYBM?t=300>`__
without any prior knowledge of the Meson syntax; it's just so obvious what it
does. The documentation is very good, especially its coverage of
`common use cases <https://mesonbuild.com/howtox.html>`__. You get
`reproducible builds for free <https://mesonbuild.com/Reproducible-builds.html>`__.
IDE integration is getting there, e.g.
`Qt Creator <https://doc.qt.io/qtcreator/creator-project-meson.html>`__
`added Meson support in summer 2020 <https://www.qt.io/blog/qt-creator-4.13-released>`__ and
`Visual Studio Code had it for a while <https://marketplace.visualstudio.com/items?itemName=asabil.meson>`__.

It has been almost two years at the time of writing since we switched to Meson
and we are very happy with how it builds |Dock|.
`The Meson community is very active <https://github.com/mesonbuild/meson/graphs/contributors>`__
and therefore we expect to continue using Meson to build |Dock| for the foreseeable future.
