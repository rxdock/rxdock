.. _quick-and-dirty-installation:

Quick and dirty installation
============================

In this section you will have short instructions to make a typical installation
of |Dock|.

To get the full documentation of all |Dock| software package and methods, please
go to the :ref:`Reference guide <reference-guide>`.

Moreover, you can also check the following information:

* :ref:`Getting started <getting-started-guide>`: installation and validation
  instructions for first-time users.
* :ref:`Validation experiments <validation-experiments>`: instructions and
  examples for re-running the validation sets we have carried out.
* :ref:`Calculating ROC curves <calculating-roc-curves>`: tutorial for
  generating ROC curves and other statistics after running |Dock| docking jobs.

Installation in 3 steps
-----------------------

We have been able to compile |Dock| in the following Linux systems:

* CentOS 5.5 64 bits
* openSUSE 11.3 32 and 64 bits
* openSUSE 12.3 32 and 64 bits
* openSUSE 13.1 32 and 64 bits
* Ubuntu 12.04 32 and 64 bits

Step 1
^^^^^^

First of all, you will need to install several packages before compiling and
running |Dock|:

* gcc and g++ compilers version > 3.3
* make
* cppunit and cppunit-devel

.. note::

   **For Ubuntu users:**

   If you are trying to use |Dock| in Ubuntu, please note that csh shell is not
   included in a default installation. We recommend to install csh in case some
   error arises (``sudo apt-get install csh``), even with all the above-stated
   dependencies installed.

Afterwards, download the source code compressed file or get it by SVN in
:ref:`Downloads <download>` section.

Step 2
^^^^^^

Then, run the following commands:

.. code-block:: bash

   $ tar -xvzf rxdock-0.1.0.tar.gz
   $ cd rxdock-0.1.0/build/

and, for 32 bits computers:

.. code-block:: bash

   $ make linux-g++

for 64 bits computers:

.. code-block:: bash

   $ make linux-g++-64

Step 3
^^^^^^

After compiling successfully, type the following command to make a test and
check that your compiled version works good and the results are correct.

.. code-block:: bash

   $ make test

If the test has succeed, you are done, enjoy using |Dock|!

Otherwise, please check your dependencies and all the previous commands or go to
:ref:`Support Section <support>` to ask for help.

Just as a concluding remark, don't forget to set the necessary environmental
variables for running |Dock| in the command line (for example, in Bash shell):

.. code-block:: bash

   $ export RBT_ROOT=/path/to/rxdock/installation/
   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RBT_ROOT/lib
   $ export PATH=$PATH:$RBT_ROOT/bin

Installation tutorial
---------------------

We can say that |Dock| has a great set of functionalities and that it is one of
the best open-source high-throughput virtual screening (HTVS) suites packed
with all the functionalities aspiring computational chemists need, but the
major hurdle in getting new users onboard and proving it to them is the
seemingly archaic way of installing and using the software. The graphical user
interface is non-existent, the software itself cannot be found maintained in
the relevant Linux (or `macOS Homebrew <https://brew.sh/>`__, or
`FreeBSD Ports <https://www.freebsd.org/ports/>`__, etc.) repositories, you have
to set up your working directory every time you start using |Dock|, you have to
build the software from its sources, and all of that sounds very taxing and
complicated to the average user, aspiring graduate students and everyone who
are not tech-savvy.

Getting started to use scientific software can be daunting for anyone,
especially beginners or newcomers to the field. Luckily, building and setting
up |Dock| is much easier done than said, and we will walk you over through the
entire process in this blog post. This blog post is intended for absolute
beginners in Linux, scientific software in general and computational chemistry.

In this tutorial, we will be using
`Ubuntu 20.04 LTS <https://ubuntu.com/blog/ubuntu-20-04-lts-arrives>`__ as an
example, but newer Ubuntu versions should work just as well.

Installing and setting up |Dock|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The good thing about command-line interface (CLI) software is that everything
can be easily replicated and everyone has the same user experience. Here, you
will be able to follow the guide simply by copy-pasting the commands in your
terminal that comes with every Ubuntu desktop installation. Just open it up and
follow the instruction further in the blog post.

Before we start, let's convince ourselves that we're running Ubuntu 20.04:

.. code-block:: bash

   $ cat /etc/os-release
   NAME="Ubuntu"
   VERSION="20.04.3 LTS (Focal Fossa)"
   ID=ubuntu
   ID_LIKE=debian
   PRETTY_NAME="Ubuntu 20.04.3 LTS"
   VERSION_ID="20.04"
   HOME_URL="https://www.ubuntu.com/"
   SUPPORT_URL="https://help.ubuntu.com/"
   BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
   PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
   VERSION_CODENAME=focal
   UBUNTU_CODENAME=focal

Alternatively, if you have `neofetch <https://github.com/dylanaraps/neofetch>`__
installed, you can use it instead:

.. code-block:: bash

   $ neofetch
               .-/+oossssoo+/-.               pnikolic@rxtxboss
           `:+ssssssssssssssssss+:`           -----------------
         -+ssssssssssssssssssyyssss+-         OS: Ubuntu 20.04.3 LTS x86_64
       .ossssssssssssssssssdMMMNysssso.       Host: KVM/QEMU (Standard PC (Q35 + ICH9, 2009) pc-q35-6.
      /ssssssssssshdmmNNmmyNMMMMhssssss/      Kernel: 5.4.0-96-generic
     +ssssssssshmydMMMMMMMNddddyssssssss+     Uptime: 55 mins
    /sssssssshNMMMyhhyyyyhmNMMMNhssssssss/    Packages: 681 (dpkg), 5 (snap)
   .ssssssssdMMMNhsssssssssshNMMMdssssssss.   Shell: bash 5.0.17
   +sssshhhyNMMNyssssssssssssyNMMMysssssss+   Resolution: 1024x768
   ossyNMMMNyMMhsssssssssssssshmmmhssssssso   Terminal: /dev/pts/0
   ossyNMMMNyMMhsssssssssssssshmmmhssssssso   CPU: AMD EPYC (with IBPB) (2) @ 3.393GHz
   +sssshhhyNMMNyssssssssssssyNMMMysssssss+   GPU: 00:01.0 Red Hat, Inc. QXL paravirtual graphic card
   .ssssssssdMMMNhsssssssssshNMMMdssssssss.   Memory: 161MiB / 3931MiB
    /sssssssshNMMMyhhyyyyhdNMMMNhssssssss/
     +sssssssssdmydMMMMMMMMddddyssssssss+
      /ssssssssssshdmNNNNmyNMMMMhssssss/
       .ossssssssssssssssssdMMMNysssso.
         -+sssssssssssssssssyyyssss+-
           `:+ssssssssssssssssss+:`
               .-/+oossssoo+/-.

The first thing we have to take care of is updates (in case the system isn't
already fully up-to-date). This is done with ``apt update`` and ``apt upgrade``
commands:

.. code-block:: bash

   $ sudo apt update
   [sudo] password for pnikolic:
   Hit:1 http://hr.archive.ubuntu.com/ubuntu focal InRelease
   Hit:2 http://hr.archive.ubuntu.com/ubuntu focal-updates InRelease
   Hit:3 http://hr.archive.ubuntu.com/ubuntu focal-backports InRelease
   Hit:4 http://hr.archive.ubuntu.com/ubuntu focal-security InRelease
   Reading package lists... Done
   Building dependency tree
   Reading state information... Done
   37 packages can be upgraded. Run 'apt list --upgradable' to see them.

.. code-block:: bash

   $ sudo apt upgrade
   Reading package lists... Done
   Building dependency tree
   Reading state information... Done
   Calculating upgrade... Done
   The following packages will be upgraded:
   alsa-ucm-conf cloud-init cloud-initramfs-copymods cloud-initramfs-dyn-netconf command-not-found
   libasound2 libasound2-data libdrm-common libdrm2 libnetplan0 libpam-modules libpam-modules-bin
   libpam-runtime libpam0g libprocps8 libssl1.1 libudisks2-0 linux-base netplan.io open-vm-tools
   openssl overlayroot procps python-apt-common python3-apt python3-commandnotfound
   python3-software-properties python3-update-manager rsync snapd software-properties-common
   ubuntu-advantage-tools udisks2 ufw update-manager-core update-notifier-common wget
   37 upgraded, 0 newly installed, 0 to remove and 0 not upgraded.
   Need to get 37.2 MB of archives.
   After this operation, 1,161 kB of additional disk space will be used.
   Do you want to continue? [Y/n] y
   Get:1 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libpam0g amd64 1.3.1-5ubuntu4.3 [55.4 kB]
   Get:2 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libpam-modules-bin amd64 1.3.1-5ubuntu4.3 [41.2 kB]
   Get:3 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libpam-modules amd64 1.3.1-5ubuntu4.3 [260 kB]
   Get:4 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 rsync amd64 3.1.3-8ubuntu0.1 [318 kB]
   Get:5 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 python-apt-common all 2.0.0ubuntu0.20.04.6 [17.1 kB]
   Get:6 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 python3-apt amd64 2.0.0ubuntu0.20.04.6 [154 kB]
   Get:7 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 update-manager-core all 1:20.04.10.9 [11.5 kB]
   Get:8 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 python3-update-manager all 1:20.04.10.9 [38.1 kB]
   Get:9 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 ubuntu-advantage-tools amd64 27.5~20.04.1 [861 kB]
   Get:10 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 update-notifier-common all 3.192.30.10 [133 kB]
   Get:11 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libprocps8 amd64 2:3.3.16-1ubuntu2.3 [33.0 kB]
   Get:12 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 procps amd64 2:3.3.16-1ubuntu2.3 [233 kB]
   Get:13 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libdrm-common all 2.4.105-3~20.04.2 [5,552 B]
   Get:14 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libdrm2 amd64 2.4.105-3~20.04.2 [32.3 kB]
   Get:15 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libssl1.1 amd64 1.1.1f-1ubuntu2.10 [1,322 kB]
   Get:16 http://hr.archive.ubuntu.com/ubuntu focal-updates/universe amd64 open-vm-tools amd64 2:11.3.0-2ubuntu0~ubuntu20.04.2 [647 kB]
   Get:17 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libpam-runtime all 1.3.1-5ubuntu4.3 [37.3 kB]
   Get:18 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libnetplan0 amd64 0.103-0ubuntu5~20.04.5 [53.4 kB]
   Get:19 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 netplan.io amd64 0.103-0ubuntu5~20.04.5 [125 kB]
   Get:20 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 openssl amd64 1.1.1f-1ubuntu2.10 [620 kB]
   Get:21 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 command-not-found all 20.04.5 [5,244 B]
   Get:22 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 python3-commandnotfound all 20.04.5 [10.2 kB]
   Get:23 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 ufw all 0.36-6ubuntu1 [147 kB]
   Get:24 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 wget amd64 1.20.3-1ubuntu2 [348 kB]
   Get:25 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libasound2 amd64 1.2.2-2.1ubuntu2.5 [335 kB]
   Get:26 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libasound2-data all 1.2.2-2.1ubuntu2.5 [20.1 kB]
   Get:27 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 alsa-ucm-conf all 1.2.2-1ubuntu0.11 [26.9 kB]
   Get:28 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libudisks2-0 amd64 2.8.4-1ubuntu2 [99.4 kB]
   Get:29 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 linux-base all 4.5ubuntu3.7 [17.6 kB]
   Get:30 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 software-properties-common all 0.99.9.8 [10.6 kB]
   Get:31 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 python3-software-properties all 0.99.9.8 [24.9 kB]
   Get:32 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 snapd amd64 2.51.1+20.04ubuntu2 [30.4 MB]
   Get:33 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 udisks2 amd64 2.8.4-1ubuntu2 [245 kB]
   Get:34 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 cloud-init all 21.4-0ubuntu1~20.04.1 [476 kB]
   Get:35 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 cloud-initramfs-copymods all 0.45ubuntu2 [4,180 B]
   Get:36 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 cloud-initramfs-dyn-netconf all 0.45ubuntu2 [6,700 B]
   Get:37 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 overlayroot all 0.45ubuntu2 [15.6 kB]
   Fetched 37.2 MB in 7s (5,471 kB/s)
   Extracting templates from packages: 100%
   Preconfiguring packages ...
   (Reading database ... 71625 files and directories currently installed.)
   Preparing to unpack .../libpam0g_1.3.1-5ubuntu4.3_amd64.deb ...
   Unpacking libpam0g:amd64 (1.3.1-5ubuntu4.3) over (1.3.1-5ubuntu4.2) ...
   Setting up libpam0g:amd64 (1.3.1-5ubuntu4.3) ...
   (Reading database ... 71625 files and directories currently installed.)
   Preparing to unpack .../libpam-modules-bin_1.3.1-5ubuntu4.3_amd64.deb ...
   Unpacking libpam-modules-bin (1.3.1-5ubuntu4.3) over (1.3.1-5ubuntu4.2) ...
   Setting up libpam-modules-bin (1.3.1-5ubuntu4.3) ...
   (Reading database ... 71625 files and directories currently installed.)
   Preparing to unpack .../libpam-modules_1.3.1-5ubuntu4.3_amd64.deb ...
   Unpacking libpam-modules:amd64 (1.3.1-5ubuntu4.3) over (1.3.1-5ubuntu4.2) ...
   Setting up libpam-modules:amd64 (1.3.1-5ubuntu4.3) ...
   (Reading database ... 71625 files and directories currently installed.)
   Preparing to unpack .../00-rsync_3.1.3-8ubuntu0.1_amd64.deb ...
   Unpacking rsync (3.1.3-8ubuntu0.1) over (3.1.3-8) ...
   Preparing to unpack .../01-python-apt-common_2.0.0ubuntu0.20.04.6_all.deb ...
   Unpacking python-apt-common (2.0.0ubuntu0.20.04.6) over (2.0.0ubuntu0.20.04.5) ...
   Preparing to unpack .../02-python3-apt_2.0.0ubuntu0.20.04.6_amd64.deb ...
   Unpacking python3-apt (2.0.0ubuntu0.20.04.6) over (2.0.0ubuntu0.20.04.5) ...
   Preparing to unpack .../03-update-manager-core_1%3a20.04.10.9_all.deb ...
   Unpacking update-manager-core (1:20.04.10.9) over (1:20.04.10.7) ...
   Preparing to unpack .../04-python3-update-manager_1%3a20.04.10.9_all.deb ...
   Unpacking python3-update-manager (1:20.04.10.9) over (1:20.04.10.7) ...
   Preparing to unpack .../05-ubuntu-advantage-tools_27.5~20.04.1_amd64.deb ...
   Unpacking ubuntu-advantage-tools (27.5~20.04.1) over (27.2.2~20.04.1) ...
   Preparing to unpack .../06-update-notifier-common_3.192.30.10_all.deb ...
   Unpacking update-notifier-common (3.192.30.10) over (3.192.30.9) ...
   Preparing to unpack .../07-libprocps8_2%3a3.3.16-1ubuntu2.3_amd64.deb ...
   Unpacking libprocps8:amd64 (2:3.3.16-1ubuntu2.3) over (2:3.3.16-1ubuntu2.2) ...
   Preparing to unpack .../08-procps_2%3a3.3.16-1ubuntu2.3_amd64.deb ...
   Unpacking procps (2:3.3.16-1ubuntu2.3) over (2:3.3.16-1ubuntu2.2) ...
   Preparing to unpack .../09-libdrm-common_2.4.105-3~20.04.2_all.deb ...
   Unpacking libdrm-common (2.4.105-3~20.04.2) over (2.4.105-3~20.04.1) ...
   Preparing to unpack .../10-libdrm2_2.4.105-3~20.04.2_amd64.deb ...
   Unpacking libdrm2:amd64 (2.4.105-3~20.04.2) over (2.4.105-3~20.04.1) ...
   Preparing to unpack .../11-libssl1.1_1.1.1f-1ubuntu2.10_amd64.deb ...
   Unpacking libssl1.1:amd64 (1.1.1f-1ubuntu2.10) over (1.1.1f-1ubuntu2.8) ...
   Preparing to unpack .../12-open-vm-tools_2%3a11.3.0-2ubuntu0~ubuntu20.04.2_amd64.deb ...
   Unpacking open-vm-tools (2:11.3.0-2ubuntu0~ubuntu20.04.2) over (2:11.0.5-4) ...
   Preparing to unpack .../13-libpam-runtime_1.3.1-5ubuntu4.3_all.deb ...
   Unpacking libpam-runtime (1.3.1-5ubuntu4.3) over (1.3.1-5ubuntu4.2) ...
   Setting up libpam-runtime (1.3.1-5ubuntu4.3) ...
   (Reading database ... 71656 files and directories currently installed.)
   Preparing to unpack .../00-libnetplan0_0.103-0ubuntu5~20.04.5_amd64.deb ...
   Unpacking libnetplan0:amd64 (0.103-0ubuntu5~20.04.5) over (0.102-0ubuntu1~20.04.2) ...
   Preparing to unpack .../01-netplan.io_0.103-0ubuntu5~20.04.5_amd64.deb ...
   Unpacking netplan.io (0.103-0ubuntu5~20.04.5) over (0.102-0ubuntu1~20.04.2) ...
   Preparing to unpack .../02-openssl_1.1.1f-1ubuntu2.10_amd64.deb ...
   Unpacking openssl (1.1.1f-1ubuntu2.10) over (1.1.1f-1ubuntu2.8) ...
   Preparing to unpack .../03-command-not-found_20.04.5_all.deb ...
   Unpacking command-not-found (20.04.5) over (20.04.4) ...
   Preparing to unpack .../04-python3-commandnotfound_20.04.5_all.deb ...
   Unpacking python3-commandnotfound (20.04.5) over (20.04.4) ...
   Preparing to unpack .../05-ufw_0.36-6ubuntu1_all.deb ...
   Unpacking ufw (0.36-6ubuntu1) over (0.36-6) ...
   Preparing to unpack .../06-wget_1.20.3-1ubuntu2_amd64.deb ...
   Unpacking wget (1.20.3-1ubuntu2) over (1.20.3-1ubuntu1) ...
   Preparing to unpack .../07-libasound2_1.2.2-2.1ubuntu2.5_amd64.deb ...
   Unpacking libasound2:amd64 (1.2.2-2.1ubuntu2.5) over (1.2.2-2.1ubuntu2.4) ...
   Preparing to unpack .../08-libasound2-data_1.2.2-2.1ubuntu2.5_all.deb ...
   Unpacking libasound2-data (1.2.2-2.1ubuntu2.5) over (1.2.2-2.1ubuntu2.4) ...
   Preparing to unpack .../09-alsa-ucm-conf_1.2.2-1ubuntu0.11_all.deb ...
   Unpacking alsa-ucm-conf (1.2.2-1ubuntu0.11) over (1.2.2-1ubuntu0.9) ...
   Preparing to unpack .../10-libudisks2-0_2.8.4-1ubuntu2_amd64.deb ...
   Unpacking libudisks2-0:amd64 (2.8.4-1ubuntu2) over (2.8.4-1ubuntu1) ...
   Preparing to unpack .../11-linux-base_4.5ubuntu3.7_all.deb ...
   Unpacking linux-base (4.5ubuntu3.7) over (4.5ubuntu3.6) ...
   Preparing to unpack .../12-software-properties-common_0.99.9.8_all.deb ...
   Unpacking software-properties-common (0.99.9.8) over (0.98.9.5) ...
   Preparing to unpack .../13-python3-software-properties_0.99.9.8_all.deb ...
   Unpacking python3-software-properties (0.99.9.8) over (0.98.9.5) ...
   Preparing to unpack .../14-snapd_2.51.1+20.04ubuntu2_amd64.deb ...
   Unpacking snapd (2.51.1+20.04ubuntu2) over (2.49.2+20.04) ...
   Preparing to unpack .../15-udisks2_2.8.4-1ubuntu2_amd64.deb ...
   Unpacking udisks2 (2.8.4-1ubuntu2) over (2.8.4-1ubuntu1) ...
   Preparing to unpack .../16-cloud-init_21.4-0ubuntu1~20.04.1_all.deb ...
   Unpacking cloud-init (21.4-0ubuntu1~20.04.1) over (21.2-3-g899bfaa9-0ubuntu2~20.04.1) ...
   Preparing to unpack .../17-cloud-initramfs-copymods_0.45ubuntu2_all.deb ...
   Unpacking cloud-initramfs-copymods (0.45ubuntu2) over (0.45ubuntu1) ...
   Preparing to unpack .../18-cloud-initramfs-dyn-netconf_0.45ubuntu2_all.deb ...
   Unpacking cloud-initramfs-dyn-netconf (0.45ubuntu2) over (0.45ubuntu1) ...
   Preparing to unpack .../19-overlayroot_0.45ubuntu2_all.deb ...
   Unpacking overlayroot (0.45ubuntu2) over (0.45ubuntu1) ...
   Setting up snapd (2.51.1+20.04ubuntu2) ...
   Installing new version of config file /etc/profile.d/apps-bin-path.sh ...
   snapd.failure.service is a disabled or a static unit, not starting it.
   snapd.snap-repair.service is a disabled or a static unit, not starting it.
   Setting up linux-base (4.5ubuntu3.7) ...
   Setting up alsa-ucm-conf (1.2.2-1ubuntu0.11) ...
   Setting up libssl1.1:amd64 (1.1.1f-1ubuntu2.10) ...
   Setting up libnetplan0:amd64 (0.103-0ubuntu5~20.04.5) ...
   Setting up ufw (0.36-6ubuntu1) ...
   Setting up netplan.io (0.103-0ubuntu5~20.04.5) ...
   Setting up libasound2-data (1.2.2-2.1ubuntu2.5) ...
   Setting up cloud-initramfs-copymods (0.45ubuntu2) ...
   Setting up python-apt-common (2.0.0ubuntu0.20.04.6) ...
   Setting up libasound2:amd64 (1.2.2-2.1ubuntu2.5) ...
   Setting up cloud-initramfs-dyn-netconf (0.45ubuntu2) ...
   Setting up openssl (1.1.1f-1ubuntu2.10) ...
   Setting up libdrm-common (2.4.105-3~20.04.2) ...
   Setting up overlayroot (0.45ubuntu2) ...
   Setting up libprocps8:amd64 (2:3.3.16-1ubuntu2.3) ...
   Setting up rsync (3.1.3-8ubuntu0.1) ...
   Setting up libudisks2-0:amd64 (2.8.4-1ubuntu2) ...
   Setting up udisks2 (2.8.4-1ubuntu2) ...
   Setting up wget (1.20.3-1ubuntu2) ...
   Setting up python3-apt (2.0.0ubuntu0.20.04.6) ...
   Setting up python3-software-properties (0.99.9.8) ...
   Setting up python3-commandnotfound (20.04.5) ...
   Setting up libdrm2:amd64 (2.4.105-3~20.04.2) ...
   Setting up open-vm-tools (2:11.3.0-2ubuntu0~ubuntu20.04.2) ...
   Installing new version of config file /etc/vmware-tools/tools.conf.example ...
   Installing new version of config file /etc/vmware-tools/vgauth.conf ...
   Removing obsolete conffile /etc/vmware-tools/vm-support ...
   Setting up python3-update-manager (1:20.04.10.9) ...
   Setting up procps (2:3.3.16-1ubuntu2.3) ...
   Setting up ubuntu-advantage-tools (27.5~20.04.1) ...
   Installing new version of config file /etc/logrotate.d/ubuntu-advantage-tools ...
   Installing new version of config file /etc/ubuntu-advantage/help_data.yaml ...
   Installing new version of config file /etc/ubuntu-advantage/uaclient.conf ...
   Created symlink /etc/systemd/system/timers.target.wants/ua-timer.timer → /lib/systemd/system/ua-timer.timer.
   Created symlink /etc/systemd/system/multi-user.target.wants/ua-license-check.path → /lib/systemd/system/ua-license-check.path.
   Setting up cloud-init (21.4-0ubuntu1~20.04.1) ...
   Installing new version of config file /etc/cloud/cloud.cfg ...
   Installing new version of config file /etc/cloud/templates/hosts.alpine.tmpl ...
   Installing new version of config file /etc/cloud/templates/hosts.debian.tmpl ...
   Installing new version of config file /etc/cloud/templates/resolv.conf.tmpl ...
   Installing new version of config file /etc/cloud/templates/sources.list.debian.tmpl ...
   Created symlink /etc/systemd/system/cloud-init.target.wants/cloud-init-hotplugd.socket → /lib/systemd/system/cloud-init-hotplugd.socket.
   Setting up software-properties-common (0.99.9.8) ...
   Setting up command-not-found (20.04.5) ...
   Setting up update-manager-core (1:20.04.10.9) ...
   Setting up update-notifier-common (3.192.30.10) ...
   Processing triggers for rsyslog (8.2001.0-1ubuntu1.1) ...
   Processing triggers for systemd (245.4-4ubuntu3.15) ...
   Processing triggers for man-db (2.9.1-1) ...
   Processing triggers for dbus (1.12.16-2ubuntu2.1) ...
   Processing triggers for install-info (6.7.0.dfsg.2-5) ...
   Processing triggers for mime-support (3.64ubuntu1) ...
   Processing triggers for libc-bin (2.31-0ubuntu9.2) ...
   Processing triggers for initramfs-tools (0.136ubuntu6.6) ...
   update-initramfs: Generating /boot/initrd.img-5.4.0-96-generic

Let's continue with installing tools for obtaining and building |Dock|,
specifically `Git <https://git-scm.com/>`__, `Meson <https://mesonbuild.com/>`__,
and `GNU C++ compiler <https://gcc.gnu.org/>`__. Additionally, `CMake <https://cmake.org/>`__
and `pkg-config <https://www.freedesktop.org/wiki/Software/pkg-config/>`__ are
used for finding dependencies on the system. The installation is done using the
``apt install`` command:

.. code-block:: bash

   $ sudo apt install git meson build-essential cmake pkg-config
   Reading package lists... Done
   Building dependency tree
   Reading state information... Done
   git is already the newest version (1:2.25.1-1ubuntu3.2).
   git set to manually installed.
   The following additional packages will be installed:
   binutils binutils-common binutils-x86-64-linux-gnu cmake-data cpp cpp-9 dpkg-dev fakeroot g++
   g++-9 gcc gcc-9 gcc-9-base libalgorithm-diff-perl libalgorithm-diff-xs-perl
   libalgorithm-merge-perl libasan5 libatomic1 libbinutils libc-dev-bin libc6-dev libcc1-0
   libcrypt-dev libctf-nobfd0 libctf0 libdpkg-perl libfakeroot libfile-fcntllock-perl libgcc-9-dev
   libgomp1 libisl22 libitm1 libjsoncpp1 liblsan0 libmpc3 libquadmath0 librhash0 libstdc++-9-dev
   libtsan0 libubsan1 linux-libc-dev make manpages-dev ninja-build
   Suggested packages:
   binutils-doc cmake-doc cpp-doc gcc-9-locales debian-keyring g++-multilib g++-9-multilib gcc-9-doc
   gcc-multilib autoconf automake libtool flex bison gdb gcc-doc gcc-9-multilib glibc-doc bzr
   libstdc++-9-doc make-doc
   The following NEW packages will be installed:
   binutils binutils-common binutils-x86-64-linux-gnu build-essential cmake cmake-data cpp cpp-9
   dpkg-dev fakeroot g++ g++-9 gcc gcc-9 gcc-9-base libalgorithm-diff-perl libalgorithm-diff-xs-perl
   libalgorithm-merge-perl libasan5 libatomic1 libbinutils libc-dev-bin libc6-dev libcc1-0
   libcrypt-dev libctf-nobfd0 libctf0 libdpkg-perl libfakeroot libfile-fcntllock-perl libgcc-9-dev
   libgomp1 libisl22 libitm1 libjsoncpp1 liblsan0 libmpc3 libquadmath0 librhash0 libstdc++-9-dev
   libtsan0 libubsan1 linux-libc-dev make manpages-dev meson ninja-build pkg-config
   0 upgraded, 48 newly installed, 0 to remove and 0 not upgraded.
   Need to get 49.0 MB of archives.
   After this operation, 220 MB of additional disk space will be used.
   Do you want to continue? [Y/n] y
   Get:1 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 binutils-common amd64 2.34-6ubuntu1.3 [207 kB]
   Get:2 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libbinutils amd64 2.34-6ubuntu1.3 [474 kB]
   Get:3 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libctf-nobfd0 amd64 2.34-6ubuntu1.3 [47.4 kB]
   Get:4 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libctf0 amd64 2.34-6ubuntu1.3 [46.6 kB]
   Get:5 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 binutils-x86-64-linux-gnu amd64 2.34-6ubuntu1.3 [1,613 kB]
   Get:6 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 binutils amd64 2.34-6ubuntu1.3 [3,380 B]
   Get:7 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libc-dev-bin amd64 2.31-0ubuntu9.2 [71.8 kB]
   Get:8 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 linux-libc-dev amd64 5.4.0-96.109 [1,114 kB]
   Get:9 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libcrypt-dev amd64 1:4.4.10-10ubuntu4 [104 kB]
   Get:10 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libc6-dev amd64 2.31-0ubuntu9.2 [2,520 kB]
   Get:11 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 gcc-9-base amd64 9.3.0-17ubuntu1~20.04 [19.1 kB]
   Get:12 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libisl22 amd64 0.22.1-1 [592 kB]
   Get:13 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libmpc3 amd64 1.1.0-1 [40.8 kB]
   Get:14 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 cpp-9 amd64 9.3.0-17ubuntu1~20.04 [7,494 kB]
   Get:15 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 cpp amd64 4:9.3.0-1ubuntu2 [27.6 kB]
   Get:16 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libcc1-0 amd64 10.3.0-1ubuntu1~20.04 [48.8 kB]
   Get:17 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libgomp1 amd64 10.3.0-1ubuntu1~20.04 [102 kB]
   Get:18 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libitm1 amd64 10.3.0-1ubuntu1~20.04 [26.2 kB]
   Get:19 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libatomic1 amd64 10.3.0-1ubuntu1~20.04 [9,284 B]
   Get:20 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libasan5 amd64 9.3.0-17ubuntu1~20.04 [394 kB]
   Get:21 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 liblsan0 amd64 10.3.0-1ubuntu1~20.04 [835 kB]
   Get:22 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libtsan0 amd64 10.3.0-1ubuntu1~20.04 [2,009 kB]
   Get:23 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libubsan1 amd64 10.3.0-1ubuntu1~20.04 [784 kB]
   Get:24 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libquadmath0 amd64 10.3.0-1ubuntu1~20.04 [146 kB]
   Get:25 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libgcc-9-dev amd64 9.3.0-17ubuntu1~20.04 [2,360 kB]
   Get:26 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 gcc-9 amd64 9.3.0-17ubuntu1~20.04 [8,241 kB]
   Get:27 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 gcc amd64 4:9.3.0-1ubuntu2 [5,208 B]
   Get:28 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 libstdc++-9-dev amd64 9.3.0-17ubuntu1~20.04 [1,714 kB]
   Get:29 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 g++-9 amd64 9.3.0-17ubuntu1~20.04 [8,405 kB]
   Get:30 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 g++ amd64 4:9.3.0-1ubuntu2 [1,604 B]
   Get:31 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 make amd64 4.2.1-1.2 [162 kB]
   Get:32 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libdpkg-perl all 1.19.7ubuntu3 [230 kB]
   Get:33 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 dpkg-dev all 1.19.7ubuntu3 [679 kB]
   Get:34 http://hr.archive.ubuntu.com/ubuntu focal-updates/main amd64 build-essential amd64 12.8ubuntu1.1 [4,664 B]
   Get:35 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 cmake-data all 3.16.3-1ubuntu1 [1,612 kB]
   Get:36 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libjsoncpp1 amd64 1.7.4-3.1ubuntu2 [75.6 kB]
   Get:37 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 librhash0 amd64 1.3.9-1 [113 kB]
   Get:38 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 cmake amd64 3.16.3-1ubuntu1 [3,669 kB]
   Get:39 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libfakeroot amd64 1.24-1 [25.7 kB]
   Get:40 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 fakeroot amd64 1.24-1 [62.6 kB]
   Get:41 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libalgorithm-diff-perl all 1.19.03-2 [46.6 kB]
   Get:42 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libalgorithm-diff-xs-perl amd64 0.04-6 [11.3 kB]
   Get:43 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libalgorithm-merge-perl all 0.08-3 [12.0 kB]
   Get:44 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 libfile-fcntllock-perl amd64 0.22-3build4 [33.1 kB]
   Get:45 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 manpages-dev all 5.05-1 [2,266 kB]
   Get:46 http://hr.archive.ubuntu.com/ubuntu focal/universe amd64 ninja-build amd64 1.10.0-1build1 [107 kB]
   Get:47 http://hr.archive.ubuntu.com/ubuntu focal/universe amd64 meson all 0.53.2-2ubuntu2 [376 kB]
   Get:48 http://hr.archive.ubuntu.com/ubuntu focal/main amd64 pkg-config amd64 0.29.1-0ubuntu4 [45.5 kB]
   Fetched 49.0 MB in 9s (5,612 kB/s)
   Extracting templates from packages: 100%
   Selecting previously unselected package binutils-common:amd64.
   (Reading database ... 71680 files and directories currently installed.)
   Preparing to unpack .../00-binutils-common_2.34-6ubuntu1.3_amd64.deb ...
   Unpacking binutils-common:amd64 (2.34-6ubuntu1.3) ...
   Selecting previously unselected package libbinutils:amd64.
   Preparing to unpack .../01-libbinutils_2.34-6ubuntu1.3_amd64.deb ...
   Unpacking libbinutils:amd64 (2.34-6ubuntu1.3) ...
   Selecting previously unselected package libctf-nobfd0:amd64.
   Preparing to unpack .../02-libctf-nobfd0_2.34-6ubuntu1.3_amd64.deb ...
   Unpacking libctf-nobfd0:amd64 (2.34-6ubuntu1.3) ...
   Selecting previously unselected package libctf0:amd64.
   Preparing to unpack .../03-libctf0_2.34-6ubuntu1.3_amd64.deb ...
   Unpacking libctf0:amd64 (2.34-6ubuntu1.3) ...
   Selecting previously unselected package binutils-x86-64-linux-gnu.
   Preparing to unpack .../04-binutils-x86-64-linux-gnu_2.34-6ubuntu1.3_amd64.deb ...
   Unpacking binutils-x86-64-linux-gnu (2.34-6ubuntu1.3) ...
   Selecting previously unselected package binutils.
   Preparing to unpack .../05-binutils_2.34-6ubuntu1.3_amd64.deb ...
   Unpacking binutils (2.34-6ubuntu1.3) ...
   Selecting previously unselected package libc-dev-bin.
   Preparing to unpack .../06-libc-dev-bin_2.31-0ubuntu9.2_amd64.deb ...
   Unpacking libc-dev-bin (2.31-0ubuntu9.2) ...
   Selecting previously unselected package linux-libc-dev:amd64.
   Preparing to unpack .../07-linux-libc-dev_5.4.0-96.109_amd64.deb ...
   Unpacking linux-libc-dev:amd64 (5.4.0-96.109) ...
   Selecting previously unselected package libcrypt-dev:amd64.
   Preparing to unpack .../08-libcrypt-dev_1%3a4.4.10-10ubuntu4_amd64.deb ...
   Unpacking libcrypt-dev:amd64 (1:4.4.10-10ubuntu4) ...
   Selecting previously unselected package libc6-dev:amd64.
   Preparing to unpack .../09-libc6-dev_2.31-0ubuntu9.2_amd64.deb ...
   Unpacking libc6-dev:amd64 (2.31-0ubuntu9.2) ...
   Selecting previously unselected package gcc-9-base:amd64.
   Preparing to unpack .../10-gcc-9-base_9.3.0-17ubuntu1~20.04_amd64.deb ...
   Unpacking gcc-9-base:amd64 (9.3.0-17ubuntu1~20.04) ...
   Selecting previously unselected package libisl22:amd64.
   Preparing to unpack .../11-libisl22_0.22.1-1_amd64.deb ...
   Unpacking libisl22:amd64 (0.22.1-1) ...
   Selecting previously unselected package libmpc3:amd64.
   Preparing to unpack .../12-libmpc3_1.1.0-1_amd64.deb ...
   Unpacking libmpc3:amd64 (1.1.0-1) ...
   Selecting previously unselected package cpp-9.
   Preparing to unpack .../13-cpp-9_9.3.0-17ubuntu1~20.04_amd64.deb ...
   Unpacking cpp-9 (9.3.0-17ubuntu1~20.04) ...
   Selecting previously unselected package cpp.
   Preparing to unpack .../14-cpp_4%3a9.3.0-1ubuntu2_amd64.deb ...
   Unpacking cpp (4:9.3.0-1ubuntu2) ...
   Selecting previously unselected package libcc1-0:amd64.
   Preparing to unpack .../15-libcc1-0_10.3.0-1ubuntu1~20.04_amd64.deb ...
   Unpacking libcc1-0:amd64 (10.3.0-1ubuntu1~20.04) ...
   Selecting previously unselected package libgomp1:amd64.
   Preparing to unpack .../16-libgomp1_10.3.0-1ubuntu1~20.04_amd64.deb ...
   Unpacking libgomp1:amd64 (10.3.0-1ubuntu1~20.04) ...
   Selecting previously unselected package libitm1:amd64.
   Preparing to unpack .../17-libitm1_10.3.0-1ubuntu1~20.04_amd64.deb ...
   Unpacking libitm1:amd64 (10.3.0-1ubuntu1~20.04) ...
   Selecting previously unselected package libatomic1:amd64.
   Preparing to unpack .../18-libatomic1_10.3.0-1ubuntu1~20.04_amd64.deb ...
   Unpacking libatomic1:amd64 (10.3.0-1ubuntu1~20.04) ...
   Selecting previously unselected package libasan5:amd64.
   Preparing to unpack .../19-libasan5_9.3.0-17ubuntu1~20.04_amd64.deb ...
   Unpacking libasan5:amd64 (9.3.0-17ubuntu1~20.04) ...
   Selecting previously unselected package liblsan0:amd64.
   Preparing to unpack .../20-liblsan0_10.3.0-1ubuntu1~20.04_amd64.deb ...
   Unpacking liblsan0:amd64 (10.3.0-1ubuntu1~20.04) ...
   Selecting previously unselected package libtsan0:amd64.
   Preparing to unpack .../21-libtsan0_10.3.0-1ubuntu1~20.04_amd64.deb ...
   Unpacking libtsan0:amd64 (10.3.0-1ubuntu1~20.04) ...
   Selecting previously unselected package libubsan1:amd64.
   Preparing to unpack .../22-libubsan1_10.3.0-1ubuntu1~20.04_amd64.deb ...
   Unpacking libubsan1:amd64 (10.3.0-1ubuntu1~20.04) ...
   Selecting previously unselected package libquadmath0:amd64.
   Preparing to unpack .../23-libquadmath0_10.3.0-1ubuntu1~20.04_amd64.deb ...
   Unpacking libquadmath0:amd64 (10.3.0-1ubuntu1~20.04) ...
   Selecting previously unselected package libgcc-9-dev:amd64.
   Preparing to unpack .../24-libgcc-9-dev_9.3.0-17ubuntu1~20.04_amd64.deb ...
   Unpacking libgcc-9-dev:amd64 (9.3.0-17ubuntu1~20.04) ...
   Selecting previously unselected package gcc-9.
   Preparing to unpack .../25-gcc-9_9.3.0-17ubuntu1~20.04_amd64.deb ...
   Unpacking gcc-9 (9.3.0-17ubuntu1~20.04) ...
   Selecting previously unselected package gcc.
   Preparing to unpack .../26-gcc_4%3a9.3.0-1ubuntu2_amd64.deb ...
   Unpacking gcc (4:9.3.0-1ubuntu2) ...
   Selecting previously unselected package libstdc++-9-dev:amd64.
   Preparing to unpack .../27-libstdc++-9-dev_9.3.0-17ubuntu1~20.04_amd64.deb ...
   Unpacking libstdc++-9-dev:amd64 (9.3.0-17ubuntu1~20.04) ...
   Selecting previously unselected package g++-9.
   Preparing to unpack .../28-g++-9_9.3.0-17ubuntu1~20.04_amd64.deb ...
   Unpacking g++-9 (9.3.0-17ubuntu1~20.04) ...
   Selecting previously unselected package g++.
   Preparing to unpack .../29-g++_4%3a9.3.0-1ubuntu2_amd64.deb ...
   Unpacking g++ (4:9.3.0-1ubuntu2) ...
   Selecting previously unselected package make.
   Preparing to unpack .../30-make_4.2.1-1.2_amd64.deb ...
   Unpacking make (4.2.1-1.2) ...
   Selecting previously unselected package libdpkg-perl.
   Preparing to unpack .../31-libdpkg-perl_1.19.7ubuntu3_all.deb ...
   Unpacking libdpkg-perl (1.19.7ubuntu3) ...
   Selecting previously unselected package dpkg-dev.
   Preparing to unpack .../32-dpkg-dev_1.19.7ubuntu3_all.deb ...
   Unpacking dpkg-dev (1.19.7ubuntu3) ...
   Selecting previously unselected package build-essential.
   Preparing to unpack .../33-build-essential_12.8ubuntu1.1_amd64.deb ...
   Unpacking build-essential (12.8ubuntu1.1) ...
   Selecting previously unselected package cmake-data.
   Preparing to unpack .../34-cmake-data_3.16.3-1ubuntu1_all.deb ...
   Unpacking cmake-data (3.16.3-1ubuntu1) ...
   Selecting previously unselected package libjsoncpp1:amd64.
   Preparing to unpack .../35-libjsoncpp1_1.7.4-3.1ubuntu2_amd64.deb ...
   Unpacking libjsoncpp1:amd64 (1.7.4-3.1ubuntu2) ...
   Selecting previously unselected package librhash0:amd64.
   Preparing to unpack .../36-librhash0_1.3.9-1_amd64.deb ...
   Unpacking librhash0:amd64 (1.3.9-1) ...
   Selecting previously unselected package cmake.
   Preparing to unpack .../37-cmake_3.16.3-1ubuntu1_amd64.deb ...
   Unpacking cmake (3.16.3-1ubuntu1) ...
   Selecting previously unselected package libfakeroot:amd64.
   Preparing to unpack .../38-libfakeroot_1.24-1_amd64.deb ...
   Unpacking libfakeroot:amd64 (1.24-1) ...
   Selecting previously unselected package fakeroot.
   Preparing to unpack .../39-fakeroot_1.24-1_amd64.deb ...
   Unpacking fakeroot (1.24-1) ...
   Selecting previously unselected package libalgorithm-diff-perl.
   Preparing to unpack .../40-libalgorithm-diff-perl_1.19.03-2_all.deb ...
   Unpacking libalgorithm-diff-perl (1.19.03-2) ...
   Selecting previously unselected package libalgorithm-diff-xs-perl.
   Preparing to unpack .../41-libalgorithm-diff-xs-perl_0.04-6_amd64.deb ...
   Unpacking libalgorithm-diff-xs-perl (0.04-6) ...
   Selecting previously unselected package libalgorithm-merge-perl.
   Preparing to unpack .../42-libalgorithm-merge-perl_0.08-3_all.deb ...
   Unpacking libalgorithm-merge-perl (0.08-3) ...
   Selecting previously unselected package libfile-fcntllock-perl.
   Preparing to unpack .../43-libfile-fcntllock-perl_0.22-3build4_amd64.deb ...
   Unpacking libfile-fcntllock-perl (0.22-3build4) ...
   Selecting previously unselected package manpages-dev.
   Preparing to unpack .../44-manpages-dev_5.05-1_all.deb ...
   Unpacking manpages-dev (5.05-1) ...
   Selecting previously unselected package ninja-build.
   Preparing to unpack .../45-ninja-build_1.10.0-1build1_amd64.deb ...
   Unpacking ninja-build (1.10.0-1build1) ...
   Selecting previously unselected package meson.
   Preparing to unpack .../46-meson_0.53.2-2ubuntu2_all.deb ...
   Unpacking meson (0.53.2-2ubuntu2) ...
   Selecting previously unselected package pkg-config.
   Preparing to unpack .../47-pkg-config_0.29.1-0ubuntu4_amd64.deb ...
   Unpacking pkg-config (0.29.1-0ubuntu4) ...
   Setting up manpages-dev (5.05-1) ...
   Setting up libfile-fcntllock-perl (0.22-3build4) ...
   Setting up libalgorithm-diff-perl (1.19.03-2) ...
   Setting up binutils-common:amd64 (2.34-6ubuntu1.3) ...
   Setting up linux-libc-dev:amd64 (5.4.0-96.109) ...
   Setting up libctf-nobfd0:amd64 (2.34-6ubuntu1.3) ...
   Setting up libgomp1:amd64 (10.3.0-1ubuntu1~20.04) ...
   Setting up libfakeroot:amd64 (1.24-1) ...
   Setting up ninja-build (1.10.0-1build1) ...
   Setting up fakeroot (1.24-1) ...
   update-alternatives: using /usr/bin/fakeroot-sysv to provide /usr/bin/fakeroot (fakeroot) in auto mode
   Setting up make (4.2.1-1.2) ...
   Setting up libquadmath0:amd64 (10.3.0-1ubuntu1~20.04) ...
   Setting up libmpc3:amd64 (1.1.0-1) ...
   Setting up libatomic1:amd64 (10.3.0-1ubuntu1~20.04) ...
   Setting up libdpkg-perl (1.19.7ubuntu3) ...
   Setting up libubsan1:amd64 (10.3.0-1ubuntu1~20.04) ...
   Setting up librhash0:amd64 (1.3.9-1) ...
   Setting up libcrypt-dev:amd64 (1:4.4.10-10ubuntu4) ...
   Setting up libisl22:amd64 (0.22.1-1) ...
   Setting up cmake-data (3.16.3-1ubuntu1) ...
   Setting up libbinutils:amd64 (2.34-6ubuntu1.3) ...
   Setting up libc-dev-bin (2.31-0ubuntu9.2) ...
   Setting up libalgorithm-diff-xs-perl (0.04-6) ...
   Setting up libcc1-0:amd64 (10.3.0-1ubuntu1~20.04) ...
   Setting up liblsan0:amd64 (10.3.0-1ubuntu1~20.04) ...
   Setting up libitm1:amd64 (10.3.0-1ubuntu1~20.04) ...
   Setting up libjsoncpp1:amd64 (1.7.4-3.1ubuntu2) ...
   Setting up gcc-9-base:amd64 (9.3.0-17ubuntu1~20.04) ...
   Setting up libalgorithm-merge-perl (0.08-3) ...
   Setting up libtsan0:amd64 (10.3.0-1ubuntu1~20.04) ...
   Setting up libctf0:amd64 (2.34-6ubuntu1.3) ...
   Setting up meson (0.53.2-2ubuntu2) ...
   Setting up libasan5:amd64 (9.3.0-17ubuntu1~20.04) ...
   Setting up pkg-config (0.29.1-0ubuntu4) ...
   Setting up cmake (3.16.3-1ubuntu1) ...
   Setting up cpp-9 (9.3.0-17ubuntu1~20.04) ...
   Setting up libc6-dev:amd64 (2.31-0ubuntu9.2) ...
   Setting up binutils-x86-64-linux-gnu (2.34-6ubuntu1.3) ...
   Setting up binutils (2.34-6ubuntu1.3) ...
   Setting up dpkg-dev (1.19.7ubuntu3) ...
   Setting up libgcc-9-dev:amd64 (9.3.0-17ubuntu1~20.04) ...
   Setting up cpp (4:9.3.0-1ubuntu2) ...
   Setting up gcc-9 (9.3.0-17ubuntu1~20.04) ...
   Setting up libstdc++-9-dev:amd64 (9.3.0-17ubuntu1~20.04) ...
   Setting up gcc (4:9.3.0-1ubuntu2) ...
   Setting up g++-9 (9.3.0-17ubuntu1~20.04) ...
   Setting up g++ (4:9.3.0-1ubuntu2) ...
   update-alternatives: using /usr/bin/g++ to provide /usr/bin/c++ (c++) in auto mode
   Setting up build-essential (12.8ubuntu1.1) ...
   Processing triggers for man-db (2.9.1-1) ...
   Processing triggers for libc-bin (2.31-0ubuntu9.2) ...

Of course, if you already have installed the packages before, apt will tell you.

Let's continue by installing the required libraries, namely
`Eigen <https://eigen.tuxfamily.org/>`__, `PCG <https://www.pcg-random.org/>`__,
`{fmt} <https://fmt.dev/>`__, and `Google Test <https://github.com/google/googletest>`__
using ``apt install``:

.. code-block:: bash

   $ sudo apt install libfmt-dev libeigen3-dev libpcg-cpp-dev googletest
   Reading package lists... Done
   Building dependency tree
   Reading state information... Done
   Suggested packages:
   libeigen3-doc libmpfrc++-dev libfmt-doc
   The following NEW packages will be installed:
   googletest libeigen3-dev libfmt-dev libpcg-cpp-dev
   0 upgraded, 4 newly installed, 0 to remove and 0 not upgraded.
   Need to get 1,617 kB of archives.
   After this operation, 12.6 MB of additional disk space will be used.
   Get:1 http://hr.archive.ubuntu.com/ubuntu focal/universe amd64 googletest all 1.10.0-2 [623 kB]
   Get:2 http://hr.archive.ubuntu.com/ubuntu focal/universe amd64 libpcg-cpp-dev all 0.98.1-2 [21.4 kB]
   Get:3 http://hr.archive.ubuntu.com/ubuntu focal/universe amd64 libeigen3-dev all 3.3.7-2 [815 kB]
   Get:4 http://hr.archive.ubuntu.com/ubuntu focal/universe amd64 libfmt-dev amd64 6.1.2+ds-2 [158 kB]
   Fetched 1,617 kB in 0s (3,660 kB/s)
   Selecting previously unselected package googletest.
   (Reading database ... 82001 files and directories currently installed.)
   Preparing to unpack .../googletest_1.10.0-2_all.deb ...
   Unpacking googletest (1.10.0-2) ...
   Selecting previously unselected package libpcg-cpp-dev.
   Preparing to unpack .../libpcg-cpp-dev_0.98.1-2_all.deb ...
   Unpacking libpcg-cpp-dev (0.98.1-2) ...
   Selecting previously unselected package libeigen3-dev.
   Preparing to unpack .../libeigen3-dev_3.3.7-2_all.deb ...
   Unpacking libeigen3-dev (3.3.7-2) ...
   Selecting previously unselected package libfmt-dev.
   Preparing to unpack .../libfmt-dev_6.1.2+ds-2_amd64.deb ...
   Unpacking libfmt-dev (6.1.2+ds-2) ...
   Setting up googletest (1.10.0-2) ...
   Setting up libeigen3-dev (3.3.7-2) ...
   Setting up libfmt-dev (6.1.2+ds-2) ...
   Setting up libpcg-cpp-dev (0.98.1-2) ...

If your output looks like this, you're all set to begin installing |Dock|. The
latest version of |Dock| can be found on
`Bitbucket <https://bitbucket.org/rxdock/rxdock>`__,
`GitHub <https://github.com/rxdock/rxdock>`__, and
`GitLab <https://gitlab.com/rxdock/rxdock>`__. The last one is the official
repository and the other two are just mirroring it. We'll be using GitLab in
the following steps.

.. code-block:: bash

   $ git clone https://gitlab.com/rxdock/rxdock.git
   Cloning into 'rxdock'...
   remote: Enumerating objects: 6181, done.
   remote: Counting objects: 100% (2505/2505), done.
   remote: Compressing objects: 100% (922/922), done.
   remote: Total 6181 (delta 1874), reused 2014 (delta 1564), pack-reused 3676
   Receiving objects: 100% (6181/6181), 33.53 MiB | 5.45 MiB/s, done.
   Resolving deltas: 100% (4610/4610), done.

You should have the directory called ``rxdock`` with all the necessary files
downloaded from GitLab. Let's enter that directory:

.. code-block:: bash

   $ cd rxdock

Here we will check out the `patched-rdock` tag which will give us what the name
says, a patched rDock. This is a good choice since documentation is still
written for rDock and does not cover the latest changes in the command-line
interface and configuration files.

.. note::

   You can choose to use the latest code (in the `master` branch) instead of
   the `patched-rdock` tag, but be aware that there are no guarantees about its
   correctness or performance: in other words, it might crash and/or produce
   incorrect results.

.. code-block:: bash

   $ git checkout patched-rdock
   Note: switching to 'patched-rdock'.

   You are in 'detached HEAD' state. You can look around, make experimental
   changes and commit them, and you can discard any commits you make in this
   state without impacting any branches by switching back to a branch.

   If you want to create a new branch to retain commits you create, you may
   do so (now or later) by using -c with the switch command. Example:

   git switch -c <new-branch-name>

   Or undo this operation with:

   git switch -

   Turn off this advice by setting config variable advice.detachedHead to false

   HEAD is now at c0c3a3a Fixed MSVC error C3016

|Dock|'s build configuration is done
:ref:`using the Meson build system <build-system>`. Here we set the build type
to ``release`` and enable the building of tests, which will be used to verify the
correctness of the resulting build.

.. code-block:: bash

   $ meson -Dbuildtype=release -Dtests=true builddir
   The Meson build system
   Version: 0.53.2
   Source dir: /home/pnikolic/rxdock
   Build dir: /home/pnikolic/rxdock/builddir
   Build type: native build
   Project name: RxDock
   Project version: 0.1.0
   C++ compiler for the host machine: c++ (gcc 9.3.0 "c++ (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0")
   C++ linker for the host machine: c++ ld.bfd 2.34
   Host machine cpu family: x86_64
   Host machine cpu: x86_64
   Found pkg-config: /usr/bin/pkg-config (0.29.1)
   Run-time dependency eigen3 found: YES 3.3.7
   Run-time dependency OpenMP found: YES 4.5
   Found CMake: /usr/bin/cmake (3.16.3)
   Run-time dependency pcg-cpp found: NO (tried pkgconfig and cmake)
   Looking for a fallback subproject for the dependency pcg-cpp
   Downloading pcg source from https://github.com/imneme/pcg-cpp/archive/v0.98.1.tar.gz
   Downloading file of unknown size.
   Downloading pcg patch from https://wrapdb.mesonbuild.com/v1/projects/pcg/0.98.1/1/get_zip
   Download size: 268
   Downloading: ..........

   |Executing subproject pcg method meson
   |
   |Project name: pcg-cpp
   |Project version: 0.98.1
   |C++ compiler for the host machine: c++ (gcc 9.3.0 "c++ (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0")
   |C++ linker for the host machine: c++ ld.bfd 2.34
   |Build targets in project: 0
   |Subproject pcg finished.

   Dependency pcg-cpp from subproject subprojects/pcg found: YES 0.98.1
   Run-time dependency cxxopts found: NO (tried pkgconfig and cmake)
   Looking for a fallback subproject for the dependency cxxopts
   Downloading cxxopts source from https://github.com/jarro2783/cxxopts/archive/v2.2.0.tar.gz
   Downloading file of unknown size.
   Downloading cxxopts patch from https://wrapdb.mesonbuild.com/v1/projects/cxxopts/2.2.0/1/get_zip
   Download size: 1022
   Downloading: ..........

   |Executing subproject cxxopts method meson
   |
   |Project name: cxxopts
   |Project version: 2.2.0
   |C++ compiler for the host machine: c++ (gcc 9.3.0 "c++ (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0")
   |C++ linker for the host machine: c++ ld.bfd 2.34
   |Build targets in project: 1
   |Subproject cxxopts finished.

   Dependency cxxopts from subproject subprojects/cxxopts found: YES 2.2.0
   Run-time dependency dt_smarts found: NO (tried pkgconfig and cmake)
   Run-time dependency dt_smiles found: NO (tried pkgconfig and cmake)
   Run-time dependency GTest found: YES (building self)
   Build targets in project: 9

   RxDock 0.1.0

   Subprojects
      cxxopts: YES
         pcg: YES

   Found ninja-1.10.0 at /usr/bin/ninja

Should Meson fail to find any of the dependencies listed above, use the
``apt search`` command to find the desired package and install it manually using
``apt install`` command. Also, should you fail to find the packages and get
stuck on this step, please use the
`issue tracker <https://gitlab.com/rxdock/rxdock/-/issues>`__ to check if a
similar issue to one that you have has already been reported and report it if
it has not.

Build |Dock| using the ``ninja`` command:

.. code-block::

   $ ninja -C builddir
   ninja: Entering directory `builddir'
   [9/154] Compiling C++ object 'rxdock@sha/src_GP_RbtGPParser.cxx.o'.
   ../src/GP/RbtGPParser.cxx: In member function 'RbtReturnType RbtGPParser::Parse1Output(RbtGPChromosomePtr, int)':
   ../src/GP/RbtGPParser.cxx:78:17: warning: unused variable 'f2' [-Wunused-variable]
      78 |   RbtReturnType f2 = commands[ncomm]->Execute();
         |                 ^~
   [31/154] Compiling C++ object 'rxdock@sha/src_lib_RbtBond.cxx.o'.
   ../src/lib/RbtBond.cxx: In constructor 'RbtBond::RbtBond(int, RbtAtomPtr&, RbtAtomPtr&, int)':
   ../src/lib/RbtBond.cxx:33:8: warning: unused variable 'bOK1' [-Wunused-variable]
      33 |   bool bOK1 = m_spAtom1->AddBond(this);
         |        ^~~~
   ../src/lib/RbtBond.cxx:34:8: warning: unused variable 'bOK2' [-Wunused-variable]
      34 |   bool bOK2 = m_spAtom2->AddBond(this);
         |        ^~~~
   ../src/lib/RbtBond.cxx: In destructor 'virtual RbtBond::~RbtBond()':
   ../src/lib/RbtBond.cxx:50:8: warning: unused variable 'bOK1' [-Wunused-variable]
      50 |   bool bOK1 = m_spAtom1->RemoveBond(this);
         |        ^~~~
   ../src/lib/RbtBond.cxx:51:8: warning: unused variable 'bOK2' [-Wunused-variable]
      51 |   bool bOK2 = m_spAtom2->RemoveBond(this);
         |        ^~~~
   [60/154] Compiling C++ object 'rxdock@sha/src_lib_RbtDockingSite.cxx.o'.
   ../src/lib/RbtDockingSite.cxx: In member function 'void RbtDockingSite::CreateGrid()':
   ../src/lib/RbtDockingSite.cxx:300:21: warning: comparison of integer expressions of different signedness: 'int' and 'unsigned int' [-Wsign-compare]
   300 |   for (int i = 0; i < m_spGrid->GetN(); i++) {
         |                   ~~^~~~~~~~~~~~~~~~~~
   [106/154] Compiling C++ object 'rxdock@sha/src_lib_RbtRealGrid.cxx.o'.
   In file included from /usr/include/eigen3/unsupported/Eigen/CXX11/Tensor:107,
                  from ../include/RbtRealGrid.h:24,
                  from ../src/lib/RbtRealGrid.cxx:18:
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h: In instantiation of 'struct Eigen::TensorEvaluator<const Eigen::Tensor<float, 3, 1>, Eigen::DefaultDevice>':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h:366:65:   recursively required from 'struct Eigen::TensorEvaluator<const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_cmp_op<float, float, Eigen::internal::cmp_GE>, const Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > >, Eigen::DefaultDevice>'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h:366:65:   required from 'struct Eigen::TensorEvaluator<const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_boolean_and_op, const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_cmp_op<float, float, Eigen::internal::cmp_GE>, const Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > >, const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_cmp_op<float, float, Eigen::internal::cmp_LT>, const Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > > >, Eigen::DefaultDevice>'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorAssign.h:96:70:   required from 'struct Eigen::TensorEvaluator<const Eigen::TensorAssignOp<Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_boolean_and_op, const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_cmp_op<float, float, Eigen::internal::cmp_GE>, const Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > >, const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_cmp_op<float, float, Eigen::internal::cmp_LT>, const Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > > > >, Eigen::DefaultDevice>'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:388:14:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::Tensor(const Eigen::TensorBase<OtherDerived, 0>&) [with OtherDerived = Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_boolean_and_op, const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_cmp_op<float, float, Eigen::internal::cmp_GE>, const Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > >, const Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_cmp_op<float, float, Eigen::internal::cmp_LT>, const Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > > >; Scalar_ = bool; int NumIndices_ = 3; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:207:41:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h:162:71: warning: ignoring attributes on template argument 'Eigen::PacketType<float, Eigen::DefaultDevice>::type' {aka '__vector(4) float'} [-Wignored-attributes]
   162 |     PacketAccess = (internal::unpacket_traits<PacketReturnType>::size > 1),
         |                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^~~~
   In file included from /usr/include/eigen3/unsupported/Eigen/CXX11/Tensor:107,
                  from ../include/RbtRealGrid.h:24,
                  from ../src/lib/RbtRealGrid.cxx:18:
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h: In instantiation of 'struct Eigen::TensorEvaluator<Eigen::Tensor<float, 3, 1>, Eigen::DefaultDevice>':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorAssign.h:100:65:   required from 'struct Eigen::TensorEvaluator<const Eigen::TensorAssignOp<Eigen::Tensor<float, 3, 1>, const Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> > >, Eigen::DefaultDevice>'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:416:14:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>& Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::operator=(const OtherDerived&) [with OtherDerived = Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> >; Scalar_ = float; int NumIndices_ = 3; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:208:56:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h:42:71: warning: ignoring attributes on template argument 'Eigen::PacketType<float, Eigen::DefaultDevice>::type' {aka '__vector(4) float'} [-Wignored-attributes]
      42 |     PacketAccess = (internal::unpacket_traits<PacketReturnType>::size > 1),
         |                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^~~~
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h: In instantiation of 'struct Eigen::TensorEvaluator<Eigen::Tensor<float, 0, 1>, Eigen::DefaultDevice>':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorAssign.h:100:65:   required from 'struct Eigen::TensorEvaluator<const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >, Eigen::DefaultDevice>'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:388:14:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::Tensor(const Eigen::TensorBase<OtherDerived, 0>&) [with OtherDerived = Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Scalar_ = float; int NumIndices_ = 0; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:352:70:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h:42:71: warning: ignoring attributes on template argument 'Eigen::PacketType<float, Eigen::DefaultDevice>::type' {aka '__vector(4) float'} [-Wignored-attributes]
   In file included from /usr/include/eigen3/unsupported/Eigen/CXX11/Tensor:141,
                  from ../include/RbtRealGrid.h:24,
                  from ../src/lib/RbtRealGrid.cxx:18:
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h: In instantiation of 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 3, 1>, const Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> > >]':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:417:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>& Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::operator=(const OtherDerived&) [with OtherDerived = Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> >; Scalar_ = float; int NumIndices_ = 3; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:208:56:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:61:17: warning: ignoring attributes on template argument 'Eigen::TensorEvaluator<const Eigen::TensorAssignOp<Eigen::Tensor<float, 3, 1>, const Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> > >, Eigen::DefaultDevice>::PacketReturnType' {aka '__vector(4) float'} [-Wignored-attributes]
      61 |       const int PacketSize = unpacket_traits<typename TensorEvaluator<Expression, DefaultDevice>::PacketReturnType>::size;
         |                 ^~~~~~~~~~
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h: In instantiation of 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >]':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:389:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::Tensor(const Eigen::TensorBase<OtherDerived, 0>&) [with OtherDerived = Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Scalar_ = float; int NumIndices_ = 0; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:352:70:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:61:17: warning: ignoring attributes on template argument 'Eigen::TensorEvaluator<const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >, Eigen::DefaultDevice>::PacketReturnType' {aka '__vector(4) float'} [-Wignored-attributes]
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h: In instantiation of 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >]':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:389:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::Tensor(const Eigen::TensorBase<OtherDerived, 0>&) [with OtherDerived = Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Scalar_ = float; int NumIndices_ = 0; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:357:70:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:61:17: warning: ignoring attributes on template argument 'Eigen::TensorEvaluator<const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >, Eigen::DefaultDevice>::PacketReturnType' {aka '__vector(4) float'} [-Wignored-attributes]
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h: In instantiation of 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 3, 1>, const Eigen::Tensor<float, 3, 1> >]':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:407:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>& Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::operator=(const Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>&) [with Scalar_ = float; int NumIndices_ = 3; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:520:63:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:61:17: warning: ignoring attributes on template argument 'Eigen::TensorEvaluator<const Eigen::TensorAssignOp<Eigen::Tensor<float, 3, 1>, const Eigen::Tensor<float, 3, 1> >, Eigen::DefaultDevice>::PacketReturnType' {aka '__vector(4) float'} [-Wignored-attributes]
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h: In instantiation of 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > >]':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:417:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>& Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::operator=(const OtherDerived&) [with OtherDerived = Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >; Scalar_ = float; int NumIndices_ = 3; int Options_ = 1; IndexType_ = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorBase.h:847:24:   required from 'Derived& Eigen::TensorBase<Derived, AccessLevel>::setConstant(const Scalar&) [with Derived = Eigen::Tensor<float, 3, 1>; int AccessLevel = 1; Eigen::TensorBase<Derived, AccessLevel>::Scalar = float]'
   ../src/lib/RbtRealGrid.cxx:200:68:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:61:17: warning: ignoring attributes on template argument 'Eigen::TensorEvaluator<const Eigen::TensorAssignOp<Eigen::Tensor<float, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> > >, Eigen::DefaultDevice>::PacketReturnType' {aka '__vector(4) float'} [-Wignored-attributes]
   In file included from /usr/include/eigen3/unsupported/Eigen/CXX11/Tensor:107,
                  from ../include/RbtRealGrid.h:24,
                  from ../src/lib/RbtRealGrid.cxx:18:
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h: In instantiation of 'const int Eigen::TensorEvaluator<const Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> >, Eigen::DefaultDevice>::PacketSize':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h:600:36:   required from 'Eigen::TensorEvaluator<const Eigen::TensorSelectOp<IfArgType, ThenArgType, ElseArgType>, Device>::PacketReturnType Eigen::TensorEvaluator<const Eigen::TensorSelectOp<IfArgType, ThenArgType, ElseArgType>, Device>::packet(Eigen::TensorEvaluator<const Eigen::TensorSelectOp<IfArgType, ThenArgType, ElseArgType>, Device>::Index) const [with int LoadMode = 16; IfArgType = const Eigen::Tensor<bool, 3, 1>; ThenArgType = const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >; ElseArgType = const Eigen::Tensor<float, 3, 1>; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorSelectOp<IfArgType, ThenArgType, ElseArgType>, Device>::PacketReturnType = __vector(4) float; Eigen::TensorEvaluator<const Eigen::TensorSelectOp<IfArgType, ThenArgType, ElseArgType>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorAssign.h:141:5:   required from 'void Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::evalPacket(Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index) [with LeftArgType = Eigen::Tensor<float, 3, 1>; RightArgType = const Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> >; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:68:11:   required from 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 3, 1>, const Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> > >]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:417:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>& Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::operator=(const OtherDerived&) [with OtherDerived = Eigen::TensorSelectOp<const Eigen::Tensor<bool, 3, 1>, const Eigen::TensorCwiseNullaryOp<Eigen::internal::scalar_constant_op<float>, const Eigen::Tensor<float, 3, 1> >, const Eigen::Tensor<float, 3, 1> >; Scalar_ = float; int NumIndices_ = 3; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:208:56:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorEvaluator.h:572:20: warning: ignoring attributes on template argument 'Eigen::PacketType<float, Eigen::DefaultDevice>::type' {aka '__vector(4) float'} [-Wignored-attributes]
   572 |   static const int PacketSize = internal::unpacket_traits<PacketReturnType>::size;
         |                    ^~~~~~~~~~
   In file included from /usr/include/eigen3/unsupported/Eigen/CXX11/Tensor:109,
                  from ../include/RbtRealGrid.h:24,
                  from ../src/lib/RbtRealGrid.cxx:18:
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h: In instantiation of 'const int Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>, Eigen::DefaultDevice>::PacketSize':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h:613:5:   required from 'Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::PacketReturnType Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::packet(Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::Index) const [with int LoadMode = 0; Op = Eigen::internal::MinReducer<float>; Dims = const Eigen::DimensionList<long int, 3>; ArgType = const Eigen::Tensor<float, 3, 1>; MakePointer_ = Eigen::MakePointer; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::PacketReturnType = __vector(4) float; Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorAssign.h:141:5:   required from 'void Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::evalPacket(Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index) [with LeftArgType = Eigen::Tensor<float, 0, 1>; RightArgType = const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:68:11:   required from 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:389:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::Tensor(const Eigen::TensorBase<OtherDerived, 0>&) [with OtherDerived = Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Scalar_ = float; int NumIndices_ = 0; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:352:70:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h:397:20: warning: ignoring attributes on template argument 'Eigen::PacketType<float, Eigen::DefaultDevice>::type' {aka '__vector(4) float'} [-Wignored-attributes]
   397 |   static const int PacketSize = internal::unpacket_traits<PacketReturnType>::size;
         |                    ^~~~~~~~~~
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h: In instantiation of 'const int Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>, Eigen::DefaultDevice>::PacketSize':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h:613:5:   required from 'Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::PacketReturnType Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::packet(Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::Index) const [with int LoadMode = 0; Op = Eigen::internal::MaxReducer<float>; Dims = const Eigen::DimensionList<long int, 3>; ArgType = const Eigen::Tensor<float, 3, 1>; MakePointer_ = Eigen::MakePointer; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::PacketReturnType = __vector(4) float; Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorAssign.h:141:5:   required from 'void Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::evalPacket(Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index) [with LeftArgType = Eigen::Tensor<float, 0, 1>; RightArgType = const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:68:11:   required from 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:389:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::Tensor(const Eigen::TensorBase<OtherDerived, 0>&) [with OtherDerived = Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Scalar_ = float; int NumIndices_ = 0; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:357:70:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h:397:20: warning: ignoring attributes on template argument 'Eigen::PacketType<float, Eigen::DefaultDevice>::type' {aka '__vector(4) float'} [-Wignored-attributes]
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h: In instantiation of 'static typename Self::CoeffReturnType Eigen::internal::InnerMostDimReducer<Self, Op, true>::reduce(const Self&, typename Self::Index, typename Self::Index, Op&) [with Self = Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>, Eigen::DefaultDevice>; Op = Eigen::internal::MinReducer<float>; typename Self::CoeffReturnType = float; typename Self::Index = long int]':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h:627:68:   required from 'Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::PacketReturnType Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::packet(Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::Index) const [with int LoadMode = 0; Op = Eigen::internal::MinReducer<float>; Dims = const Eigen::DimensionList<long int, 3>; ArgType = const Eigen::Tensor<float, 3, 1>; MakePointer_ = Eigen::MakePointer; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::PacketReturnType = __vector(4) float; Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorAssign.h:141:5:   required from 'void Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::evalPacket(Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index) [with LeftArgType = Eigen::Tensor<float, 0, 1>; RightArgType = const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:68:11:   required from 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:389:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::Tensor(const Eigen::TensorBase<OtherDerived, 0>&) [with OtherDerived = Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Scalar_ = float; int NumIndices_ = 0; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:352:70:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h:169:15: warning: ignoring attributes on template argument 'Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Eigen::internal::MinReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>, Eigen::DefaultDevice>::PacketReturnType' {aka '__vector(4) float'} [-Wignored-attributes]
   169 |     const int packetSize = internal::unpacket_traits<typename Self::PacketReturnType>::size;
         |               ^~~~~~~~~~
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h: In instantiation of 'static typename Self::CoeffReturnType Eigen::internal::InnerMostDimReducer<Self, Op, true>::reduce(const Self&, typename Self::Index, typename Self::Index, Op&) [with Self = Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>, Eigen::DefaultDevice>; Op = Eigen::internal::MaxReducer<float>; typename Self::CoeffReturnType = float; typename Self::Index = long int]':
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h:627:68:   required from 'Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::PacketReturnType Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::packet(Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::Index) const [with int LoadMode = 0; Op = Eigen::internal::MaxReducer<float>; Dims = const Eigen::DimensionList<long int, 3>; ArgType = const Eigen::Tensor<float, 3, 1>; MakePointer_ = Eigen::MakePointer; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::PacketReturnType = __vector(4) float; Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Op, Dims, XprType, MakePointer_>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorAssign.h:141:5:   required from 'void Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::evalPacket(Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index) [with LeftArgType = Eigen::Tensor<float, 0, 1>; RightArgType = const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Device = Eigen::DefaultDevice; Eigen::TensorEvaluator<const Eigen::TensorAssignOp<LhsXprType, RhsXprType>, Device>::Index = long int]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorExecutor.h:68:11:   required from 'static void Eigen::internal::TensorExecutor<Expression, Eigen::DefaultDevice, true>::run(const Expression&, const Eigen::DefaultDevice&) [with Expression = const Eigen::TensorAssignOp<Eigen::Tensor<float, 0, 1>, const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer> >]'
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/Tensor.h:389:65:   required from 'Eigen::Tensor<Scalar_, NumIndices_, Options_, IndexType>::Tensor(const Eigen::TensorBase<OtherDerived, 0>&) [with OtherDerived = Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>; Scalar_ = float; int NumIndices_ = 0; int Options_ = 1; IndexType_ = long int]'
   ../src/lib/RbtRealGrid.cxx:357:70:   required from here
   /usr/include/eigen3/unsupported/Eigen/CXX11/src/Tensor/TensorReduction.h:169:15: warning: ignoring attributes on template argument 'Eigen::TensorEvaluator<const Eigen::TensorReductionOp<Eigen::internal::MaxReducer<float>, const Eigen::DimensionList<long int, 3>, const Eigen::Tensor<float, 3, 1>, Eigen::MakePointer>, Eigen::DefaultDevice>::PacketReturnType' {aka '__vector(4) float'} [-Wignored-attributes]
   [130/154] Compiling C++ object 'rxdock@sha/src_lib_RbtVdwSF.cxx.o'.
   ../src/lib/RbtVdwSF.cxx: In member function 'void RbtVdwSF::Setup()':
   ../src/lib/RbtVdwSF.cxx:364:36: warning: 'alphai' may be used uninitialized in this function [-Wmaybe-uninitialized]
   364 |         double D = 0.345 * Ii * Ij * alphai * alphaj / (Ii + Ij);
         |                    ~~~~~~~~~~~~~~~~^~~~~~~~
   ../src/lib/RbtVdwSF.cxx:364:60: warning: 'Ii' may be used uninitialized in this function [-Wmaybe-uninitialized]
   364 |         double D = 0.345 * Ii * Ij * alphai * alphaj / (Ii + Ij);
         |                                                        ~~~~^~~~~
   [154/154] Linking target unit-test.

Run the tests, also using the ``ninja`` command, but this time with the argument
``test``:

.. code-block:: bash

   $ ninja -C builddir test
   ninja: Entering directory 'builddir'
   [0/1] Running all tests.
   1/5 rbcavity-1koc                           OK       0.87 s
   2/5 rbcavity-1yet                           OK       6.68 s
   3/5 rbcavity-1yet-test                      OK       6.29 s
   4/5 unit-test                               OK      11.84 s
   5/5 rbdock-1yet-test                        OK       1.77 s

   Ok:                    5
   Expected Fail:         0
   Fail:                  0
   Unexpected Pass:       0
   Skipped:               0
   Timeout:               0

   Full log written to /home/pnikolic/rxdock/builddir/meson-logs/testlog.txt

If everything went fine, your system is up to date, and you followed these
instructions to the letter -- congratulations! |Dock| is now successfully built
and tested on your desktop or laptop computer, and you can start using it
immediately.

To check where your |Dock| is built, go to your ``rxdock`` directory and run
``pwd``:

.. code-block:: bash

   $ pwd
   /home/pnikolic/rxdock

In order to start using |Dock|, we have to set up a working environment. Go
back to your home directory in the current terminal or open up a new terminal.

.. code-block:: bash

   $ cd

Make a new directory called ``training``, and enter it.

.. code-block:: bash

   $ mkdir training
   $ cd training

Remember the path we obtained via the ``pwd`` command? Now is the time to use
it. Execute the following instructions and modify the ``RBT_ROOT`` path to
correspond to that path, in the same manner, the code below corresponds to my
path in the shown above.

.. code-block:: bash

   $ export RBT_ROOT=/home/pnikolic/rxdock
   $ export PATH=$RBT_ROOT/bin:$PATH
   $ export PERL5LIB=$RBT_ROOT/lib:$PERL5LIB
   $ export PATH=$RBT_ROOT/builddir:$PATH
   $ export LD_LIBRARY_PATH=$RBT_ROOT/builddir:$LD_LIBRARY_PATH

Let us verify that the |Dock| is properly installed and the working environment
is properly set up by running ``rbcavity`` and ``rbdock`` commands with the
parameter ``--help`` or ``-h``:

.. code-block:: bash

   $ rbcavity -h
   ***********************************************
   The RxDock molecular docking program is licensed under GNU LGPL version 3.
   RxDock is maintained by Vedran Miletić, Patrik Nikolić, and Luka Vretenar.
   Visit https://rxdock.gitlab.io/ for more information.
   Executable:	rbcavity/0.1.0
   Library:	librxdock.so/0.1.0
   RBT_ROOT:	/home/pnikolic/rxdock
   RBT_HOME:	/home/pnikolic
   Current dir:	/home/pnikolic/rxdock
   Date:		Thu Jan 27 22:02:56 2022
   ***********************************************
   rbcavity - calculate docking cavities
   Usage:
      rbcavity [OPTION...]

      -r, --receptor-param arg      receptor param file (contains active site
                                    params)
      -W, --write-docking-cavities  write docking cavities (plus distance grid)
                                    to .as file
      -R, --read-docking-cavities   read docking cavities (plus distance grid)
                                    from .as file
      -d, --write-insightii-grids   dump InsightII grids for each cavity for
                                    visualisation
      -v, --write-psf-crd           dump target PSF/CRD files for rDock Viewer
      -l, --list-atoms-dist arg     list receptor atoms within specified distance
                                    of any cavity (in angstrom) (default: 5.0)
      -s, --print-site              print SITE descriptors (counts of exposed
                                    atoms)
      -b, --border arg              set the border around the cavities for the
                                    distance grid (in angstrom) (default: 8.0)
      -m, --write-moe-grid          write active site into a MOE grid
      -h, --help                    Print help

.. code-block:: bash

   $ rbdock -h
   ***********************************************
   The RxDock molecular docking program is licensed under GNU LGPL version 3.
   RxDock is maintained by Vedran Miletić, Patrik Nikolić, and Luka Vretenar.
   Visit https://rxdock.gitlab.io/ for more information.
   Executable:	rbdock/0.1.0
   Library:	librxdock.so/0.1.0
   RBT_ROOT:	/home/pnikolic/rxdock
   RBT_HOME:	/home/pnikolic
   Current dir:	/home/pnikolic/rxdock
   Date:		Thu Jan 27 22:03:11 2022
   ***********************************************
   rbdock - docking engine
   Usage:
      rbdock [OPTION...]

      -i, --input arg           input ligand SD file
      -o, --output arg          output file name(s) prefix
      -r, --receptor-param arg  receptor parameter file
      -p, --docking-param arg   docking protocol parameter file
      -n, --number arg          number of runs per ligand (0 = unlimited)
                                 (default: 0)
      -P, --protonate           protonate all neutral amines, guanidines,
                                 imidazoles
      -D, --deprotonate         deprotonate all carboxylic, sulphur and
                                 phosphorous acid groups
      -H, --all-hydrogens       read all hydrogens present instead of only polar
                                 hydrogens
      -t, --threshold arg       score threshold
      -C, --continue            continue if score threshold is met instead of
                                 terminating ligand
      -f, --filter arg          filter file name
      -T, --trace arg           controls output level for debugging (0 = minimal,
                                 >0 = more verbose) (default: 0)
      -s, --seed arg            random number seed to use instead of
                                 std::random_device
      -h, --help                Print help

There are :ref:`other programs <programs>` as well. You can follow the
:ref:`Docking in 3 steps <docking-in-3-steps>` part of the documentation from
here on. Just be careful with ``-was`` and ``-ras`` parameters for ``rbcavity``
which were removed and replaced with ``-W`` and ``-R``, respectively. Don't
worry too much, if you forget and use ``-was`` or ``-ras``, ``rbcavity`` will
warn you about the removal and suggest what parameters to use instead.

And that's it! You have your |Dock| properly built and your workspace is
properly set up. Just remember that you have to set up your working environment
every time you start working. That means every time you reboot your computer,
log out and then back in, close the terminal and open it again, etc. you have
to write those five lines of code above. To avoid that, you can put them into
your shell configuration (``.bashrc`` file in your home directory, in my case
``/home/pnikolic``).

Should you have any questions, do not hesitate to reach out to the |Dock|
community on the `issue tracker <https://gitlab.com/rxdock/rxdock/-/issues>`__,
but please search for similar issues first.

**Happy docking!**
