Coding standards
================

|Dock| is written using the `LLVM style C++ code
<https://clang.llvm.org/docs/ClangFormatStyleOptions.html>`__ and formatted
using `ClangFormat <https://clang.llvm.org/docs/ClangFormat.html>`__.

LLVM coding standards and additions
-----------------------------------

The coding standards are documented in `LLVM Coding Standards
<https://llvm.org/docs/CodingStandards.html>`__. With regards to that document,
|Dock| has made the following changes in the coding standards:

* Standard C++11 is used instead of C++14. In the near future the switch will
  be made to C++17.
* Implementation file extension should be ``.cxx`` instead of ``.cpp``.
* Variable names should use ``camelBack`` instead of ``CamelCase`` and member
  variables should be prefixed with ``m_`` (similar to the proposed `Variable
  Names Plan <https://llvm.org/docs/Proposals/VariableNames.html>`__ in LLVM).

See also
--------

Helpful online sources of knowledge about C++:

* `News, Status and Discussion about Standard C++ <https://isocpp.org/>`__
* `C++ reference <https://cppreference.com/>`__
* `C++ subreddit (/r/cpp) <https://www.reddit.com/r/cpp/>`__

Particularly useful books about software design, algorithms, and programming:

#. `Code Complete
   <https://www.microsoftpressstore.com/store/code-complete-9780735619678>`__
   by `Steve McConnel <https://stevemcconnell.com/>`__, published by Microsoft
   Press in 2004.

#. `Metaheuristics: From Design to Implementation
   <https://www.wiley.com/en-us/Metaheuristics%3A+From+Design+to+Implementation+-p-9780470278581>`_
   by `El-Ghazali Talbi <http://www.lifl.fr/~talbi/>`__, published by Wiley in 2009.

#. `Effective Modern C++
   <https://www.oreilly.com/library/view/effective-modern-c/9781491908419/>`_
   by `Scott Meyers <https://www.aristeia.com/>`__, published by O'Reilly
   in 2014. Also interesting and useful are "Effective C++", "More Effective
   C++", and "Effective STL" by the same author.

The future development will switch from the use of C++11 to the use of C++17. A
very good book covering the new features in C++17 is `C++17 in Detail
<https://leanpub.com/cpp17indetail>`__ by `Bartlomiej Filipek
<https://www.bfilipek.com/>`__, published by Leanpub in 2017. The blog post
`C++17 Features <https://www.bfilipek.com/2017/01/cpp17features.html>`__ is a
shorter version.
