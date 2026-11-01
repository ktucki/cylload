# Contributing

Bug reports and pull requests are welcome. A contribution should:

* build without warnings with `-Wall -Wextra -Wpedantic -Wshadow` (GCC, Clang) and `/W4`
  (MSVC), using C++17 and no external dependencies;
* keep `ctest` green and add a test for new behaviour. A closed-form or published
  reference is preferred over a regression of the program's own output;
* keep the units of the public interfaces (mm, kg, N, bar, kg·m²) and document new
  input keywords in `docs/input_format.md`.

To report a problem with a STEP file, attach the output of
`cylload step FILE --level 4` (the `volume_check` and `moment_check` columns show
faces that are not integrated consistently) or a minimal file that reproduces it.
