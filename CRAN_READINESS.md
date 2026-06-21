# CRAN Readiness

This file records the remaining tasks before submitting `rfedqr` to CRAN.
It is excluded from package builds by `.Rbuildignore`.

## Required before submission

1. Confirm that the maintainer identity in `DESCRIPTION` is final.
2. Confirm that the copyright holder in `LICENSE` is final.
3. Rebuild the package:
   ```sh
   R CMD build .
   ```
4. Run the CRAN-style check:
   ```sh
   R CMD check --as-cran rfedqr_0.1.0.tar.gz
   ```
5. Confirm that the result is clean enough for submission.
6. Submit through the official CRAN web form:
   https://cran.r-project.org/submit.html

## Current local status

The package passes CRAN-style local checks with no errors and no warnings.
The remaining local notes are:

1. New submission.
2. GitHub URL check timed out during the incoming check, although the URL was
   verified separately and returned HTTP status 200.
3. HTML validation was skipped because the local HTML Tidy installation is not
   recent enough.
