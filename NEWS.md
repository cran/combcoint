# combcoint 0.3.0

## Internal changes

- Removed the dependency on `tsDyn`; the Johansen ML VECM estimation is now
  implemented internally instead of calling `tsDyn::VECM()`.

# combcoint 0.2.0

## Internal changes
- Updated the internal model used for predicting crit values.

# combcoint 0.1.0

## First release to CRAN 

- Implements the Bayer-Hanck cointegration test.
- Includes Engle-Granger and Johansen-style alternatives for comparison.
- Compatible with `formula` interface and supports custom lag settings.
- Includes example datasets and reproducible examples in documentation.
