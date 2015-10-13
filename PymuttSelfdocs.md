
```
----- pymutt usage -----
Estimate the spectral power (and auxiliary quantities, see below)
of a one-dimenstional time series using Thomson's multi-taper algorithm.

result_dict = mtft(
   array data,      a numpy array of data
   dt = 1.0,
   kind = 2,        1 for hires, 2 for adaptive
   nwin = 5,
   npi = 3.0,
   paddenlen = 0,   if > length(data), zero-pad to this length
   dodof = 0,       produce dof and F arrays
   doweights = 0,   returns all of the taper weights
   lines = None,    desired frequencies of spectral lines at which
                    to estimate complex spectral amplitudes.
                    Produces lamp and reshaped spectrum outputs.
                    Implies dodof = 1 as well.  Must be a numpy array.
    )

mtft returns a dictionary containing both fixed and optional items.
copies of input values:
  'dt'        time series sample interval
  'nin'       number of points in original time series.
  'npi'       dimensionless width of the spectral averaging window
  'nwin'      maximum number of tapers to combine
  'kind'      spectral algorithm: 1 for hires, 2 for adaptive
  'linef'     frequencies at which to estimate spectral line amplitudes
              by default this array is empty
fixed output values:
  'df'        frequency-domain spacing of spectral estimates
  'n'         number of points in the padded input array
  'nspec'     number of points in the spectrum
  'W'         dimensionless averaging width in terms of the
              possibly-padded input series
  'power'     estimated power (not amplitude) spectrum
optional output values:
  'dof'       degrees-of-freedom vs frequency
  'F'         F-test for spectral line vs frequency
  'weights'   nwin-by-nf array of spectral weights
  'linea'     estimated line amplitudes (at frequencies in linef)
  'linevar'   estimated line amplitude variance
  'reshaped'  spectrum with lines removed and reshaped around holes
```