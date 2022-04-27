# Examples

Some examples for demonstration purposes.

## Sunspots and Melanoma

 The etiology of melanoma is complex and may include the influences of trauma, heredity and hormonal activity. In particular, exposure to solar radiation may be involved in the pathogenesis of melanoma. Melanoma is more common in fair-skinned individuals and most frequent in skin sites exposed to the sun. In white populations melanoma is more common in areas closer to the equator where the intensity of solar radiation is higher. Data from various parts of the world suggest that the incidence of melanoma is increasing. The data below, giving age-adjusted melanoma incidence, are from the Connecticut Tumor Registry from 1936-1972. Connecticut has the longest record of state population-based cancer statistics in the United States of America. The data also includes the sunspot relative number. Houghton, Munster and Viola (1978) have shown that the age-adjusted incidence rate for malignant melanoma in the state of Connecticut has risen since 1935 and that superimposed on the rise are 3-5 year periods in which the rise in the rate of incidence is excessive. These periods have a cycle of 8-11 years and follow times of maximum sunspot activity. The relationship between solar cycles and melanoma supports the hypothesis that melanoma is related to sun exposure and provides evidence that solar radiation may trigger the development of clinically apparent melanoma.

The data can be obtained through the following helper function

```@docs
get_sunspot_melanoma_data
```

### Setup

A first look at the detrended data

```@setup sunspot_melanoma
using Plots, AsymptoticPDC
default(dpi = 200)
```

```@example sunspot_melanoma
u = get_sunspot_melanoma_data()
t = u[:, 1]
melanoma, _ = detrend(u[:, 3])
sunspots, _ = detrend(u[:, 4])
p1 = plot(t, melanoma, title = "Detrended annual melanoma incidence Connecticut")
p2 = plot(t, sunspots, title = "Detrended annual sunspot number")
plot(p1, p2, layout = (2, 1), legend = :none)
```

### Estimation

Estimate the model parameters using `mcar`

```@example sunspot_melanoma
y = [melanoma sunspots]
model, _, _ = mcar(y)
```

### Connectivity

Perform a Granger causality test to get a connectivity matrix

```@example sunspot_melanoma
connectivity, pvalues = granger_causality_test(model, y);
println("Connectivity matrix:") # hide
show(stdout, "text/plain", connectivity) # hide
println("Granger causality test p-values:") # hide
show(stdout, "text/plain", pvalues) # hide
```

### PDC Analysis

#### Original PDC analysis

The time series are not standardized and the sunspots series has a larger variance which explains the counterintuitive results below

```@example sunspot_melanoma
original_pdc = pdc(model, y; metric = "euc")
pdcplot(original_pdc, cnames = ["Melanoma", "Sunspots"])
```

But when calculating the asymptotic statistics for the original PDC analysis, only the coherence from sunspots -> melanoma is considered significant (grey shaded area). 

```@example sunspot_melanoma
original_pdc_asymp = pdc(model, y; metric = "euc", α = 0.01)
pdcplot(original_pdc_asymp, cnames = ["Melanoma", "Sunspots"])
```

#### Generalized PDC analysis

The generalized PDC analysis gives a better view inside the "true" interaction. The ~11 year sun activity cycle is clearly visible

```@example sunspot_melanoma
generalized_pdc = pdc(model, y; metric = "diag", α = 0.01)
pdcplot(generalized_pdc, cnames = ["Melanoma", "Sunspots"])
```

#### Information PDC analysis

The information PDC analysis gives similar results

```@example sunspot_melanoma
info_pdc = pdc(model, y; metric = "info", α = 0.01)
pdcplot(info_pdc, cnames = ["Melanoma", "Sunspots"])
```
