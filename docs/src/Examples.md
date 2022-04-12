# Examples

## Sunspots and Melanoma

The data can be obtained through the following helper function

```@docs
get_sunspot_melanoma_data
```

### Setup

A first look at the detrended data

```@setup sunspot_melanoma
using Plots, AsymptoticPDC

```

```@example sunspot_melanoma
u = get_sunspot_melanoma_data()
t = u[:, 1]
melanoma, _ = detrend(u[:, 3])
sunspots, _ = detrend(u[:, 4])
p1 = plot(t, melanoma, lab = "Detrended annual melanoma incidence Connecticut")
p2 = plot(t, sunspots, lab = "Detrended annual sunpot number")
plot(p1, p2, layout = (2, 1))
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
println(connectivity) # hide
println("Granger causality test p-values:") # hide
println(pvalues) # hide
```

Original PDC analysis

```@example sunspot_melanoma
original_pdc = pdc(model, y; metric = "euc")
pdcplot(original_pdc, cnames = ["Melanoma", "Sunspots"])
```

Generalized PDC analysis

```@example sunspot_melanoma
generalized_pdc = pdc(model, y; metric = "diag")
pdcplot(generalized_pdc, cnames = ["Melanoma", "Sunspots"])
```