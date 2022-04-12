# Examples

## Sunspots and Melanoma

The data can be obtained through the following helper function

```@docs
get_sunspot_melanoma_data
```

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

Estimate the model parameters using `mcar`

```@example sunspot_melanoma
y = [melanoma sunspots]
model, _, _ = mcar(y)
```

Perform a Granger causality test to get a connectivity matrix

```@example sunspot_melanoma
granger_causality_test(y, model);
```

Original PDC analysis

```@example sunspot_melanoma
original_pdc = pdc(model, y; metric = "euc")
```
