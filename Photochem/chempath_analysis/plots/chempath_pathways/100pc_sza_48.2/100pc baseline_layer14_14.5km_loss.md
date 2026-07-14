# 100pc baseline | O3 loss | 14.5 km

## 1. 2.235e-06%
    O3+hv -> O+O2
    O -> O_transport
    Net:O3 -> O2 + O_transport

## 2. 6.781e-08%
    O1D+N2 -> O+N2
    O3+hv -> O1D+O2
    O -> O_transport
    Net:O3 -> O2 + O_transport

## 3. 3.157e-08%
    O3+hv -> O1D+O2
    O1D -> O1D_transport
    Net:O3 -> O2 + O1D_transport

## 4. 2.678e-08%
    O1D+O2 -> O+O2
    O3+hv -> O1D+O2
    O -> O_transport
    Net:O3 -> O2 + O_transport

## 5. 3.827e-09%
    NO+O3 -> NO2+O2
    NO2+hv -> NO+O
    O -> O_transport
    Net:O3 -> O2 + O_transport

## 6. 1.052e-10%
    O1D+CO2 -> CO2+O
    O3+hv -> O1D+O2
    O -> O_transport
    Net:O3 -> O2 + O_transport

## 7. 8.862e-11%
    OH+O3 -> HO2+O2
    NO+HO2 -> NO2+OH
    NO2+hv -> NO+O
    O -> O_transport
    Net:O3 -> O2 + O_transport

## 8. 9.787e-12%
    O1D+H2O -> OH+OH
    O3+hv -> O1D+O2
    2(OH -> OH_transport)
    Net:H2O + O3 -> O2 + 2OH_transport

## 9. 5.932e-12%
    HO2+HO2 -> H2O2+O2
    2(OH+O3 -> HO2+O2)
    H2O2+hv -> OH+OH
    Net:2O3 -> 3O2

## 10. 4.117e-12%
    OH+HO2 -> H2O+O2
    OH+O3 -> HO2+O2
    2(HCNOH+M -> HCN+OH+M)
    2(HCN -> HCN_transport)
    2(HCNOH_transport -> HCNOH)
    Net:O3 + 2HCNOH_transport -> H2O + 2O2 + 2HCN_transport

