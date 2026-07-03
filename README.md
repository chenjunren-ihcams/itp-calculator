# itp-calculator

This project is the implementation of daGOAT model-defined trinary classification presented in "**Dynamic forecasting of severe acute graft-versus-host disease after transplantation**".

## 1. SYSTEM REQUIREMENTS

Operating system that has been tested on： Microsoft Windows 10

Software version that has been tested on： R version 4.4.1 (2024-06-14)

Required R libraries: readxl ==1.4.3 e1071 == 1.7.16


## 2. INSTALLATION GUIDE

Not applicable.


## 3. INSTRUCTIONS

Mock patient data consist of `ITP_encoded.csv` (subjects' sex, age, confirmed response status, and time to confirmed response) and `patients.zip` (daily values of 20 hematological parameters of subjects).

Extract `./mock-up-data/patients.zip` to the `./mock-up-data/patients/`.

Set `root` on line 409 of `itp-calculator.R` to the absolute path of the `mock_up_data` folder under your project directory.

Run `itp-calculator.R`.

daGOAT model-defined trinary classification will be stored in the data.frame `df_risk` and can be accessed using `df_risk$risk_bind`.
