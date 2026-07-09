# itp-calculator

This project is the implementation of daGOAT model-defined trinary classification presented in "**Artificial intelligence prediction of therapy response in newly diagnosed adult severe primary immune thrombocytopenia**".

## 1. SYSTEM REQUIREMENTS

Operating system that has been tested on： Microsoft Windows 10

Software version that has been tested on： R version 4.4.1 (2024-06-14)

Required R libraries: readxl ==1.4.3 e1071 == 1.7.16


## 2. INSTALLATION GUIDE

Not applicable.


## 3. INSTRUCTIONS

Desensitized data consist of `ITP_encoded.csv` (subjects' sex and age) and `patients.zip` (daily values of 20 hematological parameters of subjects).

Extract `./desensitized-data/patients.zip` to the `./desensitized-data/patients/`.

Set `root` on line 409 of `itp-calculator.R` to the absolute path of the `desensitized-data` folder under your project directory.

Run `itp-calculator.R`.

daGOAT model-defined trinary classification will be stored in the data.frame `df_risk` and can be accessed using `df_risk$risk_bind`.
