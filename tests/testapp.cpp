#include <iostream>

int main()
{
const char* input_str = R"input(
{
  "input": {
    "project": {
            "pname": "CO_Test2",
            "rootDir": "/home/jimmy/Programs/HalIR/TestCalc/",
            "hapi_db": "",
            "pcomments": "CO Test on some different concentrations",
            "pfiles": ["F0", "F1"]
    },
    "sampleEnv": {
            "temp": 303.0,
            "tempU": "K",
            "press": 1.0,
            "pressU": "atm",
            "pathL": 300.0,
            "pathLU": "m",
            "ROI": [2000.0, 2245.0],
            "res": 0.1,
            "apod": "Boxcar",
            "fov": 0.1,
            "ftype": "Transmission",
            "bgfile": ""
    },
    "composition": [
            {"molec": "CO",
             "isotop": "Natural",
             "vmr": 5.921539600296e-05,
             "amountU": "mb",
             "prmfile": ""},
            {"molec": "H2O",
             "isotop": "Natural",
             "vmr": 5.921539600296e-05,
             "amountU": "mb",
             "prmfile": ""}
    ]
  }
})input";
  std::cout << input_str << std::endl;
  return 0;
}