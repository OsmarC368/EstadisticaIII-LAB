import csv
from pprint import pprint
from scipy.stats import f, studentized_range
from tabulate import tabulate
import numpy as np


def fisher(oxidoNitroso, humedad, temperatura, presion, t):
    #To Print Base Table
    dataTable = []
    for i in range(len(oxidoNitroso)):
        dataTable.append([round(oxidoNitroso[i], 2),
                        round(humedad[i], 2), 
                        round(temperatura[i], 2), 
                        round(presion[i], 2)
                        ])
    print("\n\n\n")
    print(tabulate(dataTable, headers=["Oxido Nitroso", "Humedad", "Temperatura", "Presion"], tablefmt='fancy_outline', colalign=('left')))
    print("============================================================\n")

    #Sumatorias de la Tabla
    Xi = sum(oxidoNitroso) + sum(humedad) + sum(temperatura) + sum(presion)
    Xi2 = sum([x**2 for x in oxidoNitroso]) + sum([x**2 for x in humedad]) +sum([x**2 for x in temperatura]) + sum([x**2 for x in presion])
    Xi2_2 = (sum(oxidoNitroso)**2) + (sum(humedad)**2) + (sum(temperatura)**2) + (sum(presion)**2)
    Nt = len(oxidoNitroso) + len(humedad) + len(temperatura) + len(presion)
    Xi2N = (sum(oxidoNitroso)**2 / len(oxidoNitroso)) + (sum(humedad)**2 / len(humedad)) + (sum(temperatura)**2 / len(temperatura)) + (sum(presion)**2 / len(presion))
    sumMean = sum(oxidoNitroso) / len(oxidoNitroso) + sum(humedad) / len(humedad) + sum(temperatura) / len(temperatura) + sum(presion) / len(presion)
    
    #Variables previas a la tabla de calc
    C = Xi**2 / Nt
    SCT = Xi2 - C
    SCTR = Xi2N - C
    SCE = Xi2 -Xi2N

    #FTab calcs
    Gl_Tratamiento = t - 1
    Gl_Error = Nt - t
    FTab = f.ppf(1-0.05, dfn=Gl_Tratamiento, dfd=Gl_Error)

    #Calcs para FCalc
    MCTR = SCTR / (t - 1)
    MCE = SCE / (Nt - t)
    FCalc = MCTR / MCE

    #∑
    #To Print Sumatoria Table
    print("\--Sumatorias de la Tabla--/")

    dataSum = [
            ["∑Xi", round(Xi, 2)],
            ["∑Xi²", round(Xi2, 2)],
            ["∑(Xi)²", format(Xi2_2, ".2f")],
            ["∑nt", Nt],
            ["∑(Xi)²/n", round(Xi2N, 2)],
            ["∑x̅", round(sumMean, 2)]
        ]
    
    print(tabulate(dataSum, headers=["Sumatoria", 'valor'], tablefmt='grid', colalign=("left", "right")))
    print("============================================================\n")

    #To Print Calc Previos
    print("\--Calculos Previos de la Tabla--/")
    dataPrev = [
            ["C", round(C, 2)],
            ["SCT", round(SCT, 2)],
            ["SCTR", round(SCTR, 2)],
            ["SCE", round(SCE, 2)],
        ]

    print(tabulate(dataPrev, headers=["Calculo", 'valor'], tablefmt='grid', colalign=("left", "right")))
    print("============================================================\n")

    #To Print Tabla para FCalc
    print("\--Calculos para FCalc--/")
    dataFCalc = [
            ["Tratamiento", f"SCTR: {round(SCTR, 2)}", f"t-1: {Gl_Tratamiento}", f"MCTR: {round(MCTR, 2)}", "MCTR / MCE"],
            ["Error", f"SCE: {round(SCE, 2)}", f"n-t: {Gl_Error}", f"MCE: {round(MCE, 2)}", f"FCalc: {round(FCalc, 2)}"],
            ["Total", f"SCT: {round(SCT, 2)}", f"n-1: {Gl_Tratamiento + Gl_Error}"],
        ]

    print(tabulate(dataFCalc, headers=["Fuente de Variacion", 'SC', "gl", "F(rv)"], tablefmt='grid', colalign=("left")))
    print("============================================================\n")

    if(FCalc < FTab):
        print("Rechazamos la Hipotesis Ha y Aceptamos la Hipotesis Ho")
    else:
        print("Rechazamos la Hipotesis Ho y Aceptamos la Hipotesis Ha")
        
    print("============================================================\n")
    print("\--Calculos con el DHS--/")
    #DHS = f.
    #print(f"DHS: {DHS}")
    means = [np.mean(oxidoNitroso), np.mean(humedad), np.mean(temperatura), np.mean(presion)]
    #dataExample = " ///// "
    tableTukey = []

    for i, x in enumerate(means):
        data = []
        test = i
        [data.append(" ///// ") for x in range(i+1)]
        while(test < len(means) - 1):
            data.append(round(means[i] - means[test+1], 2))
            test += 1
        tableTukey.append(data)
        

    # tableTukey = [
    #     [" ///// ", round(means[0] - means[1], 2), round(means[0] - means[2], 2), round(means[0] - means[3], 2)],
    #     [" ///// ", " ///// ", round(means[1] - means[2], 2), round(means[1] - means[3], 2)],
    #     [" ///// ", " ///// ", " ///// ", round(means[2] - means[3], 2)],
    #     [" ///// ", " ///// ", " ///// ", " ///// "]
    # ]
        
    print(tabulate(tableTukey, headers=["x̅1", 'x̅2', "x̅3", "x̅4"], tablefmt='grid', colalign=("left")))

    # with open("test.txt", "w") as fileToCreate:
    #     fileToCreate.write(f"""Xi: {Xi}
    #                         Xi^2: {Xi2}
    #                         (Xi)2: {Xi2_2}
    #                         n: {Nt}
    #                         xi2/n: {Xi2N}
    #                         Mean X: {sumMean}
    #                         """)
    #print(f"Ftab: {fTab}")


if __name__ == "__main__":

    oxidoNitroso = []
    humedad = []
    temperatura = []
    presion = []
    with open("datos.csv", newline="") as dataBase:
         [[oxidoNitroso.append(float(x["Oxido_nitroso"])), 
           humedad.append(float(x["Humedad(x1)"])), 
           temperatura.append(float(x["Temperatura(x2)"])), 
           presion.append(float(x["Presion(x3)"]))] for x in list(csv.DictReader(dataBase, delimiter=";"))]
         
    fisher(oxidoNitroso, humedad, temperatura, presion, 4)
