import csv
from pprint import pprint
from scipy.stats import f, studentized_range
from tabulate import tabulate
import numpy as np
from matplotlib import pyplot


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
    
    means = [np.mean(oxidoNitroso), np.mean(humedad), np.mean(temperatura), np.mean(presion)]
    tableTukey = []
    coupleMeans = {}
    for i, x in enumerate(means):
        data = []
        test = i
        [data.append(" ///// ") for x in range(i+1)]
        while(test < len(means) - 1):
            data.append(round(means[i] - means[test+1], 2))
            coupleMeans.update({f"x{i+1}x{test+2}": means[i] - means[test+1]})
            test += 1
        tableTukey.append(data)

    print(tabulate(tableTukey, headers=["x̅1", 'x̅2', "x̅3", "x̅4"], tablefmt='grid', colalign=("left")))

    ni = Nt / t
    q = studentized_range.ppf(q=0.95, k=Gl_Tratamiento + 1, df=Gl_Error)
    DHS = round(q * (MCE / ni) ** 0.5, 4)

    print(f"\nDHS: {DHS}\n")

    for var in coupleMeans:
        if coupleMeans[var] > DHS:
            print(f"Debido a que {var} > DHS es linealmente Independiente")
        else:
            print(f"Debido a que {var} < DHS es linealmente Dependiente")
        print("-----------------------------------------------------")
    
    print("============================================================\n")

    x2x4 = [["Humedad", humedad], ["Presion", presion]]
    x3x4 = [["Temperatura", temperatura], ["Presion", presion]]
    
    dataToDo = [x2x4, x3x4]
    
    [linealRegression(x) for x in dataToDo]

    multipleRegression([["Humedad", humedad], ["Temperatura", temperatura], ["Presion", presion]])


#===============================================================================
def linealRegression(variables):
    print(f"\--Correlación Lineal entre {variables[0][0]} y {variables[1][0]}--/\n")
    n = len(variables[0][1])
    
    sumX = sum(variables[0][1])
    sumY = sum(variables[1][1])
    
    sumX2 = sum([x**2 for x in variables[0][1]])    
    sumY2 = sum([x**2 for x in variables[1][1]])
    
    sumXY = sum([variables[0][1][i] * variables[1][1][i] for i, _ in enumerate(variables[0][1])])
    
    tabulateData = [
        ["n",round(n, 4)],
        ["∑X", round(sumX, 4)],
        ["∑X²", round(sumX2, 4)],
        ["∑Y", round(sumY, 4)],
        ["∑Y²", round(sumY2, 4)],
        ["∑X.Y", round(sumXY, 4)]
    ]

    print("\-Datos para la Correlación-/")
    print(tabulate(tabulateData, tablefmt='grid'))

    r = (sumXY - (sumX * sumY / n)) / np.sqrt((sumX2 - ((sumX**2) / n)) * (sumY2 - ((sumY**2) / n)))
    print(f"\nel resultado de r es {r} lo que quiere decir que:\n{correlation(r)}\n")
    
    b = ((sumX * sumY) - (n * sumX * sumY)) / ((sumX**2) - n * sumX2)
    a = (sumY - (sumX*b)) / (n)
    
    y = lambda x: a + b * x

    #randomX = np.random.normal(50.0, 1.0, 50)
    randomX = np.random.randint(1000000000, size=(20))
    #print(randomX)
    
    if  input("Desear Ver el Grafico? (1)Si\n") == "1":
        xRange = range(-500, 2000)
        #pyplot.plot(xRange, [y(i) for i in xRange])
        #pyplot.axhline(0, color="red")
        #pyplot.axvline(0, color="green")
        #pyplot.scatter(randomX, [y(i) for i in randomX])
        pyplot.scatter(variables[0][1], variables[1][1])
        #pyplot.xlim(-10, 10)
        #pyplot.ylim(-10, 10)
        pyplot.show()
    print("============================================================\n")


def correlation(r):
    if r == 0:
        return "Sin Relación"
    elif 0 < r < 0.25:
        return "Debíl Correlación Directa"
    elif 0.25 <= r < 0.75:
        return "Intermedia Correlación Directa"
    elif 0.75 <= r < 1:
        return "Fuerte Correlación Directa"
    elif r == 1:
        return "Perfecta Correlación Directa"
    elif -0.25 < r < 0:
        return "Debíl Correlación Indirecta"
    elif -0.75 < r <= -0.25:
        return "Intermedia Correlación Indirecta"
    elif -1 < r <= -0.75:
        return "Fuerte Correlación Indirecta"
    elif r == -1:
        return "Correlacion Perfecta Indirecta"
    
def multipleRegression(variables):
    print(f"\--Regresion Multiple entre {variables[0][0]}, {variables[1][0]} y {variables[2][0]}--/\n")

    n = len(variables[0][1])

    sumY = sum(variables[0][1])
    sumX1 = sum(variables[1][1])
    sumX2 = sum(variables[2][1])

    sumX1Sq = sum([x**2 for x in variables[1][1]])

    sumX1X2 = sum([variables[1][1][i] ** variables[2][1][i] for i, _ in enumerate(variables[0][1])])

    sumX1Y = sum([variables[0][1][i] ** variables[1][1][i] for i, _ in enumerate(variables[0][1])])

    sumX2Sq = sum([x**2 for x in variables[2][1]])

    sumX2Y = sum([variables[0][1][i] ** variables[2][1][i] for i, _ in enumerate(variables[0][1])])

    #varsToMatrix =[n, sumX1, sumX2, sumY]
    matrix = [
        [n, sumX1, sumX2, sumY],
        [sumX1, sumX1Sq, sumX1X2, sumX1Y],
        [sumX2, sumX1X2, sumX2Sq, sumX2Y]
    ]

    print(tabulate(matrix))

    gaussJordan(matrix)
    # for i in range(len(matrix)):
    #     if matrix[i][i] != 1:
    #         matrix[i] = [(1/matrix[i][i]) * x for x in matrix[i]]

    #for i in range(len(variables)):
    #    for x in range(i+1):
    #        matrix.append([varsToMatrix[x]])
            
def gaussJordan(matrix):
    for i in range(len(matrix)):
        # if matrix[i][i] != 1:
        matrix[i] = [(1/matrix[i][i]) * x for x in matrix[i]]
        
        for j in range(len(matrix)):
            if j == i: continue
            negValue = -matrix[j][i]
            matrix[j] = [y + (x * negValue) for x, y in zip(matrix[j], matrix[i])]
    
    
    pprint(matrix)
            


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
