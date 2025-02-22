import pandas as pd
import csv
from pprint import pprint
from scipy.stats import f, studentized_range
from tabulate import tabulate
import numpy as np
from matplotlib import pyplot

#def fisher(oxidoNitroso, humedad, temperatura, presion, t):
def fisher(data, t, tableToPrint):
    print(tableToPrint)

    print("\n\n\n")
    #headers = [x[0] for x in data]
    #print(tabulate(table, headers=headers, tablefmt='fancy_outline', colalign=('left')))
    print("============================================================\n")

    #Sumatorias de la Tabla
    Xi = sum([sum(x[1]) for x in data])
    Xi2 = sum([sum([y**2 for y in x[1]]) for x in data])
    Xi2_2 = sum([sum(x[1])**2 for x in data])
    Nt = sum([len(x[1]) for x in data])
    Xi2N = sum([ (sum(x[1]) ** 2) / len(x[1]) for x in data])
    sumMean = sum([ sum(x[1]) / len(x[1]) for x in data])
    
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
            ["∑Xi", round(Xi, 4)],
            ["∑Xi²", round(Xi2, 4)],
            ["∑(Xi)²", round(Xi2_2, 4)],
            ["∑nt", Nt],
            ["∑(Xi)²/n", round(Xi2N, 4)],
            ["∑x̅", round(sumMean, 4)]
        ]
    
    print(tabulate(dataSum, headers=["Sumatoria", 'valor'], tablefmt='grid', colalign=("left", "right")))
    print("============================================================\n")

    #To Print Calc Previos
    print("\--Calculos Previos de la Tabla--/")
    dataPrev = [
            ["C", round(C, 4)],
            ["SCT", round(SCT, 4)],
            ["SCTR", round(SCTR, 4)],
            ["SCE", round(SCE, 4)],
        ]

    print(tabulate(dataPrev, headers=["Calculo", 'valor'], tablefmt='grid', colalign=("left", "right")))
    print("============================================================\n")

    #To Print Tabla para FCalc
    print("\--Calculos para FCalc--/")
    dataFCalc = [
            ["Tratamiento", f"SCTR: {round(SCTR, 4)}", f"t-1: {Gl_Tratamiento}", f"MCTR: {round(MCTR, 4)}", "MCTR / MCE"],
            ["Error", f"SCE: {round(SCE, 4)}", f"n-t: {Gl_Error}", f"MCE: {round(MCE, 2)}", f"FCalc: {round(FCalc, 4)}"],
            ["Total", f"SCT: {round(SCT, 4)}", f"n-1: {Gl_Tratamiento + Gl_Error}"],
        ]

    print(tabulate(dataFCalc, headers=["Fuente de Variacion", 'SC', "gl", "F(rv)", "FCalc"], tablefmt='grid', colalign=("left")))
    print("============================================================\n")

    print(f"FTab = {round(FTab, 4)} y FCalc = {round(FCalc, 4)}")
    if(FCalc < FTab):
        print("Rechazamos la Hipotesis Ha y Aceptamos la Hipotesis Ho")
    else:
        print("Rechazamos la Hipotesis Ho y Aceptamos la Hipotesis Ha")
        
    print("============================================================\n")

    ans = input("Desea ver el Grafico de la Distribucion de Fisher? Si(1) ")
    if ans == "1":
        x = np.linspace(0, FTab+10, 1000)
        y = f.pdf(x, t, Gl_Error)
        pyplot.figure(figsize=(12, 6))
        pyplot.plot(x, y, label="Distribucion f de fisher")
        pyplot.axvline(FTab, color="goldenrod", label=f"FTab: {round(FTab, 4)}", linestyle='dashed')
        pyplot.axvline(FCalc, color="indigo", label=f"FCalc: {round(FCalc, 4)}", linestyle='dashed')
        pyplot.fill_between(x, y, where=(x<FTab), color="lightgreen", label="Region de Aceptación", alpha=0.5)
        pyplot.fill_between(x, y, where=(x>FTab), color="lightcoral", label="Region de Rechazo")
        pyplot.grid(True)
        pyplot.title("Distribucion f de fisher", fontsize=16)
        pyplot.legend()
        pyplot.show()

    print("\--Calculos con el DHS--/")

    ni = Nt / t
    q = studentized_range.ppf(q=0.95, k=Gl_Tratamiento + 1, df=Gl_Error)
    DHS = round(q * (MCE / ni) ** 0.5, 4)

    indVars = []

    means = [ sum(x[1]) / len(x[1]) for x in data]
    tableTukey = []
    coupleMeans = {}
    for i, x in enumerate(means):
        dataToPrint = []
        test = i
        [dataToPrint.append(" ///// ") for x in range(i+1)]
        while(test < len(means) - 1):
            dataToPrint.append(round(means[i] - means[test+1], 4))
            if means[i] - means[test+1] > DHS:
                indVars.append([data[i], data[test+1]])
            coupleMeans.update({f"x{i+1}x{test+2}": means[i] - means[test+1]})
            test += 1
        tableTukey.append(dataToPrint)

    headers = []
    [headers.append(f"x̅{i+1}") for i in range(len(data))]
    print(tabulate(tableTukey, headers=headers, tablefmt='grid', colalign=("left")))

    

    print(f"\nDHS: {DHS}\n")

    for var in coupleMeans:
        if coupleMeans[var] > DHS:
            print(f"Debido a que {var} > DHS es linealmente Independiente")
        else:
            print(f"Debido a que {var} < DHS es linealmente Dependiente")
        print("-----------------------------------------------------")
    
    print("============================================================\n")

    
    [linealRegression(x) for x in indVars]

    multipleRegression([data[1], data[2], data[3]])


#===============================================================================
def linealRegression(variables):
    print(f"\--Correlación Lineal entre {variables[0][0]} y {variables[1][0]}--/\n")
    
    x = variables[0][1]
    y = variables[1][1]

    n = len(variables[0][1])
    
    sumX = sum(x)
    sumY = sum(y)
    
    sumX2 = sum([i**2 for i in x])    
    sumY2 = sum([i**2 for i in y])
    
    sumXY = sum([x[i] * y[i] for i, _ in enumerate(x)])
    
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
    print(f"\nel resultado de r es {round(r, 4)} lo que quiere decir que:\n{correlation(r)}\n")

    
    
    #b = ((sumX * sumY) - (n * sumX * sumY)) / ((sumX**2) - n * sumX2)
    b = (sumXY - (sumX * sumY) / n) / (sumX2 - (sumX**2) / n)
    a = (sumY - (sumX*b)) / (n)
    
    yP = lambda x: a + b * x

    regressionData = [
        ["r", round(r, 4)],
        ["a", round(a, 4)],
        ["b", round(b, 4)]
    ]

    print(tabulate(regressionData, tablefmt='grid'))
    print(f"La ecuacion resultante es y = {round(a, 4)} + {round(b, 4)} * x\n")
    
    if  input("Desear Ver el Grafico? (1)Si\n") == "1":
        pyplot.figure(figsize=(10, 6))
        pyplot.scatter(x, y, color='indigo', label='datos', alpha=0.6)
        pyplot.plot(x, [a+b*x for x in x], color="darkred", label="Regresion likneal")
        pyplot.ylabel(variables[1][0])
        pyplot.xlabel(variables[0][0])
        pyplot.legend()
        pyplot.title(f"Regresion lineal para {variables[0][0]} y {variables[1][0]}")
        pyplot.grid(True)
        pyplot.show()

    ans = input("Desea ingresar algun valor x en y? (1)Si ")

    while(ans == '1'):
        numX = float(input("Introduzca el Valor de X "))
        print(f"El Resultado es {round(yP(numX), 4)}")
        ans = input("Desea ingresar algun valor x en y? (1)Si")

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

    #sumX1X2 = sum([variables[1][1][i] ** variables[2][1][i] for i in range(len(variables[0][1]))])
    sumX1X2 = sum([x * y for x, y in zip(variables[1][1], variables[2][1])])

    #sumX1Y = sum([variables[0][1][i] ** variables[1][1][i] for i, _ in enumerate(variables[0][1])])
    sumX1Y = sum([x * y for x, y in zip(variables[1][1], variables[0][1])])

    sumX2Sq = sum([x**2 for x in variables[2][1]])

    #sumX2Y = sum([variables[0][1][i] ** variables[2][1][i] for i, _ in enumerate(variables[0][1])])
    sumX2Y = sum([x * y for x, y in zip(variables[2][1], variables[0][1])])

    #varsToMatrix =[n, sumX1, sumX2, sumY]
    matrix = [
        [n, sumX1, sumX2, sumY],
        [sumX1, sumX1Sq, sumX1X2, sumX1Y],
        [sumX2, sumX1X2, sumX2Sq, sumX2Y]
    ]
    print("Matriz Ampliada")
    print(tabulate(matrix))

    gaussJordan(matrix)

            
def gaussJordan(matrix):
    for i in range(len(matrix)):
        # if matrix[i][i] != 1:
        matrix[i] = [(1/matrix[i][i]) * x for x in matrix[i]]
        
        for j in range(len(matrix)):
            if j == i: continue
            negValue = -matrix[j][i]
            matrix[j] = [y + negValue * x for x, y in zip(matrix[i], matrix[j])]
    
    print("\nMatriz Resuelta")
    print(tabulate(matrix))

    matrixResults = [
        ["B0", round(matrix[0][-1], 4)],
        ["B1", round(matrix[1][-1], 4)],
        ["B2", round(matrix[2][-1], 4)]
    ]

    print("Resultados")
    print(tabulate(matrixResults))
    print(f"Ecuacion Resultado y =  {matrixResults[0][1]} + ({matrixResults[1][1]} * x3) + ({matrixResults[2][1]} * x4)")
    yRM = lambda x3,x4: matrixResults[0][1] + (matrixResults[1][1] * x3) + (matrixResults[2][1] * x4)
    ans = input("Desea ingresar Valores a la Ecuacion? Si(1) ")

    while(ans == "1"):
        x3 = float(input("Introduzca el valor de x3 "))
        x4 = float(input("Introduzca el valor de x4 "))
        print(f"El Resultado es {round(yRM(x3,x4), 4)}")
        ans = input("Desea ingresar Valores a la Ecuacion? Si(1) ")

if __name__ == '__main__':
    datos = pd.read_csv("datos.csv", delimiter=";")
    tableData = []
    for x in datos.columns:
        aux = [x]
        aux.append(datos[x].tolist())
        tableData.append(aux)

    fisher(tableData, len(tableData), datos)
