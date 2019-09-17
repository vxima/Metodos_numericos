#Bibliotecas importadas

import numpy
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import sys
import matplotlib.pyplot as plot

#######################METODOS#################################



print("Metodo de Euler")
#Yn+1 = Yn + hf(tn, yn)
def Euler( y0 , t0 , h , n , f):#Ta funcionando
    
    t,y = symbols('t y')
    tlist=[t0]
    ylist=[y0]
    for j in range(n):
        
        fn = f.subs({t:t0, y:y0})
        y0 = y0 + h*fn
        t0 = t0 + h
        ylist.append(float(y0)) 
        tlist.append(float(t0))
    return y0 


'''funcao= input("Bote sua função ai tio\n")
f = parse_expr(funcao)
print(Euler( 1.0 , 0.0 , 0.05 , 4 , f))'''
print("Metodo de Euler_Inverso")
#Yn+1 = Yn + hf(tn+1 , Yn+1)
def Euler_Inverso ( y0 , t0 , h , n , f):#ta aproximado 
    t,y = symbols('t y')
    tlist=[t0]
    ylist=[y0]
    for j in range(n):
        auxt = t0 + h
        auxy = Euler(y0 , t0 , h , 1 , f)
        fn = f.subs({t:auxt, y:auxy})
        y0 = y0 + h*fn
        t0 = t0 + h
        ylist.append(float(y0)) 
        tlist.append(float(t0))
    return y0 

'''funcao= input("Bote sua função ai tio\n")
f = parse_expr(funcao)
print(Euler_Inverso( 1.0 , 0.0 , 0.05 , 4 , f))'''

print("Metodo de Euler_Aprimorado")
#Yn+1 = Yn + [f(tn , Yn) + f(tn+1 , Yn+1 )]h/2

def Euler_Aprimorado(y0 , t0 , h , n , f): #Ta funcionando
    t,y = symbols('t y')
    tlist=[t0]
    ylist=[y0]
    for j in range(n):
        k1 = f.subs({t:t0, y:y0})

        auxt = t0 + h
        auxy = y0 + h*k1

        k2 = f.subs({t:auxt, y:auxy})
        y0 = y0 + (h/2)*(k1+k2)
        t0 = t0 + h
        ylist.append(float(y0)) 
        tlist.append(float(t0))
    return y0 

'''funcao= input("Bote sua função ai tio\n")
f = parse_expr(funcao)
print(Euler_Aprimorado( 1.0 , 0.0 , 0.025 , 8 , f))'''

print("Runge_Kutta")
#Yn+1 = Yn + (h/6)*(Kn1 + 2Kn2 + 2Kn3 + Kn4)
#Kn1 = f(tn , Yn)
#Kn2 = f(tn + h/2 , Yn + (h/2)*Kn1) 
#Kn3 = f(tn + h/2 , Yn + (h/2)*Kn2)
#Kn4 = f(tn + h , Yn + h*Kn3) 

def Runge_Kutta(y0 , t0 , h , n , f): #Ta funcionando
    t,y = symbols('t y')
    tlist=[t0]
    ylist=[y0]
    for j in range(n):
        #K1
        k1 = f.subs({t:t0, y:y0})

        #K2
        k2_auxt = t0 + h/2
        k2_auxy = y0 + (h/2)*k1
        k2 = f.subs({t:k2_auxt, y:k2_auxy})

        #K3
        k3_auxt = t0 + h/2
        k3_auxy = y0 + (h/2)*k2
        k3 = f.subs({t:k3_auxt, y:k3_auxy})

        #K4
        k4_auxt = t0 + h
        k4_auxy = y0 + h*k3
        k4 = f.subs({t:k4_auxt, y:k4_auxy})


        y0 = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        t0 = t0 + h
        ylist.append(float(y0)) 
        tlist.append(float(t0))
    return y0 

'''funcao= input("Bote sua função ai tio\n")
f = parse_expr(funcao)
print(Runge_Kutta( 1.0 , 0.0 , 0.05 , 4 , f))'''

print("Adams_Bashforth")

def Adams_Bashforth(y0 , t0 , h , n , f , ordem): #ta imcompleto

    coefM = [
			[1.0],
			[3.0/2.0,-1.0/2.0],
			[23.0/12.0,-4.0/3.0,5.0/12.0],
			[55.0/24.0,-59.0/24.0,37.0/24.0,-3.0/8.0],
			[1901.0/720.0,-1387.0/360.0,109.0/30.0,-637.0/360.0,251.0/720.0],
			[4277.0/1440.0,-2641.0/480.0,4991.0/720.0,-3649.0/720.0,959.0/480.0,-95.0/288.0],
			[198721.0/60480.0,-18637.0/2520.0,235183.0/20160.0,-10754.0/945.0,135713.0/20160.0,-5603.0/2520.0,19087.0/60480.0],
			[16083.0/4480.0,-1152169.0/120960.0,242653.0/13440.0,-296053.0/13440.0,2102243.0/120960.0,-115747.0/13440.0,32863.0/13440.0,-5257.0/17280.0]
			]

    t,y = symbols('t y')
