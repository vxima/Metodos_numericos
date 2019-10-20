#Bibliotecas importadas

import numpy
from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import sys
import matplotlib.pyplot as matplot


####################### FUNCOES IMPLICITAS ############################


def EulerImplicito( y0 , t0 , h  , f):
    #usado para euler inverso em que se calcula o Yn+1
    t,y = symbols('t y')
    fn = f.subs({t:t0, y:y0})
    yn = y0 + h*fn
    
    return yn #retorna yn+1


def Adam_BashforthImplicito(Yant , t0 , h  , f , ordem): 
    #usado para Adams_Multon em que se calcula Yn+1
    

    coeficientes = [
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

    tn = t0
    yn = float(Yant[len(Yant) - 1])

    
    
            
    

    
    inclinacao = 0.0

    for k in range(len(coeficientes[ordem - 1])):
        auxiliar = h*coeficientes[ordem -1][k]*f.subs({t:(tn - (h*k)),y:float(Yant[len(Yant)-k-1])})
        inclinacao = inclinacao + auxiliar

    yn = yn + inclinacao
        
    return yn #retorna yn+1
######################FUNCOES QUE RETORNAM LISTAS #####################

def ListaEuler(y0, t0 , h , f , ordem):
    Ylist = []
    Ylist.append(float(y0))

    t,y = symbols('t y')
    n = ordem -1
    for i in range(n):
        fn = f.subs({t:t0, y:y0})
        y0 = y0 + h*fn
        t0 = t0 + h
        Ylist.append(y0)

    return Ylist

def ListaEuler_Inverso( y0 , t0 , h , f , ordem):
    t,y = symbols('t y')
    
    ylist=[]
    ylist.append(float(y0))
    n = ordem -1
    for j in range(n):
        auxt = t0 + h
        auxy = EulerImplicito(y0 , t0 , h , f)
        fn = f.subs({t:auxt, y:auxy})
        y0 = y0 + h*fn
        t0 = t0 + h
        ylist.append(y0) 
        
    return ylist   

def ListaEuler_aprimorado(y0 , t0 , h , f , ordem): 
    t,y = symbols('t y')
    ylist=[]
    ylist.append(float(y0)) 
    n = ordem -1
    for j in range(n):
        #K2
        k1 = f.subs({t:t0, y:y0})

        #K2
        auxt = t0 + h
        auxy = y0 + h*k1
        k2 = f.subs({t:auxt, y:auxy})


        y0 = y0 + (h/2)*(k1+k2)
        t0 = t0 + h
        ylist.append(y0) 

    return ylist

def ListaRunge_Kutta(y0 , t0 , h , f , ordem):
    t,y = symbols('t y')
    ylist=[]
    ylist.append(float(y0)) 
    n = ordem -1
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
        ylist.append(y0) 

    return ylist


#######################METODOS##############################################



#Yn+1 = Yn + hf(tn, yn)
def Euler( y0 , t0 , h , n , f):#Ta funcionando
    print('Metodo de Euler')
    print('y({}) = {}'.format(t0,y0))
    print('h = {}'.format(h))

    t,y = symbols('t y')
    tlist=[t0]
    ylist=[y0]
    for j in range(n):
        print(j,' ',y0)
        fn = f.subs({t:t0, y:y0})
        y0 = y0 + h*fn
        t0 = t0 + h
        ylist.append(float(y0)) 
        tlist.append(float(t0))
    
    matplot.xlabel('t')
    matplot.ylabel('y(t)')
    matplot.title('Metodo de Euler')
    matplot.plot(tlist , ylist , color='blue')
    matplot.show()
    return  


'''funcao= input("Bote sua função\n")
f = parse_expr(funcao)
print(Euler( 1.0 , 0.0 , 0.05 , 4 , f))'''
#Yn+1 = Yn + hf(tn+1 , Yn+1)
def Euler_Inverso ( y0 , t0 , h , n , f):#ta aproximado 
    print("Metodo de Euler_Inverso")
    print('y({}) = {}'.format(t0,y0))
    print('h = {}'.format(h))

    t,y = symbols('t y')
    tlist=[t0]
    ylist=[y0]
    for j in range(n+1):
        print(j,' ',y0)
        auxt = t0 + h
        auxy = EulerImplicito(y0 , t0 , h  , f)
        fn = f.subs({t:auxt, y:auxy})
        y0 = y0 + h*fn
        t0 = t0 + h
        ylist.append(y0) 
        tlist.append(t0)
    
    matplot.xlabel('t')
    matplot.ylabel('y(t)')
    matplot.title('Metodo de Euler Inverso')
    matplot.plot(tlist , ylist , color='red')
    matplot.show()
    return 

'''funcao= input("Bote sua função ai\n")
f = parse_expr(funcao)
print(Euler_Inverso( 1.0 , 0.0 , 0.05 , 4 , f))'''

#Yn+1 = Yn + [f(tn , Yn) + f(tn+1 , Yn+1 )]h/2

def Euler_Aprimorado(y0 , t0 , h , n , f): #Ta funcionando
    print("Metodo de Euler_Aprimorado")
    print('y({}) = {}'.format(t0,y0))
    print('h = {}'.format(h))

    t,y = symbols('t y')
    tlist=[t0]
    ylist=[y0]
    for j in range(n):
        print(j,' ',y0)
        k1 = f.subs({t:t0, y:y0})

        auxt = t0 + h
        auxy = y0 + h*k1

        k2 = f.subs({t:auxt, y:auxy})
        y0 = y0 + (h/2)*(k1+k2)
        t0 = t0 + h
        ylist.append(y0) 
        tlist.append(t0)

    matplot.xlabel('t')
    matplot.ylabel('y(t)')
    matplot.title('Metodo de Euler Aprimorado')
    matplot.plot(tlist , ylist , color='yellow')
    matplot.show()
    return 

'''funcao= input("Bote sua função ai \n")
f = parse_expr(funcao)
print(Euler_Aprimorado( 1.0 , 0.0 , 0.025 , 8 , f))'''

#Yn+1 = Yn + (h/6)*(Kn1 + 2Kn2 + 2Kn3 + Kn4)
#Kn1 = f(tn , Yn)
#Kn2 = f(tn + h/2 , Yn + (h/2)*Kn1) 
#Kn3 = f(tn + h/2 , Yn + (h/2)*Kn2)
#Kn4 = f(tn + h , Yn + h*Kn3) 

def Runge_Kutta(y0 , t0 , h , n , f): #Ta funcionando
    print("Runge_Kutta")
    print('y({}) = {}'.format(t0,y0))
    print('h = {}'.format(h))

    t,y = symbols('t y')
    tlist=[t0]
    ylist=[y0]
    for j in range(n):
        print(j,' ',y0)
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

    matplot.xlabel('t')
    matplot.ylabel('y(t)')
    matplot.title('Metodo de Runge-Kutta')
    matplot.plot(tlist , ylist , color='green')
    matplot.show()
    return 

'''#funcao que testa 
funcao= input("Bote sua função ai \n")
f = parse_expr(funcao)
print(Runge_Kutta( 1.0 , 0.0 , 0.05 , 4 , f))'''



def Adam_Bashforth(Yant , t0 , h , n , f , ordem , complemento=''): #ta funcionando
    #ta pegando para ordem = numero de pontos anteriores

    print("Metodo de Adam-Bashforth"+complemento)
    print('y({}) = {}'.format(t0,Yant[0]))
    print('h = {}'.format(h))

    coeficientes = [
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

    tn = t0
    yn = float(Yant[len(Yant) - 1]) #pega ultimo valor da lista de y anteriores

    ylist = []
    tlist = []
    
            
    for j in range(ordem):
        print(j,' ',Yant[j])
        ylist.append(Yant[j])
        tlist.append(tn)
        if(j>0):
            tn = tn + h

    for j in range(ordem , n+1 , 1): #for(i=ordem , i<=n+1 , i++)
        inclinacao = 0.0

        for k in range(len(coeficientes[ordem - 1])):
            auxiliar = h*coeficientes[ordem -1][k]*f.subs({t:(tn - (h*k)),y:float(Yant[len(Yant)-k-1])})
            inclinacao = inclinacao + auxiliar

        yn = yn + inclinacao
        tn = tn + h 
        Yant.append(yn)

        print(j,' ', yn)
        ylist.append(yn)
        tlist.append(tn)

    matplot.xlabel('t')
    matplot.ylabel('y(t)')
    matplot.title('Metodo de Adam-Bashforth' + complemento)
    matplot.plot(tlist , ylist , color='orange')
    matplot.show()
    
    return

def Adam_Multon(Yant, t0 , h , n , f , ordem , complemento=''): #ta funcionando
    #ta pegando para ordem = numero de pontos anteriores

    print('Metodo de Adam-multon'+complemento)
    print('y({}) = {}'.format(t0,Yant[0]))
    print('h = {}'.format(h))

    t,y = symbols('t y')

    coeficientes = [
			[1.0],
			[1.0/2.0,1.0/2.0],
			[5.0/12.0,2.0/3.0,-1.0/12.0],
			[3.0/8.0,19.0/24.0,-5.0/24.0,1.0/24.0],
			[251.0/720.0,323.0/360.0,-11.0/30.0,53.0/360.0,-19.0/720.0],
			[95.0/288.0,1427.0/1440.0,-133.0/240.0,241.0/720.0,-173.0/1440.0,3.0/160.0],
			[19087.0/60480.0,2713.0/2520.0,-15487.0/20160.0,586.0/945.0,-6737.0/20160.0,263.0/2520.0,-863.0/60480.0],
			[5257.0/17280.0,139849.0/120960.0,-4511.0/4480.0,123133.0/120960.0,-88547.0/120960.0,1537.0/4480.0,-11351.0/120960.0,275.0/24192.0]
			]
    

    tn = t0
    yn = float(Yant[len(Yant) - 1])

    ylist = []
    tlist = []

    for j in range(ordem):
        print(j,' ',Yant[j])
        ylist.append(Yant[j])
        tlist.append(tn)
        if(j>0):
            tn = tn + h

    for j in range(ordem , n+1 , 1): #for(i=ordem , i<=n+1 , i++)
        inclinacao = 0.0

        #parte que calcula Yn+1
        t1 = tn + h
        previsao = Adam_BashforthImplicito(Yant , tn , h , f , ordem)
        Yn1  = h*coeficientes[ordem-1][0]*f.subs({t:t1 , y:previsao})

        inclinacao = inclinacao + Yn1

        for k in range(1 , len(coeficientes[ordem - 1]), 1):
            auxiliar = h*coeficientes[ordem -1][k]*f.subs({t:(tn - (h*(k-1))),y:float(Yant[len(Yant)-(k-1)-1])})
            inclinacao = inclinacao + auxiliar


        yn = yn + inclinacao
        tn = tn + h 

        Yant.append(yn)

        print(j,' ', yn)
        ylist.append(yn)
        tlist.append(tn)
    
    matplot.xlabel('t')
    matplot.ylabel('y(t)')
    matplot.title('Metodo de Adam-Multon' + complemento)
    matplot.plot(tlist , ylist , color='purple')
    matplot.show()
    return

def Formula_inversa(Yant , t0 , h , n , f , ordem , complemento=''): #ta funcionando
    #ta pegando para ordem = numero de pontos anteriores

    print('Metodo de Formula Inversa de Diferenciacao'+complemento)
    print('y({}) = {}'.format(t0,Yant[0]))
    print('h = {}'.format(h))

    yn = y0

    coeficientes_Yn = [
            [1.0],
			[4.0/3.0,-1.0/3.0],
			[18.0/11.0,-9.0/11.0,2.0/11.0],
			[48.0/25.0,-36.0/25.0,16.0/25.0,-3.0/25.0],
			[300.0/137.0,-300.0/137.0,200.0/137.0,-75.0/137.0,12.0/137.0],
			[360.0/147.0,-450.0/147.0,400.0/147.0,-225.0/147.0,72.0/147.0,-10.0/147.0],
			]

    coeficientes_Fn = [
            [1.0],
			[2.0/3.0],
			[6.0/11.0],
			[12.0/25.0],
			[60.0/137.0],
			[60.0/147.0],
			]
    
    tn = t0
    
    t,y = symbols('t y')
    ylist = []
    tlist = []
    for j in range(ordem):
        print(j,' ',Yant[j])
        ylist.append(Yant[j])
        tlist.append(tn)
        if(j>0):
            tn = tn + h

    for j in range(ordem , n+1 , 1): #for(i=ordem , i<=n+1 , i++)

        inclinacao = 0.0

        #parte que calcula fn+1
        t1 = tn + h
        previsao = Adam_BashforthImplicito(Yant , tn , h , f , ordem)
        fn1  = h*coeficientes_Fn[ordem-1][0]*f.subs({t:t1 , y:previsao})

        inclinacao = inclinacao + fn1

        
        #parte que calcula os yn
        for k in range(len(coeficientes_Yn[ordem - 1])):
            auxiliar = coeficientes_Yn[ordem -1][k]*Yant[len(Yant)-k-1]
            inclinacao = inclinacao + auxiliar


        yn = inclinacao
        tn = tn + h 
        Yant.append(yn)

        print(j,' ', yn)
        ylist.append(yn)
        tlist.append(tn)
    
    matplot.xlabel('t')
    matplot.ylabel('y(t)')
    matplot.title('Metodo de Formula Inversa de Diferenciacao' + complemento)
    matplot.plot(tlist , ylist , color='violet')
    matplot.show()
    return 

#################### MAIN ###############################


if __name__ == "__main__":

    for arquivo in sys.stdin: # lê cada linha do arquivo
        entrada = str(arquivo).split() #divide a linha da entrada

        metodo = entrada[0] 

        if(metodo == 'euler'):
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])

            Euler(y0 , t0 , h , n , f)
            print()

        elif(metodo == 'euler_inverso'):
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])

            Euler_Inverso(y0 , t0 , h , n , f)
            print()    

        elif(metodo == 'euler_aprimorado'):
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])

            Euler_Aprimorado(y0 , t0 , h , n , f)
            print()  

        elif(metodo == 'runge_kutta'):
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])

            Runge_Kutta(y0 , t0 , h , n , f)
            print()  

        elif(metodo == 'adam_bashforth'):#ta pegando para ordem = numero de pontos anteriores

            ylist = []
            tam = len(entrada)

            
            #funcao que pega os valores anteriores 
            for j in range(int(entrada[tam -1])):
                ylist.append(entrada[j+1])
                ylist[j] = float(ylist[j])

            t0 = float(entrada[tam-5])
            h = float(entrada[tam - 4])
            n = int(entrada[tam - 3])
            f = parse_expr(entrada[tam - 2])
            ordem = int(entrada[tam -1])

            Adam_Bashforth(ylist , t0 , h , n , f , ordem )
            print()

        elif(metodo == 'adam_bashforth_by_euler'):

            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler(y0 , t0 , h , f , ordem)
            Adam_Bashforth(ylist , t0 , h , n , f , ordem , ' por Euler')
            print() 


        elif(metodo == 'adam_bashforth_by_euler_inverso'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler_Inverso(y0 , t0 , h , f , ordem)
            Adam_Bashforth(ylist , t0 , h , n , f , ordem , ' por Euler Inverso')
            print()

        elif(metodo == 'adam_bashforth_by_euler_aprimorado'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler_aprimorado(y0 , t0 , h , f , ordem)
            Adam_Bashforth(ylist , t0 , h , n , f , ordem , ' por Euler Aprimorado')
            print()

        elif(metodo == 'adam_bashforth_by_runge_kutta'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaRunge_Kutta(y0 , t0 , h , f , ordem)
            Adam_Bashforth(ylist , t0 , h , n , f , ordem , ' por Runge-Kutta ( ordem = {} )'.format(ordem))
            print()    
        

        elif(metodo == 'adam_multon'): #ta pegando para ordem = numero de pontos anteriores

            ylist = []
            tam = len(entrada)

            
            #funcao que pega os valores anteriores 
            for j in range(int(entrada[tam -1])) :
                ylist.append(entrada[j+1])
                ylist[j] = float(ylist[j])

            t0 = float(entrada[tam-5])
            h = float(entrada[tam - 4])
            n = int(entrada[tam - 3])
            f = parse_expr(entrada[tam - 2])
            ordem = int(entrada[tam -1])

            Adam_Multon(ylist , t0 , h , n , f , ordem )
            print()

        elif(metodo == 'adam_multon_by_euler'):

            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler(y0 , t0 , h , f , ordem)
            Adam_Multon(ylist , t0 , h , n , f , ordem , ' por Euler')
            print()    

        elif(metodo == 'adam_multon_by_euler_inverso'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler_Inverso(y0 , t0 , h , f , ordem)
            Adam_Multon(ylist , t0 , h , n , f , ordem , ' por Euler Inverso')
            print()

        elif(metodo == 'adam_multon_by_euler_aprimorado'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler_aprimorado(y0 , t0 , h , f , ordem)
            Adam_Multon(ylist , t0 , h , n , f , ordem , ' por Euler Aprimorado')
            print()

        elif(metodo == 'adam_multon_by_runge_kutta'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaRunge_Kutta(y0 , t0 , h , f , ordem)
            Adam_Multon(ylist , t0 , h , n , f , ordem , ' por Runge-Kutta ( ordem = {} )'.format(ordem))
            print()

        elif(metodo == 'formula_inversa'):#ta pegando para ordem = numero de pontos anteriores


            ylist = []
            tam = len(entrada)

            
            #funcao que pega os valores anteriores 
            for j in range(int(entrada[tam -1])) :
                ylist.append(entrada[j+1])
                ylist[j] = float(ylist[j])

            t0 = float(entrada[tam-5])
            h = float(entrada[tam - 4])
            n = int(entrada[tam - 3])
            f = parse_expr(entrada[tam - 2])
            ordem = int(entrada[tam -1])

            Formula_inversa(ylist , t0 , h , n , f , ordem  )
            print() 

        elif(metodo == 'formula_inversa_by_euler'):

            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler(y0 , t0 , h , f , ordem)
            Formula_inversa(ylist , t0 , h , n , f , ordem , ' por Euler')
            print() 

        elif(metodo == 'formula_inversa_by_euler_inverso'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler_Inverso(y0 , t0 , h , f , ordem)
            Formula_inversa(ylist , t0 , h , n , f , ordem , ' por Euler Inverso')
            print()  

        elif(metodo == 'formula_inversa_by_euler_aprimorado'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaEuler_aprimorado(y0 , t0 , h , f , ordem)
            Formula_inversa(ylist , t0 , h , n , f , ordem , ' por Euler Aprimorado')
            print()

        elif(metodo == 'formula_inversa_by_runge_kutta'):
        
            y0 = float(entrada[1])
            t0 = float(entrada[2])
            h = float(entrada[3])
            n = int(entrada[4])
            f = parse_expr(entrada[5])
            ordem = int(entrada[6])

            ylist = ListaRunge_Kutta(y0 , t0 , h , f , ordem)
            Formula_inversa(ylist , t0 , h , n , f , ordem , ' por Runge-Kutta ( ordem = {} )'.format(ordem))
            print()

        
 
