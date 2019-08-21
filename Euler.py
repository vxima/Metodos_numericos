#Metodo de Euler
print("Metodo de Euler")
def euler( y0 , t0 , h , n ):
    t=[t0]
    y=[y0]
    for j in range(n):
        f = float(1-t[j]+4*y[j])
        y.append(float(y[j] + h*f))
        t.append(float(t[j] + h))
    return ( y , t )
print(euler( 1.0 , 0.0 , 0.05 , 4))
